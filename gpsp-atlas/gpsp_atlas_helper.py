import numpy as np
from astropy.io import fits
from astropy.table import QTable
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator, interp1d

class GPSPAtlasHelper:
    """
    Galactic Photon Survival Probability (GPSP) Atlas Helper.

    Provides tools for loading, interpolating, and analyzing 4D GPSP FITS datasets.

    Notes
    -----
    FITS axis convention assumed: (E, D, B, L)
    NumPy axis convention assumed: (L, B, D, E)

    Parameters
    ----------
    file_path : str
        Path to the GPSP Atlas FITS file.
    """

    def __init__(self, file_path):
        self.file_path = file_path
        self.hdul = fits.open(self.file_path, memmap=True, lazy_load_hdus=False)
        self.hdr = self.hdul[0].header
        self.data = self.hdul[0].data

        # Axis metadata
        self.n_l = self.hdr['NAXIS4']
        self.crval_l = self.hdr['CRVAL4']
        self.dl = self.hdr['CDELT4']

        self.n_b = self.hdr['NAXIS3']
        self.crval_b = self.hdr['CRVAL3']
        self.db = self.hdr['CDELT3']

        self.n_d = self.hdr['NAXIS2']
        self.crval_d = self.hdr['CRVAL2']
        self.dd = self.hdr['CDELT2']

        # Physical axes
        self.l_grid = self.crval_l + np.arange(self.n_l) * self.dl
        self.b_grid = self.crval_b + np.arange(self.n_b) * self.db
        self.d_axis = self.crval_d + np.arange(self.n_d) * self.dd

        if 'ENERGIES' in self.hdul:
            self.e_axis = self.hdul['ENERGIES'].data['ENERGY']
        else:
            self.e_axis = self.hdul[1].data['ENERGY']

        self.log_e = np.log10(self.e_axis)

        # 4D interpolator
        self.interpolator = RegularGridInterpolator(
            (self.l_grid, self.b_grid, self.d_axis, self.log_e),
            self.data,
            bounds_error=False,
            fill_value=None
        )

    def close(self):
        """Close the underlying FITS file handler."""
        self.hdul.close()

    def _get_indices(self, l, b, d=None):
        """
        Convert physical coordinates to grid indices.

        Parameters
        ----------
        l : float
            Galactic longitude in degrees.
        b : float
            Galactic latitude in degrees.
        d : float, optional
            Distance in kpc.

        Returns
        -------
        tuple
            Tuple of nearest grid indices (l_idx, b_idx, [d_idx]).
        """
        l_idx = np.clip(int(round((l - self.crval_l) / self.dl)), 0, self.n_l - 1)
        b_idx = np.clip(int(round((b - self.crval_b) / self.db)), 0, self.n_b - 1)
        if d is not None:
            d_idx = np.clip(int(round((d - self.crval_d) / self.dd)), 0, self.n_d - 1)
            return l_idx, b_idx, d_idx
        return l_idx, b_idx

    def _weighted_percentile(self, values, weights, percentile):
        """
        Calculate a weighted percentile of a 1D numpy array.

        Parameters
        ----------
        values : array_like
            Data values.
        weights : array_like
            Weights for the data.
        percentile : float
            Percentile to compute (0-100).

        Returns
        -------
        float
            Weighted percentile value.
        """
        sorter = np.argsort(values)
        values = values[sorter]
        weights = weights[sorter]
        cdf = np.cumsum(weights)
        cdf /= cdf[-1]
        return np.interp(percentile / 100.0, cdf, values)

    def _get_survival_profile_interpolated(self, l, b, d):
        """
        Extract interpolated survival profile across all energies for a point.

        Parameters
        ----------
        l, b, d : float
            Galactic coordinates and distance.

        Returns
        -------
        ndarray
            Interpolated survival probability array.
        """
        pts = np.column_stack([
            np.full_like(self.log_e, l),
            np.full_like(self.log_e, b),
            np.full_like(self.log_e, d),
            self.log_e
        ])
        return self.interpolator(pts)

    def get_survival_vs_energy(self, l, b, target_distance_kpc):
        """
        Extract survival probability as a function of energy for a fixed distance.

        Parameters
        ----------
        l : float
            Galactic longitude in degrees.
        b : float
            Galactic latitude in degrees.
        target_distance_kpc : float
            Target distance in kpc.

        Returns
        -------
        tuple
            (Energy axis array in eV, Interpolated survival probabilities array).
        """
        l_idx, b_idx = self._get_indices(l, b)
        
        # Locate indices for the distance interpolation
        ih = np.searchsorted(self.d_axis, target_distance_kpc)
        il = np.clip(ih - 1, 0, len(self.d_axis) - 2)
        ih = il + 1

        # Extract the energy-dependent arrays at the bounding distance slices
        p_low = self.data[l_idx, b_idx, il, :]
        p_high = self.data[l_idx, b_idx, ih, :]
        
        # Calculate linear interpolation weight
        weight = (target_distance_kpc - self.d_axis[il]) / (self.d_axis[ih] - self.d_axis[il])
        interp_p = p_low + weight * (p_high - p_low)
        
        return self.e_axis.copy(), interp_p

    def get_survival_vs_distance(self, l, b, target_energy_ev):
        """
        Extract survival probability as a function of distance for a fixed energy.

        Parameters
        ----------
        l : float
            Galactic longitude.
        b : float
            Galactic latitude.
        target_energy_ev : float
            Target energy in eV.

        Returns
        -------
        tuple
            (Distance axis array, Interpolated survival probabilities array).
        """
        l_idx, b_idx = self._get_indices(l, b)
        log_target = np.log10(target_energy_ev)
        ih = np.searchsorted(self.e_axis, target_energy_ev)
        il = np.clip(ih - 1, 0, len(self.e_axis) - 2)
        ih = il + 1

        p_low = self.data[l_idx, b_idx, :, il]
        p_high = self.data[l_idx, b_idx, :, ih]
        weight = (log_target - self.log_e[il]) / (self.log_e[ih] - self.log_e[il])
        interp_p = p_low + weight * (p_high - p_low)
        return self.d_axis.copy(), interp_p

    def get_horizon_distances(self, l, b, energies_ev, target_p):
        """
        Calculate horizon distance where survival probability drops below a threshold.

        Parameters
        ----------
        l : float
            Galactic longitude.
        b : float
            Galactic latitude.
        energies_ev : array_like
            Array of energies in eV.
        target_p : float
            Target survival probability threshold.

        Returns
        -------
        ndarray
            Array of horizon distances matching the energy array.
        """
        energies_ev = np.atleast_1d(energies_ev)
        horizons = []
        for E in energies_ev:
            dist, p_vals = self.get_survival_vs_distance(l, b, E)
            if np.min(p_vals) > target_p:
                horizons.append(np.nan)
            else:
                d_hor = np.interp(target_p, p_vals[::-1], dist[::-1])
                horizons.append(d_hor)
        return np.array(horizons)

    def get_bird_view_map(self, target_energy_ev, b_lower=-0.5, b_upper=0.5, method='mean', percentile=50):
        """
        Extract a 2D map slice (longitude vs distance) collapsed along latitude.

        Parameters
        ----------
        target_energy_ev : float
            Target energy in eV.
        b_lower : float, optional
            Lower latitude bound, by default -0.5.
        b_upper : float, optional
            Upper latitude bound, by default 0.5.
        method : {'mean', 'median', 'min', 'max', 'percentile'}, optional
            Statistical method to collapse the latitude axis, by default 'mean'.
        percentile : float, optional
            Percentile to use if method is 'percentile', by default 50.

        Returns
        -------
        tuple
            (Longitude grid, Distance grid, 2D collapsed data map).
        """
        b1 = self._get_indices(0, b_lower)[1]
        b2 = self._get_indices(0, b_upper)[1]
        if b1 > b2:
            b1, b2 = b2, b1

        log_t = np.log10(target_energy_ev)
        ih = np.searchsorted(self.e_axis, target_energy_ev)
        il = np.clip(ih - 1, 0, len(self.e_axis) - 2)
        ih = il + 1
        w = (log_t - self.log_e[il]) / (self.log_e[ih] - self.log_e[il])

        slice_low = self.data[:, b1:b2 + 1, :, il]
        slice_high = self.data[:, b1:b2 + 1, :, ih]

        stat_funcs = {
            'mean': lambda x: np.mean(x, axis=1),
            'median': lambda x: np.median(x, axis=1),
            'min': lambda x: np.min(x, axis=1),
            'max': lambda x: np.max(x, axis=1),
            'percentile': lambda x: np.percentile(x, percentile, axis=1)
        }

        m_l = stat_funcs[method](slice_low)
        m_h = stat_funcs[method](slice_high)
        return self.l_grid, self.d_axis, m_l + w * (m_h - m_l)

    def get_galactic_plane_map(
        self,
        target_energy_ev,
        target_distance_kpc,
        b_lower=-5.0,
        b_upper=5.0,
        chunk_size=500,
        output_fits=None,
        overwrite=True,
        dtype=np.float32
    ):
        """
        Extract a Galactic longitude-latitude survival probability map
        at a fixed energy and distance.

        Parameters
        ----------
        target_energy_ev : float
            Target energy in eV.
        target_distance_kpc : float
            Target distance in kpc.
        b_lower : float, optional
            Lower Galactic latitude bound in degrees.
        b_upper : float, optional
            Upper Galactic latitude bound in degrees.
        chunk_size : int, optional
            Number of longitude pixels processed per chunk.
        output_fits : str, optional
            If provided, save the extracted map to a FITS file.
        overwrite : bool, optional
            Overwrite existing FITS file.
        dtype : numpy dtype, optional
            Internal storage dtype.

        Returns
        -------
        tuple
            (map_2d, l_grid, b_grid)

            map_2d shape = (n_b, n_l)
        """

        # --------------------------------------------
        # Latitude limits
        # --------------------------------------------
        _, b1 = self._get_indices(0.0, b_lower)
        _, b2 = self._get_indices(0.0, b_upper)

        if b1 > b2:
            b1, b2 = b2, b1

        # --------------------------------------------
        # Energy interpolation parameters
        # --------------------------------------------
        log_target = np.log10(target_energy_ev)

        ih_e = np.searchsorted(self.e_axis, target_energy_ev)
        il_e = np.clip(ih_e - 1, 0, len(self.e_axis) - 2)
        ih_e = il_e + 1

        w_e = (
            (log_target - self.log_e[il_e])
            / (self.log_e[ih_e] - self.log_e[il_e])
        )

        # --------------------------------------------
        # Distance interpolation parameters
        # --------------------------------------------
        float_d = (target_distance_kpc - self.crval_d) / self.dd

        il_d = np.clip(
            int(np.floor(float_d)),
            0,
            self.n_d - 2
        )

        ih_d = il_d + 1
        w_d = float_d - il_d

        # --------------------------------------------
        # Allocate output map
        # --------------------------------------------
        n_b_slice = b2 - b1 + 1

        map_raw = np.empty(
            (self.n_l, n_b_slice),
            dtype=dtype
        )

        n_chunks = (
            self.n_l + chunk_size - 1
        ) // chunk_size

        # --------------------------------------------
        # Chunked extraction
        # --------------------------------------------
        for chunk_idx in range(n_chunks):

            progress = (chunk_idx + 1) / n_chunks * 100
            print(
                f"\rProcessing: {progress:5.1f}% "
                f"(chunk {chunk_idx+1}/{n_chunks})",
                end="",
                flush=True
            )

            l_start = chunk_idx * chunk_size
            l_end = min(
                (chunk_idx + 1) * chunk_size,
                self.n_l
            )

            sub_cube = np.asarray(
                self.data[
                    l_start:l_end,
                    b1:b2 + 1,
                    il_d:ih_d + 1,
                    il_e:ih_e + 1
                ],
                dtype=dtype
            )

            s00 = sub_cube[:, :, 0, 0]
            s01 = sub_cube[:, :, 0, 1]
            s10 = sub_cube[:, :, 1, 0]
            s11 = sub_cube[:, :, 1, 1]

            low_d = s00 + w_e * (s01 - s00)
            high_d = s10 + w_e * (s11 - s10)

            map_raw[l_start:l_end, :] = (
                low_d + w_d * (high_d - low_d)
            )

        # Loop finished
        print()

        # rows = latitude, cols = longitude
        map_2d = map_raw.T

        l_grid = self.l_grid.copy()
        b_grid = self.b_grid[b1:b2 + 1].copy()

        # --------------------------------------------
        # Optional FITS output
        # --------------------------------------------
        if output_fits is not None:

            hdu = fits.PrimaryHDU(map_2d)

            hdu.header['CRVAL1'] = float(self.crval_l)
            hdu.header['CDELT1'] = float(self.dl)
            hdu.header['CRPIX1'] = 1.0
            hdu.header['CTYPE1'] = 'GLON-CAR'

            hdu.header['CRVAL2'] = float(b_grid[0])
            hdu.header['CDELT2'] = float(self.db)
            hdu.header['CRPIX2'] = 1.0
            hdu.header['CTYPE2'] = 'GLAT-CAR'

            hdu.header['ENERGYEV'] = float(target_energy_ev)
            hdu.header['DISTKPC'] = float(target_distance_kpc)
            hdu.header['BUNIT'] = 'Probability'

            hdu.writeto(
                output_fits,
                overwrite=overwrite
            )

        return map_2d, l_grid, b_grid

    def get_2d_phase_space(self, l, b):
        """
        Extract the 2D energy-distance phase space for a given line of sight.

        Parameters
        ----------
        l : float
            Galactic longitude.
        b : float
            Galactic latitude.

        Returns
        -------
        tuple
            (Meshgrid array [Energy, Distance], 2D phase space data array).
        """
        l_idx, b_idx = self._get_indices(l, b)
        return np.meshgrid(self.e_axis, self.d_axis), self.data[l_idx, b_idx, :, :].copy()

    def _generate_region_points(self, l, b, shape='point', r_in=0.0, r_out=0.0, phi_start=0.0, phi_stop=360.0, 
                                width_l=None, width_b=None, n_points=20, radial_sampling='equal_area', center_density_power=1.0):
        """Generate a set of (l, b) coordinates based on geometry configurations."""
        points = []
        if shape == 'point':
            return np.array([(l, b)])
        elif shape == 'box':
            if width_l is None or width_b is None:
                raise ValueError("width_l and width_b required.")
            n_side = int(np.ceil(np.sqrt(n_points)))
            l_vals = np.linspace(l - width_l / 2, l + width_l / 2, n_side)
            b_vals = np.linspace(b - width_b / 2, b + width_b / 2, n_side)
            for ll in l_vals:
                for bb in b_vals:
                    points.append((ll, bb))
            return np.array(points)
        elif shape == 'extended':
            if r_out <= 0:
                raise ValueError("r_out must be > 0.")
            n_theta = max(8, int(np.sqrt(n_points) * 2))
            theta_vals = np.linspace(np.deg2rad(phi_start), np.deg2rad(phi_stop), n_theta, endpoint=False)
            n_r = max(2, int(np.ceil(n_points / n_theta)))

            if radial_sampling == 'equal_area':
                uvals = np.linspace(0, 1, n_r)
                r_vals = np.sqrt(uvals * (r_out**2 - r_in**2) + r_in**2)
            elif radial_sampling == 'uniform_radius':
                r_vals = np.linspace(r_in, r_out, n_r)
            elif radial_sampling == 'center_enhanced':
                uvals = np.linspace(0, 1, n_r)
                r_vals = r_in + (r_out - r_in) * (uvals**center_density_power)
            else:
                raise ValueError("Invalid radial_sampling.")

            if r_in == 0:
                points.append((l, b))
            for r in r_vals:
                if r == 0: continue
                for theta in theta_vals:
                    points.append((l + r * np.cos(theta), b + r * np.sin(theta)))
            return np.array(points)
        else:
            raise ValueError(f"Unknown shape: {shape}")

    def _compute_weights(self, points, center_l, center_b, weighting='flat', gaussian_sigma=None):
        """Compute spatial weights for a set of generated coordinates."""
        if weighting == 'flat':
            weights = np.ones(len(points))
        elif weighting == 'gaussian':
            if gaussian_sigma is None:
                raise ValueError("gaussian_sigma required.")
            dl = points[:, 0] - center_l
            db = points[:, 1] - center_b
            weights = np.exp(-0.5 * (dl**2 + db**2) / gaussian_sigma**2)
        else:
            raise ValueError(f"Unknown weighting: {weighting}")
        weights /= np.sum(weights)
        return weights

    def get_region_survival_profile(self, l, b, d, d_err=None, shape='point', r_in=0.0, r_out=0.0, 
                                    phi_start=0.0, phi_stop=360.0, width_l=None, width_b=None, n_points=20, 
                                    weighting='flat', gaussian_sigma=None, average_mode='percentile', 
                                    percentile=68, include_distance_uncertainty=True, radial_sampling='equal_area', 
                                    center_density_power=1.0, distance_sampling='mc', n_distance_samples=200, 
                                    distance_sigma_clip=3.0, random_seed=42, return_individual=True):
        """
        Calculate an averaged survival probability profile across an extended region and/or distance uncertainty.

        Parameters
        ----------
        l, b, d : float
            Center coordinates and distance.
        d_err : float, optional
            Distance error constraint.
        shape : {'point', 'box', 'extended'}, optional
            Spatial geometry of the target, by default 'point'.
        r_in, r_out : float, optional
            Inner and outer radii for extended shapes.
        phi_start, phi_stop : float, optional
            Angular constraints for extended shapes.
        width_l, width_b : float, optional
            Box dimensions for box shapes.
        n_points : int, optional
            Number of spatial sampling points, by default 20.
        weighting : {'flat', 'gaussian'}, optional
            Weighting schema, by default 'flat'.
        gaussian_sigma : float, optional
            Sigma for gaussian weighting.
        average_mode : {'percentile', 'minmax'}, optional
            Method to establish confidence intervals, by default 'percentile'.
        percentile : float, optional
            Confidence interval percentage, by default 68.
        include_distance_uncertainty : bool, optional
            Integrate over distance uncertainty, by default True.
        radial_sampling : str, optional
            Radial spacing distribution algorithm, by default 'equal_area'.
        center_density_power : float, optional
            Center density scale, by default 1.0.
        distance_sampling : {'mc', 'grid'}, optional
            Method for distance sampling, by default 'mc'.
        n_distance_samples : int, optional
            Number of MC distance samplings, by default 200.
        distance_sigma_clip : float, optional
            Sigma cutoff for distance integration, by default 3.0.
        random_seed : int, optional
            RNG seed for reproducibility, by default 42.
        return_individual : bool, optional
            Return sub-profiles in final dictionary, by default True.

        Returns
        -------
        dict
            Contains energies_ev, average, low/high bounds, points, weights, distance_samples.
        """
        points = self._generate_region_points(
            l=l, b=b, shape=shape, r_in=r_in, r_out=r_out, phi_start=phi_start, phi_stop=phi_stop, 
            width_l=width_l, width_b=width_b, n_points=n_points, radial_sampling=radial_sampling, 
            center_density_power=center_density_power
        )
        weights = self._compute_weights(points, center_l=l, center_b=b, weighting=weighting, gaussian_sigma=gaussian_sigma)
        rng = np.random.default_rng(random_seed)

        if include_distance_uncertainty and d_err is not None and d_err > 0:
            if distance_sampling == 'grid':
                mask = (self.d_axis >= d - d_err) & (self.d_axis <= d + d_err)
                d_vals = self.d_axis[mask]
                if len(d_vals) == 0:
                    d_vals = np.array([d])
                d_weights = np.ones(len(d_vals)) / len(d_vals)
            elif distance_sampling == 'mc':
                d_vals = rng.normal(loc=d, scale=d_err, size=n_distance_samples)
                mask = np.abs(d_vals - d) <= distance_sigma_clip * d_err
                d_vals = np.clip(d_vals[mask], np.min(self.d_axis), np.max(self.d_axis))
                d_weights = np.ones(len(d_vals)) / len(d_vals)
            else:
                raise ValueError("distance_sampling must be 'grid' or 'mc'")
        else:
            d_vals = np.array([d])
            d_weights = np.array([1.0])

        all_profiles, all_weights = [], []
        for d_sample, dw in zip(d_vals, d_weights):
            for (ll, bb), w in zip(points, weights):
                prof = self._get_survival_profile_interpolated(ll, bb, d_sample)
                all_profiles.append(prof)
                all_weights.append(w * dw)

        all_profiles = np.array(all_profiles)
        all_weights = np.array(all_weights)
        average_profile = np.average(all_profiles, axis=0, weights=all_weights)

        if average_mode == 'percentile':
            low_q = (100 - percentile) / 2
            high_q = 100 - low_q
            low = np.array([self._weighted_percentile(all_profiles[:, i], all_weights, low_q) for i in range(all_profiles.shape[1])])
            high = np.array([self._weighted_percentile(all_profiles[:, i], all_weights, high_q) for i in range(all_profiles.shape[1])])
        elif average_mode == 'minmax':
            low = np.min(all_profiles, axis=0)
            high = np.max(all_profiles, axis=0)
        else:
            raise ValueError("average_mode must be 'percentile' or 'minmax'")

        result = {
            "energies_ev": self.e_axis.copy(),
            "average": average_profile,
            "low": low,
            "high": high,
            "points": points,
            "weights": weights,
            "distance_samples": d_vals,
            "distance_weights": d_weights
        }
        if return_individual:
            result["profiles"] = all_profiles
            result["profile_weights"] = all_weights

        return result

    def process_spectrum_with_gpsp(self, gpsp_result, mode='apply', ecsv_file=None, 
                                   energy=None, flux=None, flux_err=None, flux_errn=None, flux_errp=None, 
                                   energy_column='e_ref', flux_column='dnde', flux_err_column='dnde_err', 
                                   flux_errn_column='dnde_errn', flux_errp_column='dnde_errp', 
                                   calc_errors=True, error_comb_mode='linear', min_survival_prob=1e-6, 
                                   return_table=False, output_ecsv=None):
        """
        Apply or correct spectral fluxes using a computed GPSP profile.

        Parameters
        ----------
        gpsp_result : dict
            Result dictionary returned by get_region_survival_profile.
        mode : {'apply', 'correct'}, optional
            Direction of transformation, by default 'apply'.
        ecsv_file : str, optional
            Path to input ECSV table.
        energy, flux : array_like, optional
            Direct array inputs if not using ecsv_file.
        flux_err, flux_errn, flux_errp : array_like, optional
            Symmetric or asymmetric error arrays.
        energy_column, flux_column : str, optional
            Target column names for ECSV file.
        flux_err_column, flux_errn_column, flux_errp_column : str, optional
            Target error column names for ECSV file.
        calc_errors : bool, optional
            Propagate uncertainties, by default True.
        error_comb_mode : {'linear', 'quadrature'}, optional
            Error propagation method, by default 'linear'.
        min_survival_prob : float, optional
            Floor value for survival probability division, by default 1e-6.
        return_table : bool, optional
            Compile output into an astropy QTable, by default False.
        output_ecsv : str, optional
            Path to dump out final processed ECSV file.

        Returns
        -------
        dict
            Dictionary containing processed energy, fluxes, bounds, and optionally the constructed QTable.
        """
        if ecsv_file is not None:
            t = QTable.read(ecsv_file, format='ascii.ecsv')
            energy = t[energy_column]
            flux = t[flux_column]
            if flux_err_column in t.colnames: flux_err = t[flux_err_column]
            if flux_errn_column in t.colnames: flux_errn = t[flux_errn_column]
            if flux_errp_column in t.colnames: flux_errp = t[flux_errp_column]
        else:
            if energy is None or flux is None:
                raise ValueError("Provide either ecsv_file or (energy, flux).")

        energy_ev = energy.to(u.eV).value
        gpsp_e = gpsp_result["energies_ev"]
        
        interp_med = interp1d(np.log10(gpsp_e), gpsp_result["average"], bounds_error=False, fill_value='extrapolate')
        interp_low = interp1d(np.log10(gpsp_e), gpsp_result["low"], bounds_error=False, fill_value='extrapolate')
        interp_high = interp1d(np.log10(gpsp_e), gpsp_result["high"], bounds_error=False, fill_value='extrapolate')

        p_med = np.clip(interp_med(np.log10(energy_ev)), min_survival_prob, None)
        p_low = np.clip(interp_low(np.log10(energy_ev)), min_survival_prob, None)
        p_high = np.clip(interp_high(np.log10(energy_ev)), min_survival_prob, None)

        if mode == 'apply':
            flux_out = flux * p_med
            flux_sys_low = flux * p_low
            flux_sys_high = flux * p_high
        elif mode == 'correct':
            flux_out = flux / p_med
            flux_sys_low = flux / p_high
            flux_sys_high = flux / p_low
        else:
            raise ValueError("mode must be 'apply' or 'correct'")

        result = {
            "energy": energy,
            "flux": flux_out,
            "flux_sys_low": flux_sys_low,
            "flux_sys_high": flux_sys_high,
            "survival_probability": p_med,
            "survival_low": p_low,
            "survival_high": p_high
        }

        if calc_errors:
            scale = p_med if mode == 'apply' else 1.0 / p_med

            # Symmetric statistical errors
            if flux_err is not None:
                flux_err_stat = flux_err * scale
                delta_sys_low = flux_out - flux_sys_low
                delta_sys_high = flux_sys_high - flux_out

                if error_comb_mode == 'linear':
                    flux_total_low = flux_out - (delta_sys_low + flux_err_stat)
                    flux_total_high = flux_out + (delta_sys_high + flux_err_stat)
                elif error_comb_mode == 'quadrature':
                    flux_total_low = flux_out - np.sqrt(delta_sys_low**2 + flux_err_stat**2)
                    flux_total_high = flux_out + np.sqrt(delta_sys_high**2 + flux_err_stat**2)
                else:
                    raise ValueError("error_comb_mode must be 'linear' or 'quadrature'")
                
                flux_total_low = flux_total_low.clip(min=0 * flux_total_low.unit)
                result["flux_err"] = flux_err_stat
                result["flux_total_low"] = flux_total_low
                result["flux_total_high"] = flux_total_high

            # Asymmetric statistical errors
            elif flux_errn is not None and flux_errp is not None:
                flux_errn_stat = flux_errn * scale
                flux_errp_stat = flux_errp * scale
                delta_sys_low = flux_out - flux_sys_low
                delta_sys_high = flux_sys_high - flux_out

                if error_comb_mode == 'linear':
                    total_err_low = delta_sys_low + flux_errn_stat
                    total_err_high = delta_sys_high + flux_errp_stat
                elif error_comb_mode == 'quadrature':
                    total_err_low = np.sqrt(delta_sys_low**2 + flux_errn_stat**2)
                    total_err_high = np.sqrt(delta_sys_high**2 + flux_errp_stat**2)
                else:
                    raise ValueError("error_comb_mode must be 'linear' or 'quadrature'")

                flux_total_low = (flux_out - total_err_low).clip(min=0 * flux_out.unit)
                result["flux_errn"] = flux_errn_stat
                result["flux_errp"] = flux_errp_stat
                result["flux_total_low"] = flux_total_low
                result["flux_total_high"] = flux_out + total_err_high

        if return_table:
            table = QTable()
            table["energy"] = energy
            table["flux"] = flux_out
            if "flux_err" in result: table["flux_err"] = result["flux_err"]
            if "flux_errn" in result: table["flux_errn"] = result["flux_errn"]
            if "flux_errp" in result: table["flux_errp"] = result["flux_errp"]
            
            table["flux_sys_low"] = flux_sys_low
            table["flux_sys_high"] = flux_sys_high
            if "flux_total_low" in result: table["flux_total_low"] = result["flux_total_low"]
            if "flux_total_high" in result: table["flux_total_high"] = result["flux_total_high"]
            
            table["survival_probability"] = p_med
            table["survival_low"] = p_low
            table["survival_high"] = p_high
            result["table"] = table

            if output_ecsv is not None:
                table.write(output_ecsv, format='ascii.ecsv', overwrite=True)

        return result
