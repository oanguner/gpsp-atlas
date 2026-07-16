import numpy as np
from astropy.io import fits
from astropy.table import QTable
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator, interp1d

class GPSPLIVAtlasHelper:
    """
    Galactic Photon Survival Probability LIV (GPSP-LIV) Atlas Helper.
    
    A utility class to load, interpolate, and process a 5D FITS data cube containing 
    gamma-ray survival probabilities to evaluate Lorentz Invariance Violation (LIV).

    FITS axis order: (Lambda, E, D, B, L)
    NumPy array order: (L, B, D, E, Lambda)
    """

    def __init__(self, file_path):
        """
        Initialize the helper and load the 5D FITS data cube.

        Args:
            file_path (str): Path to the GPSP-LIV FITS file.
        """
        self.file_path = file_path
        self.hdul = fits.open(self.file_path, memmap=True, lazy_load_hdus=False)
        self.hdr, self.data = self.hdul[0].header, self.hdul[0].data

        # Axis metadata
        self.n_l, self.crval_l, self.dl = self.hdr['NAXIS5'], self.hdr['CRVAL5'], self.hdr['CDELT5']
        self.n_b, self.crval_b, self.db = self.hdr['NAXIS4'], self.hdr['CRVAL4'], self.hdr['CDELT4']
        self.n_d, self.crval_d, self.dd = self.hdr['NAXIS3'], self.hdr['CRVAL3'], self.hdr['CDELT3']
        self.n_e = self.hdr['NAXIS2']
        self.n_lam, self.crval_lam, self.dlam = self.hdr['NAXIS1'], self.hdr['CRVAL1'], self.hdr['CDELT1']

        # Physical axes
        self.l_grid = self.crval_l + np.arange(self.n_l) * self.dl
        self.b_grid = self.crval_b + np.arange(self.n_b) * self.db
        self.d_axis = self.crval_d + np.arange(self.n_d) * self.dd
        self.lam_axis = self.crval_lam + np.arange(self.n_lam) * self.dlam

        self.e_axis = self.hdul['ENERGIES'].data['ENERGY'] if 'ENERGIES' in self.hdul else self.hdul[1].data['ENERGY']
        self.log_e = np.log10(self.e_axis)

        # 5D interpolator
        self.interpolator = RegularGridInterpolator(
            (self.l_grid, self.b_grid, self.d_axis, self.log_e, self.lam_axis),
            self.data, bounds_error=False, fill_value=None
        )

    def close(self):
        """Close the underlying FITS file handler."""
        self.hdul.close()

    def _get_indices(self, l, b, d=None):
        """Convert physical coordinates to grid indices."""
        l_idx = np.clip(int(round((l - self.crval_l) / self.dl)), 0, self.n_l - 1)
        b_idx = np.clip(int(round((b - self.crval_b) / self.db)), 0, self.n_b - 1)
        if d is not None:
            d_idx = np.clip(int(round((d - self.crval_d) / self.dd)), 0, self.n_d - 1)
            return l_idx, b_idx, d_idx
        return l_idx, b_idx

    def _get_lambda_index(self, lambda_):
        """Convert LIV lambda scale to grid index."""
        return np.clip(int(round((lambda_ - self.crval_lam) / self.dlam)), 0, self.n_lam - 1)

    def _weighted_percentile(self, values, weights, percentile):
        """Compute a weighted percentile for an array of values."""
        sorter = np.argsort(values)
        values, weights = values[sorter], weights[sorter]
        cdf = np.cumsum(weights)
        cdf /= cdf[-1]
        return np.interp(percentile / 100.0, cdf, values)

    def _get_survival_profile_interpolated(self, l, b, d, lambda_):
        """Evaluate the 5D interpolator for a specific parameter combination."""
        pts = np.column_stack([
            np.full_like(self.log_e, l), np.full_like(self.log_e, b),
            np.full_like(self.log_e, d), self.log_e, np.full_like(self.log_e, lambda_)
        ])
        return self.interpolator(pts)

    def get_survival_vs_distance(self, l, b, target_energy_ev, lambda_):
        """
        Extract the survival probability as a function of distance.

        Args:
            l (float): Galactic longitude in degrees.
            b (float): Galactic latitude in degrees.
            target_energy_ev (float): Target photon energy in eV.
            lambda_ (float): LIV scale parameter.

        Returns:
            tuple: (distance_axis (ndarray), survival_probabilities (ndarray))
        """
        l_idx, b_idx = self._get_indices(l, b)
        lam_idx = self._get_lambda_index(lambda_)
        log_target = np.log10(target_energy_ev)

        ih = np.searchsorted(self.e_axis, target_energy_ev)
        il = np.clip(ih - 1, 0, len(self.e_axis) - 2)
        ih = il + 1

        p_low = self.data[l_idx, b_idx, :, il, lam_idx]
        p_high = self.data[l_idx, b_idx, :, ih, lam_idx]
        weight = (log_target - self.log_e[il]) / (self.log_e[ih] - self.log_e[il])
        
        return self.d_axis.copy(), p_low + weight * (p_high - p_low)

    def get_horizon_distances(self, l, b, energies_ev, lambda_, target_p):
        """
        Calculate the distance at which survival probability drops to a target value.

        Args:
            l (float): Galactic longitude in degrees.
            b (float): Galactic latitude in degrees.
            energies_ev (float or ndarray): Energy values in eV.
            lambda_ (float): LIV scale parameter.
            target_p (float): Target survival probability threshold (e.g., 0.367 for 1/e).

        Returns:
            numpy.ndarray: Horizon distances in kpc (NaN if threshold is never crossed).
        """
        energies_ev = np.atleast_1d(energies_ev)
        horizons = []
        for E in energies_ev:
            dist, p_vals = self.get_survival_vs_distance(l, b, E, lambda_)
            if np.min(p_vals) > target_p:
                horizons.append(np.nan)
            else:
                horizons.append(np.interp(target_p, p_vals[::-1], dist[::-1]))
        return np.array(horizons)

    def get_bird_view_map(self, target_energy_ev, lambda_, b_lower=-0.5, b_upper=0.5, method='mean', percentile=50):
        """
        Generate a 2D map (Longitude vs. Distance) of survival probabilities.

        Args:
            target_energy_ev (float): Energy to map in eV.
            lambda_ (float): LIV scale parameter.
            b_lower (float, optional): Lower Galactic latitude bound. Defaults to -0.5.
            b_upper (float, optional): Upper Galactic latitude bound. Defaults to 0.5.
            method (str, optional): Statistical method across latitude ('mean', 'median', 'min', 'max', 'percentile').
            percentile (float, optional): Percentile value if method='percentile'. Defaults to 50.

        Returns:
            tuple: (longitude_grid, distance_axis, 2D_map_array)
        """
        b1, b2 = sorted([self._get_indices(0, b_lower)[1], self._get_indices(0, b_upper)[1]])
        lam_idx = self._get_lambda_index(lambda_)
        log_t = np.log10(target_energy_ev)

        ih = np.searchsorted(self.e_axis, target_energy_ev)
        il = np.clip(ih - 1, 0, len(self.e_axis) - 2)
        ih = il + 1

        w = (log_t - self.log_e[il]) / (self.log_e[ih] - self.log_e[il])
        slice_low = self.data[:, b1:b2 + 1, :, il, lam_idx]
        slice_high = self.data[:, b1:b2 + 1, :, ih, lam_idx]

        stat_funcs = {
            'mean': lambda x: np.mean(x, axis=1),
            'median': lambda x: np.median(x, axis=1),
            'min': lambda x: np.min(x, axis=1),
            'max': lambda x: np.max(x, axis=1),
            'percentile': lambda x: np.percentile(x, percentile, axis=1)
        }

        m_l, m_h = stat_funcs[method](slice_low), stat_funcs[method](slice_high)
        return self.l_grid, self.d_axis, m_l + w * (m_h - m_l)

    def get_2d_phase_space(self, l, b, lambda_):
        """
        Retrieve the Energy vs. Distance phase space matrix for a specific direction.

        Args:
            l (float): Galactic longitude in degrees.
            b (float): Galactic latitude in degrees.
            lambda_ (float): LIV scale parameter.

        Returns:
            tuple: (meshgrid_grids, phase_space_array)
        """
        l_idx, b_idx = self._get_indices(l, b)
        lam_idx = self._get_lambda_index(lambda_)
        return np.meshgrid(self.e_axis, self.d_axis), self.data[l_idx, b_idx, :, :, lam_idx].copy()

    def _generate_region_points(self, l, b, shape='point', r_in=0.0, r_out=0.0, phi_start=0.0, phi_stop=360.0, 
                                width_l=None, width_b=None, n_points=20, radial_sampling='equal_area', center_density_power=1.0):
        """Generate spatial sampling points for extended regions or boxes."""
        if shape == 'point':
            return np.array([(l, b)])
        elif shape == 'box':
            if width_l is None or width_b is None: raise ValueError("width_l and width_b required.")
            n_side = int(np.ceil(np.sqrt(n_points)))
            l_vals = np.linspace(l - width_l / 2, l + width_l / 2, n_side)
            b_vals = np.linspace(b - width_b / 2, b + width_b / 2, n_side)
            return np.array([(ll, bb) for ll in l_vals for bb in b_vals])
        elif shape == 'extended':
            if r_out <= 0: raise ValueError("r_out must be > 0.")
            n_theta = max(8, int(np.sqrt(n_points) * 2))
            theta_vals = np.linspace(np.deg2rad(phi_start), np.deg2rad(phi_stop), n_theta, endpoint=False)
            n_r = max(2, int(np.ceil(n_points / n_theta)))
            uvals = np.linspace(0, 1, n_r)

            if radial_sampling == 'equal_area':
                r_vals = np.sqrt(uvals * (r_out**2 - r_in**2) + r_in**2)
            elif radial_sampling == 'uniform_radius':
                r_vals = np.linspace(r_in, r_out, n_r)
            elif radial_sampling == 'center_enhanced':
                r_vals = r_in + (r_out - r_in) * (uvals**center_density_power)
            else:
                raise ValueError("Invalid radial_sampling.")

            points = [(l, b)] if r_in == 0 else []
            for r in r_vals:
                if r > 0:
                    points.extend([(l + r * np.cos(t), b + r * np.sin(t)) for t in theta_vals])
            return np.array(points)
        raise ValueError(f"Unknown shape: {shape}")

    def _compute_weights(self, points, center_l, center_b, weighting='flat', gaussian_sigma=None):
        """Compute spatial weights for the sampled region points."""
        if weighting == 'flat':
            weights = np.ones(len(points))
        elif weighting == 'gaussian':
            if gaussian_sigma is None: raise ValueError("gaussian_sigma required.")
            weights = np.exp(-0.5 * ((points[:, 0] - center_l)**2 + (points[:, 1] - center_b)**2) / gaussian_sigma**2)
        else:
            raise ValueError(f"Unknown weighting: {weighting}")
        return weights / np.sum(weights)

    def get_region_survival_profile(self, l, b, d, lambda_, d_err=None, shape='point', r_in=0.0, r_out=0.0, 
                                    phi_start=0.0, phi_stop=360.0, width_l=None, width_b=None, n_points=20, 
                                    weighting='flat', gaussian_sigma=None, average_mode='percentile', percentile=68, 
                                    include_distance_uncertainty=True, radial_sampling='equal_area', 
                                    center_density_power=1.0, distance_sampling='mc', n_distance_samples=200, 
                                    distance_sigma_clip=3.0, random_seed=42, return_individual=True):
        """
        Calculate an averaged survival profile across a spatial region and distance uncertainty.

        Args:
            l (float): Center Galactic longitude in degrees.
            b (float): Center Galactic latitude in degrees.
            d (float): Source distance in kpc.
            lambda_ (float): LIV scale parameter.
            d_err (float, optional): Distance error in kpc.
            shape (str): Region shape ('point', 'box', 'extended').
            r_in (float): Inner radius for extended shapes.
            r_out (float): Outer radius for extended shapes.
            n_points (int): Approximate number of sampling points.
            average_mode (str): Method to aggregate profiles ('percentile' or 'minmax').
            include_distance_uncertainty (bool): If True, samples distance distribution.
            distance_sampling (str): Distance sampling method ('grid' or 'mc').
            return_individual (bool): If True, return full array of sampled profiles.

        Returns:
            dict: Dictionary containing averaged profiles, uncertainty bands, and sampling details.
        """
        points = self._generate_region_points(l, b, shape, r_in, r_out, phi_start, phi_stop, width_l, width_b, n_points, radial_sampling, center_density_power)
        weights = self._compute_weights(points, l, b, weighting, gaussian_sigma)
        rng = np.random.default_rng(random_seed)

        if include_distance_uncertainty and d_err and d_err > 0:
            if distance_sampling == 'grid':
                mask = (self.d_axis >= d - d_err) & (self.d_axis <= d + d_err)
                d_vals = self.d_axis[mask] if np.any(mask) else np.array([d])
            elif distance_sampling == 'mc':
                d_vals = rng.normal(loc=d, scale=d_err, size=n_distance_samples)
                d_vals = d_vals[np.abs(d_vals - d) <= distance_sigma_clip * d_err]
                d_vals = np.clip(d_vals, np.min(self.d_axis), np.max(self.d_axis))
            else:
                raise ValueError("distance_sampling must be 'grid' or 'mc'")
            d_weights = np.ones(len(d_vals)) / len(d_vals)
        else:
            d_vals, d_weights = np.array([d]), np.array([1.0])

        all_profiles, all_weights = [], []
        for d_sample, dw in zip(d_vals, d_weights):
            for (ll, bb), w in zip(points, weights):
                all_profiles.append(self._get_survival_profile_interpolated(ll, bb, d_sample, lambda_))
                all_weights.append(w * dw)

        all_profiles, all_weights = np.array(all_profiles), np.array(all_weights)
        average_profile = np.average(all_profiles, axis=0, weights=all_weights)

        if average_mode == 'percentile':
            low_q, high_q = (100 - percentile) / 2, 100 - (100 - percentile) / 2
            low = np.array([self._weighted_percentile(all_profiles[:, i], all_weights, low_q) for i in range(all_profiles.shape[1])])
            high = np.array([self._weighted_percentile(all_profiles[:, i], all_weights, high_q) for i in range(all_profiles.shape[1])])
        elif average_mode == 'minmax':
            low, high = np.min(all_profiles, axis=0), np.max(all_profiles, axis=0)
        else:
            raise ValueError("average_mode must be 'percentile' or 'minmax'")

        result = {
            "lambda": lambda_, "energies_ev": self.e_axis.copy(), "average": average_profile, 
            "low": low, "high": high, "points": points, "weights": weights, 
            "distance_samples": d_vals, "distance_weights": d_weights
        }
        if return_individual:
            result.update({"profiles": all_profiles, "profile_weights": all_weights})
        return result

    def process_spectrum_with_gpsp_liv(self, gpsp_result, mode='apply', spect_type='model', ecsv_file=None, energy=None, 
                                       flux=None, flux_err=None, flux_errn=None, flux_errp=None, energy_column='e_ref', 
                                       flux_column='dnde', flux_err_column='dnde_err', flux_errn_column='dnde_errn', 
                                       flux_errp_column='dnde_errp', calc_errors=True, error_comb_mode='linear', 
                                       min_survival_prob=1e-6, return_table=False, output_ecsv=None):
        """
        Apply or correct a gamma-ray spectrum using computed GPSP-LIV survival profiles.

        Args:
            gpsp_result (dict): Output dictionary from `get_region_survival_profile`.
            mode (str): Action to perform ('apply' to attenuate, 'correct' to de-absorb).
            spect_type (str): Origin type of spectrum ('model' or 'data').
            ecsv_file (str, optional): Path to input Astropy ECSV file.
            energy (astropy.units.Quantity, optional): Array of energies.
            flux (astropy.units.Quantity, optional): Array of flux values.
            calc_errors (bool): Whether to propagate uncertainties.
            error_comb_mode (str): How to combine errors ('linear' or 'quadrature').
            min_survival_prob (float): Threshold to prevent division by zero in correction.
            return_table (bool): If True, returns an astropy QTable in the result.
            output_ecsv (str, optional): Path to save the processed spectrum.

        Returns:
            dict: Processed spectrum data, survival fractions, and propagated errors.
        """
        
        # Read spectrum
        if ecsv_file:
            t = QTable.read(ecsv_file, format='ascii.ecsv')
            energy, flux = t[energy_column], t[flux_column]
            flux_err = t[flux_err_column] if flux_err_column in t.colnames else flux_err
            flux_errn = t[flux_errn_column] if flux_errn_column in t.colnames else flux_errn
            flux_errp = t[flux_errp_column] if flux_errp_column in t.colnames else flux_errp
        elif energy is None or flux is None:
            raise ValueError("Provide either ecsv_file or (energy, flux).")

        # GPSP interpolation
        energy_ev, gpsp_e = energy.to(u.eV).value, gpsp_result["energies_ev"]
        
        interp_med = interp1d(np.log10(gpsp_e), gpsp_result["average"], bounds_error=False, fill_value='extrapolate')
        interp_low = interp1d(np.log10(gpsp_e), gpsp_result["low"], bounds_error=False, fill_value='extrapolate')
        interp_high = interp1d(np.log10(gpsp_e), gpsp_result["high"], bounds_error=False, fill_value='extrapolate')

        p_med = np.clip(interp_med(np.log10(energy_ev)), min_survival_prob, None)
        p_low = np.clip(interp_low(np.log10(energy_ev)), min_survival_prob, None)
        p_high = np.clip(interp_high(np.log10(energy_ev)), min_survival_prob, None)

        # Apply / correct
        if mode == 'apply':
            flux_out, flux_sys_low, flux_sys_high = flux * p_med, flux * p_low, flux * p_high
        elif mode == 'correct':
            flux_out, flux_sys_low, flux_sys_high = flux / p_med, flux / p_high, flux / p_low
        else:
            raise ValueError("mode must be 'apply' or 'correct'")

        result = {
            "lambda": gpsp_result.get("lambda", None), "energy": energy, "flux": flux_out,
            "flux_sys_low": flux_sys_low, "flux_sys_high": flux_sys_high,
            "survival_probability": p_med, "survival_low": p_low, "survival_high": p_high
        }

        # Statistical errors
        if calc_errors:
            scale = p_med if mode == 'apply' else 1.0 / p_med
            d_sys_low, d_sys_high = flux_out - flux_sys_low, flux_sys_high - flux_out

            if flux_err is not None:
                err_stat = flux_err * scale
                if error_comb_mode == 'linear':
                    f_tot_low, f_tot_high = flux_out - (d_sys_low + err_stat), flux_out + (d_sys_high + err_stat)
                elif error_comb_mode == 'quadrature':
                    f_tot_low, f_tot_high = flux_out - np.sqrt(d_sys_low**2 + err_stat**2), flux_out + np.sqrt(d_sys_high**2 + err_stat**2)
                else:
                    raise ValueError("error_comb_mode must be 'linear' or 'quadrature'")
                
                result.update({"flux_err": err_stat, "flux_total_low": f_tot_low.clip(min=0 * f_tot_low.unit), "flux_total_high": f_tot_high})

            elif flux_errn is not None and flux_errp is not None:
                errn_stat, errp_stat = flux_errn * scale, flux_errp * scale
                if error_comb_mode == 'linear':
                    err_tot_low, err_tot_high = d_sys_low + errn_stat, d_sys_high + errp_stat
                elif error_comb_mode == 'quadrature':
                    err_tot_low, err_tot_high = np.sqrt(d_sys_low**2 + errn_stat**2), np.sqrt(d_sys_high**2 + errp_stat**2)
                else:
                    raise ValueError("error_comb_mode must be 'linear' or 'quadrature'")
                
                f_tot_low, f_tot_high = flux_out - err_tot_low, flux_out + err_tot_high
                result.update({
                    "flux_errn": errn_stat, "flux_errp": errp_stat,
                    "flux_total_low": f_tot_low.clip(min=0 * f_tot_low.unit), "flux_total_high": f_tot_high
                })

        # Output table
        if return_table:
            table = QTable()
            table["energy"], table["flux"] = energy, flux_out
            for key in ["flux_err", "flux_errn", "flux_errp", "flux_sys_low", "flux_sys_high", "flux_total_low", "flux_total_high"]:
                if key in result: table[key] = result[key]
            
            table["survival_probability"], table["survival_low"], table["survival_high"] = p_med, p_low, p_high
            result["table"] = table
            if output_ecsv: table.write(output_ecsv, format='ascii.ecsv', overwrite=True)

        return result
