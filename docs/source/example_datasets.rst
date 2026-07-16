Example Datasets
================

The following datasets are provided for reproducing the examples shown throughout the GPSP documentation.

Galactic Center PeVatron (H.E.S.S.)
-----------------------------------

* Flux points measured by H.E.S.S. for the GC PeVatron Diffuse Emission region (:download:`gc_pacman_hess.dat <example_datasets/gc_pacman_hess.dat>`). The H.E.S.S. flux points were taken from "Abramowski, A., et al. 2016, Nature, 531, 476, doi: 10.1038/nature17147"

Expected gamma-ray spectra from hadronic proton-proton interactions followed by neutral-pion decay are obtained by fitting H.E.S.S. flux data to Exponential Cutoff Power-Law proton spectral model with NaimaSpectralModel implementation of Gammapy. Three representative fixed proton cutoff energies of 3 PeV, 10 PeV, and 30 PeV are considered. The gamma-ray spectral models can be downloaded below.

* Expected gamma-ray spectrum from pion decay for fixed proton cutoff of :math:`E_{p}` = 3 PeV (:download:`gc_pp_model_3PeV.ecsv <example_datasets/gc_pp_model_3PeV.ecsv>`)
* Expected gamma-ray spectrum from pion decay for fixed proton cutoff of :math:`E_{p}` = 10 PeV (:download:`gc_pp_model_10PeV.ecsv <example_datasets/gc_pp_model_10PeV.ecsv>`)
* Expected gamma-ray spectrum from pion decay for fixed proton cutoff of :math:`E_{p}` = 30 PeV (:download:`gc_pp_model_30PeV.ecsv <example_datasets/gc_pp_model_30PeV.ecsv>`)

Cygnus X-3 UHE Emission Region (LHAASO)
---------------------------------------

The measured and attenuation-corrected LHAASO flux points listed below were taken from Cao et al. (2025). The attenuation-corrected flux points were
derived using the Popescu et al. (2017) ISRF model.

* Flux points measured by LHAASO for the Cygnus X-3 region (:download:`Cygnus_X-3_LHAASO_FluxData_Uncorrected.dat <example_datasets/Cygnus_X-3_LHAASO_FluxData_Uncorrected.dat>`).
* Attenuation corrected flux points for the Cygnus X-3 region (:download:`Cygnus_X-3_LHAASO_FluxData_Corrected_Popescu2017.dat <example_datasets/Cygnus_X-3_LHAASO_FluxData_Corrected_Popescu2017.dat>`).

The best fit power-law (PL) model for the Cygnus X-3 spectrum is obtained by fitting the flux data to PL model with Gammapy package. The reference energy was fixed to :math:`E_{0}` = 50 TeV. The best fit parameters were found to be :math:`N_{0}` = (2.7 :math:`\pm` 0.8) :math:`\times` 10 :math:`TeV^{-1}` :math:`cm^{-2}` :math:`s^{-1}` and :math:`\Gamma` = 2.20 :math:`\pm` 0.21, fully consistent with the values reported by the LHAASO Collaboration.

* The best fit PL model for Cygnus X-3 spectrum (:download:`Cygnus_X-3_BestFitPL_Uncorrected.ecsv <example_datasets/Cygnus_X-3_BestFitPL_Uncorrected.ecsv>`)
