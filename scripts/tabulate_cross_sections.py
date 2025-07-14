from pathlib import Path

import h5py
import numpy as np
from pscatter_cross_section import diff_cross, diff_cross_polarized
from scipy.integrate import simpson

### Global configurations
USE_LOG_RANGE = True
GRIDSIZE = 250
OMEGA_MIN = 5.0e-2
OMEGA_MAX = 1.0e4
POL_KEYS = ["++++", "+++-", "++--"]

if __name__ == "__main__":

    if USE_LOG_RANGE:
        ### Logarithmic spacing chosen for omega = ZMF energy of each photon
        log_range = np.linspace(np.log10(OMEGA_MIN),
                                np.log10(OMEGA_MAX),
                                GRIDSIZE)
        energy = np.power(10.0, log_range)
    else:
        ### Linear spacing chosen for omega = ZMF energy of each photon
        energy = np.linspace(OMEGA_MIN, OMEGA_MAX, GRIDSIZE)
    
    s = 4.0 * energy**2 ### Mandelstam s
    angle = np.linspace(0.01, 3.13, GRIDSIZE) ### (diverges at 0 or pi)
    p_axis = np.linspace(0.0, 1.0, GRIDSIZE)
    angle_mesh, energy_mesh = np.meshgrid(angle, energy)

    ### Allocating storage
    ### Calculation of differential and total cross section for unpolarized photons
    dsig_unpol = diff_cross(energy_mesh, angle_mesh)
    tot_sig_unpol = 2.0*np.pi * simpson(y=dsig_unpol*np.sin(angle_mesh),
                                        x=angle)
    ### Calculation of differential and total cross section for polarized photons
    dsig_pol = diff_cross_polarized(energy_mesh, angle_mesh)
    tot_sig_pol = {key: 2.0*np.pi * simpson(y=dsig_pol[key]*np.sin(angle_mesh), x=angle)
                   for key in POL_KEYS}

    ### Invert CDF for MC sampling
    def compute_inv_cdf(ds):
        cdf = np.zeros_like(ds)
        cdf_inv = np.zeros_like(ds)
        for i, angle_array in enumerate(ds):
            data = angle_array * np.sin(angle)
            cdf[i, 0] = data[0]
            for j in range(1, len(data)):
                cdf[i, j] = cdf[i, j-1] + data[j]
            cdf[i, :] = cdf[i, :] / cdf[i, len(data) - 1]
            cdf_inv[i, :] = np.interp(p_axis, cdf[i, :], angle)
        return cdf_inv

    cdf_inv_unpol = compute_inv_cdf(dsig_unpol)
    cdf_inv_pol = {key: compute_inv_cdf(dsig_pol[key]) for key in POL_KEYS}

    ### Writing tabulations to HDF5 file
    data_path = Path("data")
    hf = h5py.File(data_path / 'PhotonScatter_total.h5', 'w')
    hf.create_dataset('s', data=s)
    hf.create_dataset('unpolarized', data=tot_sig_unpol)

    for key in POL_KEYS:
        hf.create_dataset(f'polarized_{key}', data=tot_sig_pol[key])

    hf = h5py.File(data_path / 'PhotonScatter_diff.h5', 'w')
    hf.create_dataset('s', data=s)
    hf.create_dataset('p', data=p_axis)
    hf.create_dataset('angle_unpolarized', data=cdf_inv_unpol)
    for key in POL_KEYS:
        hf.create_dataset(f'angle_polarized_{key}', data=cdf_inv_pol[key])
