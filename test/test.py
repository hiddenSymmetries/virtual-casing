#!/usr/bin/env python

import virtual_casing as vc
import numpy as np

def test(digits, NFP, half_period, Nt, Np, surf_type, src_Nt, src_Np, trg_Nt, trg_Np):
    X = vc.VirtualCasingTestData.surface_coordinates(NFP, half_period, Nt, Np, surf_type)

    Bext, Bint = vc.VirtualCasingTestData.magnetic_field_data(NFP, half_period, Nt, Np, X, trg_Nt, trg_Np)
    Bext_, Bint_ = vc.VirtualCasingTestData.magnetic_field_data(NFP, half_period, Nt, Np, X, src_Nt, src_Np)
    B = np.array(Bext_) + np.array(Bint_)

    vcasing = vc.VirtualCasing()
    vcasing.setup(digits, NFP, half_period, Nt, Np, X, src_Nt, src_Np, trg_Nt, trg_Np)

    Bext_ = vcasing.compute_external_B(B)
    Berr = np.array(Bext) - np.array(Bext_)
    max_val = np.abs(B).max()
    max_err = np.abs(Berr).max()
    print(f"For digits = {digits}, maximum relative error = {max_err/max_val}")

if __name__ == '__main__':
    for i in range(3, 13, 3):
        test(i, 5, true, 10, 20, vc.SurfType.W7X_, 6*i, 32*i, 20, 40)
