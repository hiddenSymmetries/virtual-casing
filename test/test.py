#!/usr/bin/env python

import unittest
import numpy as np
import virtual_casing as vc

class VirtualCasingTests(unittest.TestCase):
    def subtest_w7x(self, digits, nfp, half_period, nphi, ntheta, surf_type, src_nphi, src_ntheta, trg_nphi, trg_ntheta, expected_error):
        X = vc.VirtualCasingTestData.surface_coordinates(nfp, half_period, nphi, ntheta, surf_type)

        B_external_true_trg_resolution, B_internal_true_trg_resolution = vc.VirtualCasingTestData.magnetic_field_data(nfp, half_period, nphi, ntheta, X, trg_nphi, trg_ntheta)
        B_external_true_src_resolution, B_internal_true_src_resolution = vc.VirtualCasingTestData.magnetic_field_data(nfp, half_period, nphi, ntheta, X, src_nphi, src_ntheta)
        B_total_true_src_resolution = np.array(B_external_true_src_resolution) + np.array(B_internal_true_src_resolution)

        vcasing = vc.VirtualCasing()
        vcasing.setup(digits, nfp, half_period, nphi, ntheta, X, src_nphi, src_ntheta, trg_nphi, trg_ntheta)

        B_external_computed = vcasing.compute_external_B(B_total_true_src_resolution)
        B_err = np.array(B_external_computed) - np.array(B_external_true_trg_resolution)
        max_val = np.abs(B_total_true_src_resolution).max()
        max_err = np.abs(B_err).max()
        max_rel_err = max_err / max_val
        print(f"For digits = {digits} and half_period = {half_period}, maximum relative error = {max_rel_err}")
        self.assertLess(max_rel_err, expected_error)

    def test_w7x(self):
        """
        Run the test shown in Figure 2 of Malhotra et al, Plasma Physics
        and Controlled Fusion 62, 024004 (2020).
        """
        expected_errors = [2e-4, 3e-7, 2e-9, 3e-11]
        for j in range(4):
            digits = 3 + 3 * j
            for half_period in [True, False]:
                with self.subTest(digits=digits, half_period=half_period):
                    # Double the number of toroidal grid points if we do not exploit stellarator symmetry:
                    nphi_factor = 1 if half_period else 2

                    self.subtest_w7x(digits, 5, half_period, nphi_factor * 10, 20,
                                     vc.SurfType.W7X_,
                                     nphi_factor * 6 * digits, 32 * digits,
                                     nphi_factor * 20, 40,
                                     expected_errors[j])

    def test_test_data_phi_grid(self):
        """
        For VirtualCasingTestData.surface_coordinates, make sure the phi
        grid behaves as expected for both half_period=True and False.
        """
        for half_period in [True, False]:
            for surf_type in [vc.SurfType.W7X_]:
                nfp = 5
                nphi = 30
                ntheta = 25
                x1d = vc.VirtualCasingTestData.surface_coordinates(nfp, half_period, nphi, ntheta, surf_type)
                x1d = np.array(x1d)

                # Unpack 1D array results:
                x3d = np.zeros((nphi, ntheta, 3))
                for jxyz in range(3):
                    x3d[:, :, jxyz] = x1d[jxyz * nphi * ntheta: (jxyz + 1) * nphi * ntheta].reshape((nphi, ntheta), order='C')

                # Check order:
                index = 0
                for jxyz in range(3):
                    for jphi in range(nphi):
                        for jtheta in range(ntheta):
                            np.testing.assert_allclose(x1d[index], x3d[jphi, jtheta, jxyz])
                            index += 1

                phi_from_vc = np.arctan2(x3d[:, :, 1], x3d[:, :, 0])

                if half_period:
                    phi1d = np.linspace(0, np.pi / nfp, nphi, endpoint=False)
                    phi1d += 0.5 * (phi1d[1] - phi1d[0])
                else:
                    phi1d = np.linspace(0, 2 * np.pi / nfp, nphi, endpoint=False)
                theta1d = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)

                theta2d, phi2d = np.meshgrid(theta1d, phi1d)

                np.testing.assert_allclose(phi2d, phi_from_vc, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
