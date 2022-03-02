import virtual_casing as vc
import numpy as np

def test(Nt, Np, digits, surf_type):
    X = vc.VirtualCasingTestData.surface_coordinates(Nt, Np, surf_type)

    Bext, Bint = vc.VirtualCasingTestData.magnetic_field_data(Nt, Np, X)

    B = np.array(Bext) + np.array(Bint)

    vcasing = vc.VirtualCasing()
    vcasing.set_surface(Nt, Np, X)
    vcasing.set_accuracy(digits)

    Bext_ = vcasing.compute_external_B(list(B))
    Berr = np.array(Bext) - np.array(Bext_)
    max_val = np.abs(B).max()
    max_err = np.abs(Berr).max()
    print(f"Maximum relative error: {max_err/max_val}")

if __name__ == '__main__':
    for i in range(1, 12):
        test(40*i, 10*i, i, vc.SurfType.AxisymNarrow)
