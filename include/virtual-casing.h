#ifndef _BIEST_VIRTUAL_CASING_TEST_DATA_H_
#define _BIEST_VIRTUAL_CASING_TEST_DATA_H_

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Create virtual-casing context and return its pointer.
 */
void* VirtualCasingCreateContextF();
void* VirtualCasingCreateContextD();

/**
 * Destroy virtual-casing context.
 *
 * @param[in,out] virtual-casing context pointer.
 */
void VirtualCasingDestroyContextF(void** ctx);
void VirtualCasingDestroyContextD(void** ctx);

/**
 * Setup the VirtualCasing instance.
 *
 * @param[in] digits number of decimal digits of accuracy.
 *
 * @param[in] NFP number of toroidal field periods. The surface as well as
 * the magnetic field must have this toroidal periodic symmetry.
 *
 * @param[in] half_period whether the surface and data are defined on half field period.
 *
 * @param[in] Nt surface discretization order in toroidal direction (in one field period).
 *
 * @param[in] Np surface discretization order in poloidal direction.
 *
 * @param[in] X the surface coordinates in the order {x11, x12, ..., x1Np,
 * x21, x22, ... , xNtNp, y11, ... , z11, ...}.
 *
 * @param[in] src_Nt input B-field discretization order in toroidal direction (in one field period).
 *
 * @param[in] src_Np input B-field discretization order in poloidal direction.
 *
 * @param[in] trg_Nt output Bext-field discretization order in toroidal direction (in one field period).
 *
 * @param[in] trg_Np output Bext-field discretization order in poloidal direction.
 *
 * @param[in,out] virtual-casing context pointer.
 *
 * The grids for the input and output data differ in an important
 * way depending on whether or not stellarator symmetry is
 * exploited, i.e. on half_period. For this discussion, consider
 * the toroidal angle phi and poloidal angle theta to have period
 * 1 (not 2*pi).
 *
 * If you do not exploit stellarator symmetry, then half_period
 * is set to false.  In this case the grids in the toroidal angles
 * begin at phi=0.  The grid spacing is 1 / (NFP * Nt), and there
 * is no point at the symmetry plane phi = 1 / NFP.
 *
 * If you do wish to exploit stellarator symmetry, set half_period
 * to true. In this case the toroidal grids are each shifted by
 * half a grid point, so there is no grid point at phi = 0. The
 * phi grid for the surface shape has points at 0.5 / (NFP * Nt),
 * 1.5 / (NFP * Nt), ..., (Nt - 0.5) / (NFP * Nt). The phi grid
 * for the output B_external has grid points at 0.5 / (NFP *
 * trg_Nt), 1.5 / (NFP * trg_Nt), ..., (trg_Nt - 0.5) / (NFP *
 * Nt), and similarly for the input B field with trg -> src.
 *
 * The rationale for these conventions is that in both the
 * stellarator-symmetric and non-stellarator-symmetric cases,
 * integration over the surface can be achieved with spectral
 * accuracy on these grids using uniform weights.
 *
 * Regardless of half_period, the poloidal grid always ranges
 * uniformly over [0, 1), with the first grid point at theta = 0,
 * and no grid point at theta = 1.
 *
 * The resolution parameters for the surface shape (Nt, Np), input
 * magnetic field (src_Nt, src_Np), and output external field
 * (trg_Nt, trg_Np) do not need to be related to each other in any
 * particular way. For example, there is no performance penalty if
 * src_Nt and trg_Nt are relatively prime.
 */
void VirtualCasingSetupF(int digits, int NFP, bool half_period, long Nt, long Np, const float* X, long src_Nt, long src_Np, long trg_Nt, long trg_Np, void* ctx);
void VirtualCasingSetupD(int digits, int NFP, bool half_period, long Nt, long Np, const double* X, long src_Nt, long src_Np, long trg_Nt, long trg_Np, void* ctx);

/**
 * Recover the Bext component from the total field B = Bext + Bint by
 * applying the virtual-casing principle:
 * Bext = B/2 + gradG[B . n] + BiotSavart[n x B]
 *
 * @param[out] Bext the component of magnetic field on the surface due to
 * currents in the exterior of the surface, computed using the
 * virtual-casing principle.
 *
 * @param[in] B the total magnetic field on the surface due to all currents.
 * B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ... , BxNtNp, By11, ... , Bz11, ...},
 * where Nt and Np are the number of discretizations in toroidal and
 * poloidal directions.
 */

/**
 * Recover the Bext component from the total field B = Bext + Bint by
 * applying the virtual-casing principle:
 * Bext = B/2 + gradG[B . n] + BiotSavart[n x B]
 *
 * @param[out] Bext the component of magnetic field on the surface due to
 * currents in the exterior of the surface, computed using the
 * virtual-casing principle.
 *
 * @param[in] B the total magnetic field on the surface due to all currents.
 * B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ... , BxNtNp, By11, ... , Bz11, ...},
 * where Nt and Np are the number of discretizations in toroidal and
 * poloidal directions.
 *
 * @param[in] Nt input B-field discretization order in toroidal direction (in one field period).
 *
 * @param[in] Np input B-field discretization order in poloidal direction.
 *
 * @param[in] virtual-casing context pointer.
 */
void VirtualCasingComputeBextF(float* Bext, const float* B, long Nt, long Np, const void* ctx);
void VirtualCasingComputeBextD(double* Bext, const double* B, long Nt, long Np, const void* ctx);

/**
 * Generate B field data to be used for testing virtual-casing.
 *
 * Generate B field data to be used with class VirtualCasing.
 *
 * @param[out] Bext magnetic field generated by an internal current loop.
 *
 * @param[out] Bint magnetic field generated by an external current loop.
 *
 * @param[in] NFP number of toroidal field periods.
 *
 * @param[in] half_period whether the result should be on half field period.
 *
 * @param[in] trg_Nt output B-field discretization order in toroidal direction (in one field period).
 *
 * @param[in] trg_Np output B-field discretization order in poloidal direction.
 *
 * @param[in] X the surface coordinates in the order {x11, x12, ..., x1Np,
 * x21, x22, ... , xNtNp, y11, ... , z11, ...}.
 *
 * @param[in] Nt surface discretization order in toroidal direction (in one field period).
 *
 * @param[in] Np surface discretization order in poloidal direction.
 *
 * @return Bext and Bint, magnetic fields generated by an internal current loop and
 *  an external current loop respectively.
 */
void GenerateVirtualCasingTestDataF(float* Bext, float* Bint, int NFP, bool half_period, long trg_Nt, long trg_Np, const float* X, long Nt, long Np);
void GenerateVirtualCasingTestDataD(double* Bext, double* Bint, int NFP, bool half_period, long trg_Nt, long trg_Np, const double* X, long Nt, long Np);

#ifdef __cplusplus
}
#endif

#endif
