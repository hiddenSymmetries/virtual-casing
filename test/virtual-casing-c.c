#include <virtual-casing.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void test(int digits, int NFP, bool half_period, long Nt, long Np, long src_Nt, long src_Np, long trg_Nt, long trg_Np) {
  // Construct the surface
  double* X = (double*)malloc(3*Nt*Np*sizeof(double));
  for (long t = 0; t < Nt; t++) { // toroidal direction
    double phi = 2*M_PI * (half_period?(2*t+1)/(double)4:(double)t) / Nt/NFP;
    for (long p = 0; p < Np; p++) { // poloidal direction
      double x = (2 + 0.5*cos(2*M_PI*p/Np)) * cos(phi);
      double y = (2 + 0.5*cos(2*M_PI*p/Np)) * sin(phi);
      double z = 0.5*sin(2*M_PI*p/Np);
      X[0*Nt*Np+t*Np+p] = x;
      X[1*Nt*Np+t*Np+p] = y;
      X[2*Nt*Np+t*Np+p] = z;
    }
  }

  // Generate B fields for testing virtual-casing principle
  double* Bext_trg = (double*)malloc(3*trg_Nt*trg_Np*sizeof(double));
  double* Bint_trg = (double*)malloc(3*trg_Nt*trg_Np*sizeof(double));
  double* Bext = (double*)malloc(3*src_Nt*src_Np*sizeof(double));
  double* Bint = (double*)malloc(3*src_Nt*src_Np*sizeof(double));
  double* B = (double*)malloc(3*src_Nt*src_Np*sizeof(double));
  GenerateVirtualCasingTestDataD(Bext_trg, Bint_trg, NFP, half_period, trg_Nt, trg_Np, X, Nt, Np);
  GenerateVirtualCasingTestDataD(Bext, Bint, NFP, half_period, src_Nt, src_Np, X, Nt, Np);
  for (long i = 0; i < 3*src_Nt*src_Np; i++) B[i] = Bext[i] + Bint[i];

  // Setup
  void* virtual_casing = VirtualCasingCreateContextD();
  VirtualCasingSetupD(digits, NFP, half_period, Nt, Np, X, src_Nt, src_Np, trg_Nt, trg_Np, virtual_casing);

  // Compute Bext field
  double* Bext_ = (double*)malloc(3*trg_Nt*trg_Np*sizeof(double));
  VirtualCasingComputeBextD(Bext_, B, src_Nt, src_Np, virtual_casing);


  // print error
  double max_err = 0, max_val = 0;
  for (long i = 0; i < 3*trg_Nt*trg_Np; i++) {
    double Berr = Bext_trg[i] - Bext_[i];
    max_err = (max_err > fabs(Berr) ? max_err : fabs(Berr));
  }
  for (long i = 0; i < 3*src_Nt*src_Np; i++) {
    max_val = (max_val > fabs(B[i]) ? max_val : fabs(B[i]));
  }
  printf("Maximum relative error: %e\n", max_err/max_val);

  // Destroy context
  VirtualCasingDestroyContextD(&virtual_casing);

  // Free memory
  free(X);
  free(B);
  free(Bint);
  free(Bext);
  free(Bint_trg);
  free(Bext_trg);
  free(Bext_);
}

int main() {
  for (long digits = 3; digits <= 12; digits+=3) {
    test(digits, 5, false, 1, 4, 2*digits, 7*digits, 20, 20);
    test(digits, 5,  true, 1, 4, 1*digits, 7*digits, 10, 20);
  }
  return 0;
}

