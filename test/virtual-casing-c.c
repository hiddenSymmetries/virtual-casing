#include <virtual-casing.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void test(long Nt, long Np, int digits) {
  // Construct the surface
  double* X = (double*)malloc(3*Nt*Np*sizeof(double));
  for (long t = 0; t < Nt; t++) { // toroidal direction
    for (long p = 0; p < Np; p++) { // poloidal direction
      double x = (2 + 0.5*cos(2*M_PI*p/Np)) * cos(2*M_PI*t/Nt);
      double y = (2 + 0.5*cos(2*M_PI*p/Np)) * sin(2*M_PI*t/Nt);
      double z = 0.5*sin(2*M_PI*p/Np);
      X[0*Nt*Np+t*Np+p] = x;
      X[1*Nt*Np+t*Np+p] = y;
      X[2*Nt*Np+t*Np+p] = z;
    }
  }

  // Generate B fields for testing virtual-casing principle
  double* Bext = (double*)malloc(3*Nt*Np*sizeof(double));
  double* Bint = (double*)malloc(3*Nt*Np*sizeof(double));
  double* B = (double*)malloc(3*Nt*Np*sizeof(double));
  GenerateVirtualCasingTestDataD(&Bext[0], &Bint[0], Nt, Np, &X[0]);
  for (long i = 0; i < 3*Nt*Np; i++) B[i] = Bext[i] + Bint[i];

  // Setup
  void* virtual_casing = VirtualCasingCreateContextD();
  VirtualCasingSetSurfaceD(Nt, Np, X, virtual_casing);
  VirtualCasingSetAccuracyD(digits, virtual_casing);

  // Compute Bext field
  double* Bext_ = (double*)malloc(3*Nt*Np*sizeof(double));
  VirtualCasingComputeBextD(Bext_, B, Nt, Np, virtual_casing);


  // print error
  double max_err = 0, max_val = 0;
  for (long i = 0; i < 3*Nt*Np; i++) {
    double Berr = Bext[i] - Bext_[i];
    max_val = (max_val > fabs(B[i]) ? max_val : fabs(B[i]));
    max_err = (max_err > fabs(Berr) ? max_err : fabs(Berr));
  }
  printf("Maximum relative error: %e\n", max_err/max_val);

  // Destroy context
  VirtualCasingDestroyContextD(&virtual_casing);

  // Free memory
  free(X);
  free(B);
  free(Bint);
  free(Bext);
  free(Bext_);
}

int main() {
  for (long digits = 1; digits < 12; digits++) {
    test(20*digits, 5*digits, digits);
  }
  return 0;
}

