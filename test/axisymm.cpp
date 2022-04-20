#include <virtual-casing.hpp>

template <class Real> void test(const Real R0, const Real a, const Real kappa, int digits, int NFP, long Nt, long Np, long src_Nt, long src_Np, long trg_Nt, long trg_Np) {
  // Construct the surface
  std::vector<Real> X(3*Nt*Np);
  //X = VirtualCasingTestData<Real>::SurfaceCoordinates(NFP, Nt, Np, biest::SurfType::W7X_);
  for (long t = 0; t < Nt; t++) { // toroidal direction
    for (long p = 0; p < Np; p++) { // poloidal direction
      Real theta = 2*M_PI*p/Np;
      Real phi   = 2*M_PI*t/(NFP*Nt);

      Real R = sqrt(R0*R0 + 2*a*R0*cos(theta));
      Real x = R * cos(phi);
      Real y = R * sin(phi);
      Real z = kappa*a*R0/R*sin(theta);

      X[(0*Nt+t)*Np+p] = x;
      X[(1*Nt+t)*Np+p] = y;
      X[(2*Nt+t)*Np+p] = z;
    }
  }

  // Generate B field
  std::vector<Real> B(3*src_Nt*src_Np);
  for (long t = 0; t < src_Nt; t++) { // toroidal direction
    for (long p = 0; p < src_Np; p++) { // poloidal direction
      Real theta = 2*M_PI*p/src_Np;
      Real phi   = 2*M_PI*t/(NFP*src_Nt);

      B[(0*src_Nt+t)*src_Np+p] = -sin(theta)/3 * cos(phi);
      B[(1*src_Nt+t)*src_Np+p] = -sin(theta)/3 * sin(phi);
      B[(2*src_Nt+t)*src_Np+p] = (1.0/3)*cos(theta) + (1.0/9)*sin(theta)*sin(theta) / (1+2.0/3*cos(theta));
    }
  }

  // Setup
  VirtualCasing<Real> virtual_casing;
  virtual_casing.Setup(digits, NFP, Nt, Np, X, src_Nt, src_Np, trg_Nt, trg_Np);

  // Compute Bext field
  std::vector<Real> Bext_ = virtual_casing.ComputeBext(B);

  for (long k = 0; k < 3; k++) { // Print Bext_
    for (long t = 0; t < trg_Nt; t++) {
      for (long p = 0; p < trg_Np; p++) {
        printf("%10.5f", Bext_[(k*trg_Nt+t)*trg_Np+p]);
      }
      std::cout<<'\n';
    }
  }

  // print error
  //auto Berr = Bext;
  //Real max_err = 0, max_val = 0;
  //for (long i = 0; i < (long)Berr.size(); i++) Berr[i] -= Bext_[i];
  //for (const auto& x:B   ) max_val = std::max<Real>(max_val,fabs(x));
  //for (const auto& x:Berr) max_err = std::max<Real>(max_err,fabs(x));
  //std::cout<<"Maximum relative error: "<<max_err/max_val<<'\n';
}

int main() {
  long digits = 10;
  double R0 = 1, a = 1./3, kappa = 1;

  test<double>(R0, a, kappa, digits, 4, 1, 10, 1, 10, 1, 10);

  return 0;
}
