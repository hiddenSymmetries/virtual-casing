#include <virtual-casing.hpp>

template <class Real> void test(int digits, int NFP, long Nt, long Np, biest::SurfType surf_type, long src_Nt, long src_Np, long trg_Nt, long trg_Np) {
  // Construct the surface
  std::vector<Real> X(3*Nt*Np);
  X = VirtualCasingTestData<Real>::SurfaceCoordinates(NFP, Nt, Np, surf_type);
  //for (long t = 0; t < Nt; t++) { // toroidal direction
  //  for (long p = 0; p < Np; p++) { // poloidal direction
  //    Real x = (2 + 0.5*cos(2*M_PI*p/Np)) * cos(2*M_PI*t/(NFP*Nt));
  //    Real y = (2 + 0.5*cos(2*M_PI*p/Np)) * sin(2*M_PI*t/(NFP*Nt));
  //    Real z = 0.5*sin(2*M_PI*p/Np);
  //    X[(0*Nt+t)*Np+p] = x;
  //    X[(1*Nt+t)*Np+p] = y;
  //    X[(2*Nt+t)*Np+p] = z;
  //  }
  //}

  // Generate B fields for testing virtual-casing principle
  std::vector<Real> B, Bext;
  { // Set B, Bext
    std::vector<Real> Bint_, Bext_;
    std::tie(Bext, std::ignore) = VirtualCasingTestData<Real>::BFieldData(NFP, Nt, Np, X, trg_Nt, trg_Np);
    std::tie(Bext_, Bint_) = VirtualCasingTestData<Real>::BFieldData(NFP, Nt, Np, X, src_Nt, src_Np);
    const auto B_ = sctl::Vector<Real>(Bint_) + sctl::Vector<Real>(Bext_);
    B.assign(B_.begin(), B_.end());
  }

  // Setup
  VirtualCasing<Real> virtual_casing;
  virtual_casing.Setup(digits, NFP, Nt, Np, X, src_Nt, src_Np, trg_Nt, trg_Np);

  // Compute Bext field
  auto Bext_ = virtual_casing.ComputeBext(B);

  // print error
  auto Berr = Bext;
  Real max_err = 0, max_val = 0;
  for (long i = 0; i < (long)Berr.size(); i++) Berr[i] -= Bext_[i];
  for (const auto& x:B   ) max_val = std::max<Real>(max_val,fabs(x));
  for (const auto& x:Berr) max_err = std::max<Real>(max_err,fabs(x));
  std::cout<<"Maximum relative error: "<<max_err/max_val<<'\n';
}

int main() {
  for (long digits = 3; digits <= 12; digits+=3) {
    //test<double>(digits, 5, 1, 4, biest::SurfType::AxisymNarrow, 2*digits, 7*digits, 20, 20);
    test<double>(digits, 5, 20, 20, biest::SurfType::W7X_, 12*digits, 32*digits, 40, 40);
  }
  return 0;
}
