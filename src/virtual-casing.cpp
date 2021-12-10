#include <virtual-casing.hpp>
#include <virtual-casing-test-data.hpp>

template <class Real> void test(long Nt, long Np, int digits, biest::SurfType surf_type = biest::SurfType::AxisymNarrow) {
  sctl::Vector<Real> X, Bext, Bint, B;
  GenerateVirtualCasingTestData<Real>(X, Bext, Bint, Nt, Np, surf_type);
  B = Bext + Bint;

  VirtualCasing<Real> virtual_casing;
  virtual_casing.SetSurface(Nt, Np, X);
  virtual_casing.SetAccuracy(digits);

  sctl::Vector<Real> Berr, Bext_;
  virtual_casing.ComputeBext(Bext_, B);
  Berr = Bext - Bext_;

  Real max_err = 0, max_val = 0;
  for (const auto& x:B   ) max_val = std::max<Real>(max_val,fabs(x));
  for (const auto& x:Berr) max_err = std::max<Real>(max_err,fabs(x));
  std::cout<<"Maximum relative error: "<<max_err/max_val<<'\n';
}

int main() {
  for (long digits = 1; digits < 12; digits++) {
    test<double>(700, 140, digits, biest::SurfType::AxisymNarrow);
    //test<double>(1100, 140, digits, biest::SurfType::W7X_);
  }
  return 0;
}

