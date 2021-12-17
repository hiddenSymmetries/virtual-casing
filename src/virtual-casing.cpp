#include <virtual-casing.h>
#include <virtual-casing.hpp>

#ifdef __cplusplus
extern "C" {
#endif

void GenerateVirtualCasingTestDataF(float* Bext, float* Bint, long Nt, long Np, const float* X) {
  sctl::Vector<float> Bext_(3*Nt*Np, sctl::Ptr2Itr<float>(Bext, 3*Nt*Np), false);
  sctl::Vector<float> Bint_(3*Nt*Np, sctl::Ptr2Itr<float>(Bint, 3*Nt*Np), false);
  const sctl::Vector<float> X_(3*Nt*Np, (sctl::Iterator<float>)sctl::Ptr2ConstItr<float>(X, 3*Nt*Np), false);
  VirtualCasingTestData<float>::BFieldData(Bext_, Bint_, Nt, Np, X_);
}

void GenerateVirtualCasingTestDataD(double* Bext, double* Bint, long Nt, long Np, const double* X) {
  sctl::Vector<double> Bext_(3*Nt*Np, sctl::Ptr2Itr<double>(Bext, 3*Nt*Np), false);
  sctl::Vector<double> Bint_(3*Nt*Np, sctl::Ptr2Itr<double>(Bint, 3*Nt*Np), false);
  const sctl::Vector<double> X_(3*Nt*Np, (sctl::Iterator<double>)sctl::Ptr2ConstItr<double>(X, 3*Nt*Np), false);
  VirtualCasingTestData<double>::BFieldData(Bext_, Bint_, Nt, Np, X_);
}

void* VirtualCasingCreateContextF() {
  return new VirtualCasing<float>;
}
void* VirtualCasingCreateContextD() {
  return new VirtualCasing<double>;
}

void VirtualCasingDestroyContextF(void** ctx) {
  delete (VirtualCasing<float>*)*ctx;
  ctx[0] = nullptr;
}
void VirtualCasingDestroyContextD(void** ctx) {
  delete (VirtualCasing<double>*)*ctx;
  ctx[0] = nullptr;
}

void VirtualCasingSetSurfaceF(long Nt, long Np, const float* X, void* ctx) {
  VirtualCasing<float>& virtual_casing = *(VirtualCasing<float>*)ctx;
  const sctl::Vector<float> X_(3*Nt*Np, (sctl::Iterator<float>)sctl::Ptr2ConstItr<double>(X, 3*Nt*Np), false);
  virtual_casing.SetSurface(Nt, Np, X_);
}
void VirtualCasingSetSurfaceD(long Nt, long Np, const double* X, void* ctx) {
  VirtualCasing<double>& virtual_casing = *(VirtualCasing<double>*)ctx;
  const sctl::Vector<double> X_(3*Nt*Np, (sctl::Iterator<double>)sctl::Ptr2ConstItr<double>(X, 3*Nt*Np), false);
  virtual_casing.SetSurface(Nt, Np, X_);
}

void VirtualCasingSetAccuracyF(int digits, void* ctx) {
  VirtualCasing<float>& virtual_casing = *(VirtualCasing<float>*)ctx;
  virtual_casing.SetAccuracy(digits);
}
void VirtualCasingSetAccuracyD(int digits, void* ctx) {
  VirtualCasing<double>& virtual_casing = *(VirtualCasing<double>*)ctx;
  virtual_casing.SetAccuracy(digits);
}

void VirtualCasingComputeBextF(float* Bext, const float* B, long Nt, long Np, const void* ctx) {
  sctl::Vector<float> Bext_(3*Nt*Np, sctl::Ptr2Itr<float>(Bext, 3*Nt*Np), false);
  const sctl::Vector<float> B_(3*Nt*Np, (sctl::Iterator<float>)sctl::Ptr2ConstItr<float>(B, 3*Nt*Np), false);
  const VirtualCasing<float>& virtual_casing = *(VirtualCasing<float>*)ctx;
  virtual_casing.ComputeBext(Bext_, B_);
}
void VirtualCasingComputeBextD(double* Bext, const double* B, long Nt, long Np, const void* ctx) {
  sctl::Vector<double> Bext_(3*Nt*Np, sctl::Ptr2Itr<double>(Bext, 3*Nt*Np), false);
  const sctl::Vector<double> B_(3*Nt*Np, (sctl::Iterator<double>)sctl::Ptr2ConstItr<double>(B, 3*Nt*Np), false);
  const VirtualCasing<double>& virtual_casing = *(VirtualCasing<double>*)ctx;
  virtual_casing.ComputeBext(Bext_, B_);
}

#ifdef __cplusplus
}
#endif
