#include <virtual-casing.h>
#include <virtual-casing.hpp>
#include <tuple>

#ifdef __cplusplus
extern "C" {
#endif

void GenerateVirtualCasingTestDataF(float* Bext, float* Bint, int NFP, bool half_period, long trg_Nt, long trg_Np, const float* X, long Nt, long Np) {
  std::vector<float> Bext_, Bint_;
  const std::vector<float> X_(X, X+3*Nt*Np);
  std::tie(Bext_, Bint_) = VirtualCasingTestData<float>::BFieldData(NFP, half_period, Nt, Np, X_, trg_Nt, trg_Np);
  SCTL_ASSERT((long)Bext_.size() == 3*trg_Nt*trg_Np);
  SCTL_ASSERT((long)Bint_.size() == 3*trg_Nt*trg_Np);
  for (long i = 0; i < 3*trg_Nt*trg_Np; i++) {
    Bext[i] = Bext_[i];
    Bint[i] = Bint_[i];
  }
}

void GenerateVirtualCasingTestDataD(double* Bext, double* Bint, int NFP, bool half_period, long trg_Nt, long trg_Np, const double* X, long Nt, long Np) {
  std::vector<double> Bext_, Bint_;
  const std::vector<double> X_(X, X+3*Nt*Np);
  std::tie(Bext_, Bint_) = VirtualCasingTestData<double>::BFieldData(NFP, half_period, Nt, Np, X_, trg_Nt, trg_Np);
  SCTL_ASSERT((long)Bext_.size() == 3*trg_Nt*trg_Np);
  SCTL_ASSERT((long)Bint_.size() == 3*trg_Nt*trg_Np);
  for (long i = 0; i < 3*trg_Nt*trg_Np; i++) {
    Bext[i] = Bext_[i];
    Bint[i] = Bint_[i];
  }
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

void VirtualCasingSetupF(int digits, int NFP, bool half_period, long Nt, long Np, const float* X, long src_Nt, long src_Np, long trg_Nt, long trg_Np, void* ctx) {
  VirtualCasing<float>& virtual_casing = *(VirtualCasing<float>*)ctx;
  const std::vector<float> X_(X, X+3*Nt*Np);
  virtual_casing.Setup(digits, NFP, half_period, Nt, Np, X_, src_Nt, src_Np, trg_Nt, trg_Np);
}
void VirtualCasingSetupD(int digits, int NFP, bool half_period, long Nt, long Np, const double* X, long src_Nt, long src_Np, long trg_Nt, long trg_Np, void* ctx) {
  VirtualCasing<double>& virtual_casing = *(VirtualCasing<double>*)ctx;
  const std::vector<double> X_(X, X+3*Nt*Np);
  virtual_casing.Setup(digits, NFP, half_period, Nt, Np, X_, src_Nt, src_Np, trg_Nt, trg_Np);
}

void VirtualCasingComputeBextF(float* Bext, const float* B, long Nt, long Np, const void* ctx) {
  const std::vector<float> B_(B, B+3*Nt*Np);
  const VirtualCasing<float>& virtual_casing = *(VirtualCasing<float>*)ctx;
  std::vector<float> Bext_ = virtual_casing.ComputeBext(B_);
  for (long i = 0; i < (long)Bext_.size(); i++) Bext[i] = Bext_[i];
}
void VirtualCasingComputeBextD(double* Bext, const double* B, long Nt, long Np, const void* ctx) {
  const std::vector<double> B_(B, B+3*Nt*Np);
  const VirtualCasing<double>& virtual_casing = *(VirtualCasing<double>*)ctx;
  std::vector<double> Bext_ = virtual_casing.ComputeBext(B_);
  for (long i = 0; i < (long)Bext_.size(); i++) Bext[i] = Bext_[i];
}

#ifdef __cplusplus
}
#endif
