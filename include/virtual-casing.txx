
template <class Real, sctl::Integer COORD_DIM, sctl::Integer KDIM0, sctl::Integer KDIM1> class BIOpWrapper {
  public:

    BIOpWrapper(const sctl::Comm& comm) : biop(nullptr), comm_(comm) {}

    ~BIOpWrapper() {
      if (biop) biop_delete(&biop);
    }

    void SetupSingular(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, sctl::Integer digits) {
      Real cond = 1;
      { // Set cond
        sctl::Vector<Real> dX;
        SCTL_ASSERT(Svec.Dim() == 1);
        sctl::StaticArray<sctl::Long,2> SurfDim{Svec[0].NTor(),Svec[0].NPol()};
        biest::SurfaceOp<Real> SurfOp(comm_, SurfDim[0], SurfDim[1]);
        SurfOp.Grad2D(dX, Svec[0].Coord());

        sctl::Matrix<Real> M(2,2), U, S, Vt;
        const sctl::Long N = SurfDim[0] * SurfDim[1];
        for (sctl::Long i = 0; i < N; i++) {
          for (sctl::Integer k0 = 0; k0 < 2; k0++) {
            for (sctl::Integer k1 = 0; k1 < 2; k1++) {
              Real dot_prod = 0;
              for (sctl::Integer j = 0; j < COORD_DIM; j++) {
                dot_prod += dX[(j*2+k0)*N+i] * dX[(j*2+k1)*N+i] / SurfDim[k0] / SurfDim[k1];
              }
              M[k0][k1] = dot_prod;
            }
          }

          M.SVD(U,S,Vt);
          Real cond2 = std::max<Real>(S[0][0],S[1][1]) / std::min<Real>(S[0][0],S[1][1]);
          cond = std::max<Real>(cond, sctl::sqrt<Real>(cond2));
        }
      }
      if (cond > 4) {
        SCTL_WARN("The surface mesh is highly anisotropic! Quadrature generation will be very slow. Consider using a better surface discretization.");
        std::cout<<"Mesh anisotropy = "<<cond<<'\n';
      }
      SetupSingular_<1>(Svec, ker, digits*1.2*cond);

      if (0) {
      if (biop) biop_delete(&biop);
      sctl::Integer SurfDim = std::min(Svec[0].NPol(), Svec[0].NTor());

      if (digits >= 18 && SurfDim > 36*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 36, RDIM = 54*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >= 17 && SurfDim > 34*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 34, RDIM = 51*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >= 16 && SurfDim > 32*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 32, RDIM = 48*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >= 15 && SurfDim > 30*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 30, RDIM = 45*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >= 14 && SurfDim > 28*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 28, RDIM = 42*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >= 13 && SurfDim > 26*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 26, RDIM = 39*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >= 12 && SurfDim > 24*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 24, RDIM = 36*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >= 11 && SurfDim > 22*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 22, RDIM = 33*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >= 10 && SurfDim > 20*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 20, RDIM = 30*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >=  9 && SurfDim > 18*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 18, RDIM = 27*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >=  8 && SurfDim > 16*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 16, RDIM = 24*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >=  7 && SurfDim > 14*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 14, RDIM = 21*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >=  6 && SurfDim > 12*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 12, RDIM = 18*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >=  5 && SurfDim > 10*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM = 10, RDIM = 15*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >=  4 && SurfDim > 8*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM =  8, RDIM = 12*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >=  3 && SurfDim > 6*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM =  6, RDIM =  9*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >=  2 && SurfDim > 6*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM =  6, RDIM =  6*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else if (digits >=  1 && SurfDim > 6*2+1) {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM =  6, RDIM =  3*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      } else {
        static constexpr sctl::Integer UPSAMPLE = 1, PDIM =  6, RDIM =  1*3;
        biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
        biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
        biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;
      }

      biop = biop_build(Svec, ker, comm_);
      }
    }

    void Eval(sctl::Vector<Real>& U, const sctl::Vector<Real>& F) {
      biop_eval(U, F, biop);
    }

  private:

    template <sctl::Integer UPSAMPLE, sctl::Integer PDIM, sctl::Integer RDIM> void SetupSingular0(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker) {
      if (biop) biop_delete(&biop);

      biop_build = BIOpBuild<UPSAMPLE,PDIM,RDIM>;
      biop_delete = BIOpDelete<UPSAMPLE,PDIM,RDIM>;
      biop_eval = BIOpEval<UPSAMPLE,PDIM,RDIM>;

      biop = biop_build(Svec, ker, comm_);
    }
    template <sctl::Integer UPSAMPLE> void SetupSingular_(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, sctl::Integer PDIM_) {
      const sctl::Integer SurfDim = std::min(Svec[0].NPol(), Svec[0].NTor());
      SCTL_ASSERT(Svec.Dim() == 1);
      if (PDIM_ >= 64 && SurfDim*UPSAMPLE > 64*2) {
        static constexpr sctl::Integer PDIM = 64, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 60 && SurfDim*UPSAMPLE > 60*2) {
        static constexpr sctl::Integer PDIM = 60, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 56 && SurfDim*UPSAMPLE > 56*2) {
        static constexpr sctl::Integer PDIM = 56, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 52 && SurfDim*UPSAMPLE > 52*2) {
        static constexpr sctl::Integer PDIM = 52, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 48 && SurfDim*UPSAMPLE > 48*2) {
        static constexpr sctl::Integer PDIM = 48, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 44 && SurfDim*UPSAMPLE > 44*2) {
        static constexpr sctl::Integer PDIM = 44, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 40 && SurfDim*UPSAMPLE > 40*2) {
        static constexpr sctl::Integer PDIM = 40, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 36 && SurfDim*UPSAMPLE > 36*2) {
        static constexpr sctl::Integer PDIM = 36, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 32 && SurfDim*UPSAMPLE > 32*2) {
        static constexpr sctl::Integer PDIM = 32, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 28 && SurfDim*UPSAMPLE > 28*2) {
        static constexpr sctl::Integer PDIM = 28, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 24 && SurfDim*UPSAMPLE > 24*2) {
        static constexpr sctl::Integer PDIM = 24, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 20 && SurfDim*UPSAMPLE > 20*2) {
        static constexpr sctl::Integer PDIM = 20, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 16 && SurfDim*UPSAMPLE > 16*2) {
        static constexpr sctl::Integer PDIM = 16, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 12 && SurfDim*UPSAMPLE > 12*2) {
        static constexpr sctl::Integer PDIM = 12, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >=  8 && SurfDim*UPSAMPLE >  8*2) {
        static constexpr sctl::Integer PDIM =  8, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (SurfDim*UPSAMPLE > 12) {
        static constexpr sctl::Integer PDIM =  6, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else {
        static constexpr sctl::Integer PDIM =  6, RDIM = PDIM*1.5;
        if (SurfDim > 6) SetupSingular0< 2,PDIM,RDIM>(Svec, ker);
        else if (SurfDim > 4) SetupSingular0< 3,PDIM,RDIM>(Svec, ker);
        else if (SurfDim == 4) SetupSingular0< 4,PDIM,RDIM>(Svec, ker);
        else if (SurfDim == 3) SetupSingular0< 6,PDIM,RDIM>(Svec, ker);
        else if (SurfDim == 2) SetupSingular0< 7,PDIM,RDIM>(Svec, ker);
        else if (SurfDim == 1) SetupSingular0<13,PDIM,RDIM>(Svec, ker);
        else SCTL_ASSERT(false);
      }
    }

    template <sctl::Integer UPSAMPLE, sctl::Integer PDIM, sctl::Integer RDIM> static void* BIOpBuild(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, const sctl::Comm& comm) {
      using BIOp = biest::BoundaryIntegralOp<Real,KDIM0,KDIM1,UPSAMPLE,PDIM,RDIM>;
      BIOp* biop = new BIOp(comm);
      biop[0].SetupSingular(Svec, ker);
      return biop;
    }
    template <sctl::Integer UPSAMPLE, sctl::Integer PDIM, sctl::Integer RDIM> static void BIOpDelete(void** self) {
      using BIOp = biest::BoundaryIntegralOp<Real,KDIM0,KDIM1,UPSAMPLE,PDIM,RDIM>;
      delete (BIOp*)self[0];
      self[0] = nullptr;
    }
    template <sctl::Integer UPSAMPLE, sctl::Integer PDIM, sctl::Integer RDIM> static void BIOpEval(sctl::Vector<Real>& U, const sctl::Vector<Real>& F, void* self) {
      using BIOp = biest::BoundaryIntegralOp<Real,KDIM0,KDIM1,UPSAMPLE,PDIM,RDIM>;
      ((BIOp*)self)[0](U, F);
    }

    void* biop;
    void* (*biop_build)(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, const sctl::Comm& comm);
    void (*biop_delete)(void** self);
    void (*biop_eval)(sctl::Vector<Real>& U, const sctl::Vector<Real>& F, void* self);
    sctl::Comm comm_;
};

template <class Real> VirtualCasing<Real>::VirtualCasing() : comm_(sctl::Comm::Self()), LaplaceFxdU(comm_), Svec(1), digits_(10), dosetup(true) {
}

template <class Real> void VirtualCasing<Real>::SetSurface(sctl::Long Nt, sctl::Integer Np, const sctl::Vector<Real>& X) {
  dosetup = true;
  SCTL_ASSERT(Nt*Np*COORD_DIM == X.Dim());
  Svec[0] = biest::Surface<Real>(Nt, Np);
  Svec[0].Coord() = X;
}

template <class Real> void VirtualCasing<Real>::SetAccuracy(sctl::Integer digits) {
  if (digits != digits_) dosetup = true;
  digits_ = digits;
}

template <class Real> void VirtualCasing<Real>::ComputeBext(sctl::Vector<Real>& Bext, const sctl::Vector<Real>& B) const {
  if (dosetup) {
    //BiotSavartFxU.SetupSingular(Svec, biest::BiotSavart3D<Real>::FxU(), digits_);
    LaplaceFxdU.SetupSingular(Svec, biest::Laplace3D<Real>::FxdU(), digits_);

    biest::SurfaceOp<Real> SurfOp(comm_, Svec[0].NTor(), Svec[0].NPol());
    SurfOp.Grad2D(dX, Svec[0].Coord());
    SurfOp.SurfNormalAreaElem(&normal, nullptr, dX, &Svec[0].Coord());
    dosetup = false;
  }

  sctl::Vector<Real> BdotN, J, Bext_;
  DotProd(BdotN, B, normal);
  CrossProd(J, normal, B);
  LaplaceFxdU.Eval(Bext_, BdotN);
  { //BiotSavartFxU.Eval(Bext, J);
    const sctl::Long N = J.Dim() / COORD_DIM;
    sctl::Vector<Real> gradG_J(N * COORD_DIM * COORD_DIM); gradG_J = 0;
    sctl::Vector<Real> gradG_J0(N*COORD_DIM, gradG_J.begin() + N*COORD_DIM*0, false);
    sctl::Vector<Real> gradG_J1(N*COORD_DIM, gradG_J.begin() + N*COORD_DIM*1, false);
    sctl::Vector<Real> gradG_J2(N*COORD_DIM, gradG_J.begin() + N*COORD_DIM*2, false);
    LaplaceFxdU.Eval(gradG_J0, sctl::Vector<Real>(N, J.begin() + N*0, false));
    LaplaceFxdU.Eval(gradG_J1, sctl::Vector<Real>(N, J.begin() + N*1, false));
    LaplaceFxdU.Eval(gradG_J2, sctl::Vector<Real>(N, J.begin() + N*2, false));

    if (Bext.Dim() != N * COORD_DIM) Bext.ReInit(N * COORD_DIM);
    for (sctl::Long i = 0; i < N; i++) {
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        const sctl::Integer k1 = (k+1)%COORD_DIM;
        const sctl::Integer k2 = (k+2)%COORD_DIM;
        Bext[k*N+i] = gradG_J[(k1*COORD_DIM+k2)*N+i] - gradG_J[(k2*COORD_DIM+k1)*N+i];
      }
    }
  }
  Bext += Bext_ + 0.5 * B;
}

template <class Real> void VirtualCasing<Real>::DotProd(sctl::Vector<Real>& AdotB, const sctl::Vector<Real>& A, const sctl::Vector<Real>& B) {
  sctl::Long N = A.Dim() / COORD_DIM;
  SCTL_ASSERT(A.Dim() == COORD_DIM * N);
  SCTL_ASSERT(B.Dim() == COORD_DIM * N);
  if (AdotB.Dim() != N) AdotB.ReInit(N);
  for (sctl::Long i = 0; i < N; i++) {
    Real AdotB_ = 0;
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      AdotB_ += A[k*N+i] * B[k*N+i];
    }
    AdotB[i] = AdotB_;
  }
};

template <class Real> void VirtualCasing<Real>::CrossProd(sctl::Vector<Real>& AcrossB, const sctl::Vector<Real>& A, const sctl::Vector<Real>& B) {
  sctl::Long N = A.Dim() / COORD_DIM;
  SCTL_ASSERT(A.Dim() == COORD_DIM * N);
  SCTL_ASSERT(B.Dim() == COORD_DIM * N);
  if (AcrossB.Dim() != COORD_DIM * N) AcrossB.ReInit(COORD_DIM * N);
  for (sctl::Long i = 0; i < N; i++) {
    sctl::StaticArray<Real,COORD_DIM> A_, B_, AcrossB_;
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      A_[k] = A[k*N+i];
      B_[k] = B[k*N+i];
    }
    AcrossB_[0] = A_[1] * B_[2] - B_[1] * A_[2];
    AcrossB_[1] = A_[2] * B_[0] - B_[2] * A_[0];
    AcrossB_[2] = A_[0] * B_[1] - B_[0] * A_[1];
    for (sctl::Integer k = 0; k < COORD_DIM; k++) {
      AcrossB[k*N+i] = AcrossB_[k];
    }
  }
};

template <class Real> void VirtualCasingTestData<Real>::SurfaceCoordinates(sctl::Vector<Real>& X, int Nt, int Np, biest::SurfType surf_type) {
  biest::Surface<Real> S(Nt,Np, surf_type);
  X = S.Coord();
}

template <class Real> void VirtualCasingTestData<Real>::BFieldData(sctl::Vector<Real>& Bext, sctl::Vector<Real>& Bint, int Nt, int Np, const sctl::Vector<Real>& X) {
  constexpr sctl::Integer COORD_DIM = 3;
  auto WriteVTK_ = [](std::string fname, const sctl::Vector<sctl::Vector<Real>>& coords, const sctl::Vector<sctl::Vector<Real>>& values) {
    biest::VTKData data;
    typedef biest::VTKData::VTKReal VTKReal;
    auto& point_coord =data.point_coord ;
    auto& point_value =data.point_value ;
    auto& line_connect=data.line_connect;
    auto& line_offset =data.line_offset ;
    constexpr sctl::Integer COORD_DIM = biest::VTKData::COORD_DIM;

    SCTL_ASSERT(coords.Dim() == values.Dim());
    for (sctl::Long j = 0; j < coords.Dim(); j++) { // set point_coord, line_connect
      const auto& coord = coords[j];
      const auto& value = values[j];
      sctl::Long N = coord.Dim() / COORD_DIM;
      sctl::Long dof = value.Dim() / N;
      SCTL_ASSERT(value.Dim() == dof * N);
      for (sctl::Long i = 0; i < N; i++) {
        line_connect.push_back(point_coord.size()/COORD_DIM);
        point_coord.push_back((VTKReal)coord[0*N+i]);
        point_coord.push_back((VTKReal)coord[1*N+i]);
        point_coord.push_back((VTKReal)coord[2*N+i]);
        for (sctl::Long k = 0; k < dof; k++) {
          point_value.push_back((VTKReal)value[k*N+i]);
        }
      }
      line_offset.push_back(line_connect.size());
    }
    data.WriteVTK(fname.c_str(), sctl::Comm::Self());
  };
  auto eval_BiotSavart = [](const sctl::Vector<Real>& Xt, const sctl::Vector<sctl::Vector<Real>>& source, const sctl::Vector<sctl::Vector<Real>>& density) {
    const auto& kernel = biest::BiotSavart3D<Real>::FxU();
    sctl::Long Nt = Xt.Dim() / COORD_DIM;
    SCTL_ASSERT(Xt.Dim() == COORD_DIM * Nt);
    SCTL_ASSERT(source.Dim() == density.Dim());

    sctl::Vector<Real> B(COORD_DIM*Nt);
    B = 0;
    for (sctl::Long i = 0; i < source.Dim(); i++) {
      const auto& Xs = source[i];
      const auto& Fs = density[i];
      sctl::Long Ns = Xs.Dim() / COORD_DIM;
      SCTL_ASSERT(Xs.Dim() == COORD_DIM * Ns);
      SCTL_ASSERT(Fs.Dim() == COORD_DIM * Ns);
      sctl::Vector<Real> SrcNormal(COORD_DIM*Ns);
      kernel(Xs,SrcNormal,Fs, Xt,B);
    }
    return B;
  };
  auto add_source_loop = [](sctl::Vector<sctl::Vector<Real>>& source, sctl::Vector<sctl::Vector<Real>>& density, std::initializer_list<Real> coord, std::initializer_list<Real> normal, Real radius) {
    auto cross_norm = [](const sctl::Vector<Real>& A, const sctl::Vector<Real>& B) {
      sctl::Vector<Real> C(COORD_DIM);
      C[0] = A[1]*B[2] - B[1]*A[2];
      C[1] = A[2]*B[0] - B[2]*A[0];
      C[2] = A[0]*B[1] - B[0]*A[1];
      Real r = sctl::sqrt<Real>(C[0]*C[0]+C[1]*C[1]+C[2]*C[2]);
      return C*(1/r);
    };
    sctl::Vector<Real> coord_(COORD_DIM), normal_(COORD_DIM), e0(COORD_DIM), e1(COORD_DIM);
    coord_[0] = coord.begin()[0];
    coord_[1] = coord.begin()[1];
    coord_[2] = coord.begin()[2];
    normal_[0] = normal.begin()[0];
    normal_[1] = normal.begin()[1];
    normal_[2] = normal.begin()[2];
    Real normal_scal = 1/sctl::sqrt<Real>(normal_[0]*normal_[0]+normal_[1]*normal_[1]+normal_[2]*normal_[2]);
    normal_ *= normal_scal;

    e0[0] = drand48();
    e0[1] = drand48();
    e0[2] = drand48();
    e0 = cross_norm(e0,normal_)*radius;
    e1 = cross_norm(e0,normal_)*radius;

    sctl::Long N = 10000;
    sctl::Vector<Real> X(COORD_DIM * N);
    sctl::Vector<Real> dX(COORD_DIM * N);
    for (sctl::Long i = 0; i < N; i++) {
      Real t = 2*sctl::const_pi<Real>()*i/N;
      sctl::Vector<Real> r = coord_ + e0*sctl::sin<Real>(t) + e1*sctl::cos<Real>(t);
      sctl::Vector<Real> dr = e0*sctl::cos<Real>(t) - e1*sctl::sin<Real>(t);
      X[0*N+i] = r[0];
      X[1*N+i] = r[1];
      X[2*N+i] = r[2];
      dX[0*N+i] = dr[0];
      dX[1*N+i] = dr[1];
      dX[2*N+i] = dr[2];
    }
    source.PushBack(X);
    density.PushBack(dX);
  };
  auto DotProd = [](sctl::Vector<Real>& AdotB, const sctl::Vector<Real>& A, const sctl::Vector<Real>& B) {
    sctl::Long N = A.Dim() / COORD_DIM;
    SCTL_ASSERT(A.Dim() == COORD_DIM * N);
    SCTL_ASSERT(B.Dim() == COORD_DIM * N);
    if (AdotB.Dim() != N) AdotB.ReInit(N);
    for (sctl::Long i = 0; i < N; i++) {
      Real AdotB_ = 0;
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        AdotB_ += A[k*N+i] * B[k*N+i];
      }
      AdotB[i] = AdotB_;
    }
  };

  sctl::Comm comm = sctl::Comm::Self();
  sctl::Vector<biest::Surface<Real>> Svec(1);
  Svec[0] = biest::Surface<Real>(Nt,Np);
  Svec[0].Coord() = X;

  sctl::Vector<sctl::Vector<Real>> source0, density0;
  sctl::Vector<sctl::Vector<Real>> source1, density1;
  { // Set inside sources (source0, density0)
    sctl::Long N = 10000;
    sctl::Vector<Real> X(COORD_DIM*N), dX(COORD_DIM*N);
    { // Set X, dX
      sctl::Long Nt = 100, Np = 100;
      sctl::Vector<Real> coord(COORD_DIM*Nt);
      { // Set coord
        auto S = biest::Surface<Real>(Nt,Np);
        biest::SurfaceOp<Real>::Upsample(Svec[0].Coord(), Svec[0].NTor(), Svec[0].NPol(), S.Coord(), Nt, Np);

        sctl::Vector<Real> normal, dX;
        biest::SurfaceOp<Real> SurfOp(comm, Nt, Np);
        SurfOp.Grad2D(dX, S.Coord());
        SurfOp.SurfNormalAreaElem(&normal, nullptr, dX, &S.Coord());
        S.Coord() += -2.17*normal;

        coord = 0;
        for (sctl::Long t = 0; t < Nt; t++) {
          for (sctl::Long p = 0; p < Np; p++) {
            coord[0*Nt+t] += S.Coord()[(0*Nt+t)*Np+p]/Np;
            coord[1*Nt+t] += S.Coord()[(1*Nt+t)*Np+p]/Np;
            coord[2*Nt+t] += S.Coord()[(2*Nt+t)*Np+p]/Np;
          }
        }
      }

      sctl::Vector<Real> dX_;
      biest::SurfaceOp<Real>::Upsample(coord,Nt,1, X,N,1);
      biest::SurfaceOp<Real> SurfOp(comm,N,1);
      SurfOp.Grad2D(dX_, X);
      for (sctl::Long i = 0; i < N; i++) {
        for (sctl::Integer k = 0; k < COORD_DIM; k++) {
          dX[k*N+i] = dX_[(2*k+0)*N+i];
        }
      }
    }
    source0.PushBack(X);
    density0.PushBack(dX*0.05);
  }
  { // Set outside sources (source1, density1)
    const sctl::Long N = source0[0].Dim()/COORD_DIM;
    const sctl::StaticArray<Real,COORD_DIM> X{source0[0][0*N],source0[0][1*N],source0[0][2*N]};
    const sctl::StaticArray<Real,COORD_DIM> Xn{density0[0][0*N],density0[0][1*N],density0[0][2*N]};
    const Real R = sctl::sqrt<Real>(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
    add_source_loop(source1, density1, {X[0],X[1],X[2]}, {Xn[0],Xn[1],Xn[2]}, R);
  }

  Bint = eval_BiotSavart(Svec[0].Coord(), source0, density0);
  Bext = eval_BiotSavart(Svec[0].Coord(), source1, density1);

  // Visualization
  //WriteVTK_("loop0", source0, density0);
  //WriteVTK_("loop1", source1, density1);
  //biest::WriteVTK("S", Svec, Bext+Bint);

  SCTL_UNUSED(WriteVTK_);
  SCTL_UNUSED(DotProd);
}
