
template <class Real> static void CompleteVecField(sctl::Vector<Real>& X, const sctl::Integer NFP, const sctl::Long Nt, const sctl::Long Np, const sctl::Vector<Real>& Y) {
  static constexpr sctl::Integer COORD_DIM = 3;
  SCTL_ASSERT(Y.Dim() == COORD_DIM*Nt*Np);
  if (X.Dim() != COORD_DIM*NFP*Nt*Np) X.ReInit(COORD_DIM*NFP*Nt*Np);
  for (sctl::Long k = 0; k < NFP; k++) {
    const Real cost = sctl::cos<Real>(2*sctl::const_pi<Real>()*k/NFP);
    const Real sint = sctl::sin<Real>(2*sctl::const_pi<Real>()*k/NFP);
    for (sctl::Long i = 0; i < Nt*Np; i++) {
      const Real x0 = Y[0*Nt*Np+i];
      const Real y0 = Y[1*Nt*Np+i];
      const Real z0 = Y[2*Nt*Np+i];

      const Real x = x0*cost - y0*sint;
      const Real y = x0*sint + y0*cost;
      const Real z = z0;

      X[(0*NFP+k)*Nt*Np+i] = x;
      X[(1*NFP+k)*Nt*Np+i] = y;
      X[(2*NFP+k)*Nt*Np+i] = z;
    }
  }
}

template <class Real> static void Resample(sctl::Vector<Real>& X1, const sctl::Long Nt1, const sctl::Long Np1, const sctl::Vector<Real>& X0, const sctl::Long Nt0, const sctl::Long Np0) {
  const sctl::Long skip_tor = (sctl::Long)std::ceil(Nt0/(Real)Nt1);
  const sctl::Long skip_pol = (sctl::Long)std::ceil(Np0/(Real)Np1);
  const sctl::Long dof = X0.Dim() / (Nt0 * Np0);
  SCTL_ASSERT(X0.Dim() == dof*Nt0*Np0);

  sctl::Vector<Real> XX;
  biest::SurfaceOp<Real>::Upsample(X0, Nt0, Np0, XX, Nt1*skip_tor, Np1*skip_pol);

  if (X1.Dim() != dof * Nt1*Np1) X1.ReInit(dof * Nt1*Np1);
  for (sctl::Long k = 0; k < dof; k++) {
    for (sctl::Long i = 0; i < Nt1; i++) {
      for (sctl::Long j = 0; j < Np1; j++) {
        X1[(k*Nt1+i)*Np1+j] = XX[((k*Nt1+i)*skip_tor*Np1+j)*skip_pol];
      }
    }
  }
}

template <class Real, sctl::Integer COORD_DIM, sctl::Integer KDIM0, sctl::Integer KDIM1> class BIOpWrapper {
  public:

    BIOpWrapper(const sctl::Comm& comm) : biop(nullptr), comm_(comm) {
      NFP_ = 0;
      trg_Nt_ = 0;
      trg_Np_ = 0;
      quad_Nt_ = 0;
      quad_Np_ = 0;
    }

    ~BIOpWrapper() {
      if (biop) biop_delete(&biop);
    }

    void SetupSingular(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, const sctl::Integer digits, const sctl::Integer NFP, const sctl::Long src_Nt, const sctl::Long src_Np, const sctl::Long trg_Nt, const sctl::Long trg_Np, const sctl::Long qNt = 0, const sctl::Long qNp = 0) {
      SCTL_ASSERT(Svec[0].NTor() % NFP == 0);
      trg_Nt_ = trg_Nt;
      trg_Np_ = trg_Np;
      NFP_ = NFP;

      Real optim_aspect_ratio = 0, cond = 1;
      sctl::StaticArray<Real,2> aspect_ratio{100,0};
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
          aspect_ratio[0] = std::min<Real>(aspect_ratio[0], sctl::sqrt<Real>(M[0][0]/M[1][1]));
          aspect_ratio[1] = std::max<Real>(aspect_ratio[1], sctl::sqrt<Real>(M[0][0]/M[1][1]));

          //M.SVD(U,S,Vt);
          //Real cond2 = std::max<Real>(S[0][0],S[1][1]) / std::min<Real>(S[0][0],S[1][1]);
          //cond = std::max<Real>(cond, sctl::sqrt<Real>(cond2));
        }
      }
      { // Set tor_upsample, cond
        optim_aspect_ratio = sctl::sqrt<Real>(aspect_ratio[0]*aspect_ratio[1]) * Svec[0].NTor() / Svec[0].NPol();
        cond = sctl::sqrt<Real>(aspect_ratio[1]/aspect_ratio[0]);
      }
      if (cond > 4) {
        SCTL_WARN("The surface mesh is highly anisotropic! Quadrature generation will be very slow. Consider using a better surface discretization.");
        std::cout<<"Mesh anisotropy = "<<cond<<'\n';
      }

      const sctl::Integer PDIM = digits*cond*1.6;
      if (qNt > 0 && qNp > 0) { // Set quad_Nt_, quad_Np_
        quad_Nt_ = qNt;
        quad_Np_ = qNp;
      } else {
        const sctl::Long surf_Nt = Svec[0].NTor();
        const sctl::Long surf_Np = Svec[0].NPol();

        quad_Np_ =       trg_Np_  * (sctl::Long)std::ceil(std::max<sctl::Long>(std::max<sctl::Long>(surf_Np, src_Np),                    2*PDIM+1) / (Real)      trg_Np_ );
        quad_Nt_ = (NFP_*trg_Nt_) * (sctl::Long)std::ceil(std::max<Real>((Real)std::max<sctl::Long>(surf_Nt, src_Nt), optim_aspect_ratio*quad_Np_) / (Real)(NFP_*trg_Nt_));

        for (sctl::Integer i = 0; i < 3; i++) { // adaptive refinement using double-layer test
          sctl::Long quad_Nt = (sctl::Long)std::ceil(quad_Nt_ / (Real)(NFP_*surf_Nt)) * (NFP_*surf_Nt);
          sctl::Long quad_Np = (sctl::Long)std::ceil(quad_Np_ / (Real)surf_Np) * surf_Np;

          BIOpWrapper<Real, COORD_DIM, 1, 1> dbl_op(comm_);
          dbl_op.SetupSingular(Svec, biest::Laplace3D<Real>::DxU(), digits, NFP_, 0, 0, surf_Nt/NFP_, surf_Np, quad_Nt, quad_Np);
          sctl::Vector<Real> U, F(quad_Nt * quad_Np); F = 1;
          dbl_op.Eval(U, F);

          Real err = 0;
          for (const auto& a : U) err = std::max<Real>(err, fabs(a-0.5));
          Real scal = std::max<Real>(1, (digits+1)/(std::log(err)/std::log((Real)0.1))); // assuming exponential/geometric convergence
          quad_Nt_ = scal * quad_Nt;
          quad_Np_ = scal * quad_Np;
          if (err < sctl::pow<Real>(0.1,digits) || scal < 1.5) break;
        }

        // quad_Nt_/quad_Np_ ~ optim_aspect_ratio
        quad_Np_ =       trg_Np_  * (sctl::Long)std::round((quad_Nt_/optim_aspect_ratio) / (Real)      trg_Np_ );
        quad_Nt_ = (NFP_*trg_Nt_) * (sctl::Long)std::round((optim_aspect_ratio*quad_Np_) / (Real)(NFP_*trg_Nt_));
      }
      if (!(quad_Np_ > 2*PDIM  ) ||
          !(quad_Nt_ >= Svec[0].NTor()) ||
          !(quad_Np_ >= Svec[0].NPol()) ||
          !(quad_Nt_ >= src_Nt ) ||
          !(quad_Np_ >= src_Np ) ||
          !((quad_Nt_/(NFP_*trg_Nt_))*(NFP_*trg_Nt_) == quad_Nt_) ||
          !((quad_Np_/      trg_Np_ )*      trg_Np_  == quad_Np_)) {
        std::cout<<"digits             = "<<digits            <<'\n';
        std::cout<<"cond               = "<<cond              <<'\n';
        std::cout<<"PDIM               = "<<PDIM              <<'\n';
        std::cout<<"NFP_               = "<<NFP_              <<'\n';
        std::cout<<"surf_Nt            = "<<Svec[0].NTor()    <<'\n';
        std::cout<<"surf_Np            = "<<Svec[0].NPol()    <<'\n';
        std::cout<<"src_Nt             = "<<src_Nt            <<'\n';
        std::cout<<"src_Np             = "<<src_Np            <<'\n';
        std::cout<<"trg_Nt_            = "<<trg_Nt_           <<'\n';
        std::cout<<"trg_Np_            = "<<trg_Np_           <<'\n';
        std::cout<<"quad_Nt_           = "<<quad_Nt_          <<'\n';
        std::cout<<"quad_Np_           = "<<quad_Np_          <<'\n';
        std::cout<<"optim_aspect_ratio = "<<optim_aspect_ratio<<'\n';
        SCTL_ASSERT(false);
      }
      { // Setup singular
        sctl::Vector<Real> XX;
        Resample(XX, quad_Nt_, quad_Np_, Svec[0].Coord(), Svec[0].NTor(), Svec[0].NPol());

        Svec_.ReInit(1);
        Svec_[0] = biest::Surface<Real>(quad_Nt_, quad_Np_);
        Svec_[0].Coord() = XX;

        sctl::Vector<sctl::Vector<sctl::Long>> trg_idx(1);
        trg_idx[0].ReInit(trg_Nt_*trg_Np_);
        const sctl::Long skip_Nt = quad_Nt_ / (NFP_ * trg_Nt_);
        const sctl::Long skip_Np = quad_Np_ / trg_Np_;
        for (sctl::Long i = 0; i < trg_Nt_; i++) {
          for (sctl::Long j = 0; j < trg_Np_; j++) {
            trg_idx[0][i*trg_Np_+j] = (i*skip_Nt*trg_Np_+j)*skip_Np;
            SCTL_ASSERT(trg_idx[0][i*trg_Np_+j] < quad_Nt_*quad_Np_);
          }
        }
        SetupSingular_(Svec_, ker, PDIM, trg_idx);
      }
    }

    void Eval(sctl::Vector<Real>& U, const sctl::Vector<Real>& F) const {
      SCTL_ASSERT(F.Dim() == KDIM0 * quad_Nt_ * quad_Np_);
      biop_eval(U, F, biop);
    }

    const sctl::Long QuadNt() const { return quad_Nt_; }
    const sctl::Long QuadNp() const { return quad_Np_; }

  private:

    template <sctl::Integer PDIM, sctl::Integer RDIM> static void* BIOpBuild(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, const sctl::Comm& comm, const sctl::Vector<sctl::Vector<sctl::Long>>& trg_idx) {
      using BIOp = biest::BoundaryIntegralOp<Real,KDIM0,KDIM1,1,PDIM,RDIM>;
      BIOp* biop = new BIOp(comm);
      biop[0].SetupSingular(Svec, ker, trg_idx);
      return biop;
    }
    template <sctl::Integer PDIM, sctl::Integer RDIM> static void BIOpDelete(void** self) {
      using BIOp = biest::BoundaryIntegralOp<Real,KDIM0,KDIM1,1,PDIM,RDIM>;
      delete (BIOp*)self[0];
      self[0] = nullptr;
    }
    template <sctl::Integer PDIM, sctl::Integer RDIM> static void BIOpEval(sctl::Vector<Real>& U, const sctl::Vector<Real>& F, void* self) {
      using BIOp = biest::BoundaryIntegralOp<Real,KDIM0,KDIM1,1,PDIM,RDIM>;
      ((BIOp*)self)[0](U, F);
    }

    template <sctl::Integer PDIM, sctl::Integer RDIM> void SetupSingular0(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, const sctl::Vector<sctl::Vector<sctl::Long>>& trg_idx) {
      if (biop) biop_delete(&biop);

      biop_build = BIOpBuild<PDIM,RDIM>;
      biop_delete = BIOpDelete<PDIM,RDIM>;
      biop_eval = BIOpEval<PDIM,RDIM>;

      biop = biop_build(Svec, ker, comm_, trg_idx);
    }
    void SetupSingular_(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, const sctl::Integer PDIM_, const sctl::Vector<sctl::Vector<sctl::Long>>& trg_idx) {
      SCTL_ASSERT(Svec.Dim() == 1);
      if (PDIM_ >= 64) {
        static constexpr sctl::Integer PDIM = 64, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 60) {
        static constexpr sctl::Integer PDIM = 60, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 56) {
        static constexpr sctl::Integer PDIM = 56, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 52) {
        static constexpr sctl::Integer PDIM = 52, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 48) {
        static constexpr sctl::Integer PDIM = 48, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 44) {
        static constexpr sctl::Integer PDIM = 44, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 40) {
        static constexpr sctl::Integer PDIM = 40, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 36) {
        static constexpr sctl::Integer PDIM = 36, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 32) {
        static constexpr sctl::Integer PDIM = 32, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 28) {
        static constexpr sctl::Integer PDIM = 28, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 24) {
        static constexpr sctl::Integer PDIM = 24, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 20) {
        static constexpr sctl::Integer PDIM = 20, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 16) {
        static constexpr sctl::Integer PDIM = 16, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >= 12) {
        static constexpr sctl::Integer PDIM = 12, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else if (PDIM_ >=  8) {
        static constexpr sctl::Integer PDIM =  8, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      } else {
        static constexpr sctl::Integer PDIM =  6, RDIM = PDIM*1.6;
        SetupSingular0<PDIM,RDIM>(Svec, ker, trg_idx);
      }
    }

    void* biop;
    void* (*biop_build)(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, const sctl::Comm& comm, const sctl::Vector<sctl::Vector<sctl::Long>>& trg_idx);
    void (*biop_delete)(void** self);
    void (*biop_eval)(sctl::Vector<Real>& U, const sctl::Vector<Real>& F, void* self);
    sctl::Integer NFP_;
    sctl::Long trg_Nt_, trg_Np_;
    sctl::Long quad_Nt_, quad_Np_;
    sctl::Vector<biest::Surface<Real>> Svec_;
    sctl::Comm comm_;
};

template <class Real> VirtualCasing<Real>::VirtualCasing() : comm_(sctl::Comm::Self()), LaplaceFxdU(comm_), Svec(1), NFP_(0), digits_(10), dosetup(true) {
}

template <class Real> void VirtualCasing<Real>::Setup(const sctl::Integer digits, const sctl::Integer NFP, const sctl::Long surf_Nt, const sctl::Long surf_Np, const std::vector<Real>& X, const sctl::Long src_Nt, const sctl::Long src_Np, const sctl::Long trg_Nt, const sctl::Long trg_Np) {
  dosetup = true;
  SCTL_ASSERT(surf_Nt*surf_Np*COORD_DIM == (sctl::Long)X.size());
  Svec[0] = biest::Surface<Real>(NFP*surf_Nt, surf_Np);
  CompleteVecField(Svec[0].Coord(), NFP, surf_Nt, surf_Np, sctl::Vector<Real>(X));

  NFP_ = NFP;
  src_Nt_ = src_Nt;
  src_Np_ = src_Np;
  trg_Nt_ = trg_Nt;
  trg_Np_ = trg_Np;
  digits_ = digits;
}

template <class Real> std::vector<Real> VirtualCasing<Real>::ComputeBext(const std::vector<Real>& B0) const {
  SCTL_ASSERT((sctl::Long)B0.size() == COORD_DIM * src_Nt_ * src_Np_);
  if (dosetup) {
    //BiotSavartFxU.SetupSingular(Svec, biest::BiotSavart3D<Real>::FxU(), digits_, NFP_, src_Nt_, src_Np_, trg_Nt_, trg_Np_);
    LaplaceFxdU.SetupSingular(Svec, biest::Laplace3D<Real>::FxdU(), digits_, NFP_, src_Nt_, src_Np_, trg_Nt_, trg_Np_);
    quad_Nt_ = LaplaceFxdU.QuadNt();
    quad_Np_ = LaplaceFxdU.QuadNp();

    sctl::Vector<Real> XX;
    Resample(XX, quad_Nt_, quad_Np_, Svec[0].Coord(), Svec[0].NTor(), Svec[0].NPol());

    biest::SurfaceOp<Real> SurfOp(comm_, quad_Nt_, quad_Np_);
    SurfOp.Grad2D(dX, XX);
    SurfOp.SurfNormalAreaElem(&normal, nullptr, dX, &XX);

    dosetup = false;
  }

  sctl::Vector<Real> B0_, B;
  CompleteVecField(B0_, NFP_, src_Nt_, src_Np_, sctl::Vector<Real>(B0));
  Resample(B, quad_Nt_, quad_Np_, B0_, NFP_*src_Nt_, src_Np_);

  std::vector<Real> Bext;
  { // Bext = BiotSavartFxU.Eval(normal x B);
    const sctl::Long N = trg_Nt_ * trg_Np_;
    sctl::Vector<Real> gradG_J(N * COORD_DIM * COORD_DIM); gradG_J = 0;
    sctl::Vector<Real> gradG_J0(N*COORD_DIM, gradG_J.begin() + N*COORD_DIM*0, false);
    sctl::Vector<Real> gradG_J1(N*COORD_DIM, gradG_J.begin() + N*COORD_DIM*1, false);
    sctl::Vector<Real> gradG_J2(N*COORD_DIM, gradG_J.begin() + N*COORD_DIM*2, false);

    sctl::Vector<Real> J;
    CrossProd(J, normal, B);
    LaplaceFxdU.Eval(gradG_J0, sctl::Vector<Real>(J.Dim()/3, J.begin() + (J.Dim()/3)*0, false));
    LaplaceFxdU.Eval(gradG_J1, sctl::Vector<Real>(J.Dim()/3, J.begin() + (J.Dim()/3)*1, false));
    LaplaceFxdU.Eval(gradG_J2, sctl::Vector<Real>(J.Dim()/3, J.begin() + (J.Dim()/3)*2, false));

    if ((sctl::Long)Bext.size() != N * COORD_DIM) Bext.resize(N * COORD_DIM);
    for (sctl::Long i = 0; i < N; i++) {
      for (sctl::Integer k = 0; k < COORD_DIM; k++) {
        const sctl::Integer k1 = (k+1)%COORD_DIM;
        const sctl::Integer k2 = (k+2)%COORD_DIM;
        Bext[k*N+i] = gradG_J[(k1*COORD_DIM+k2)*N+i] - gradG_J[(k2*COORD_DIM+k1)*N+i];
      }
    }
  }
  { // Bext += G[B . normal] + B/2
    sctl::Vector<Real> B_;
    Resample(B_, NFP_*trg_Nt_, trg_Np_, B0_, NFP_*src_Nt_, src_Np_);

    sctl::Vector<Real> BdotN, Bext_;
    DotProd(BdotN, B, normal);
    LaplaceFxdU.Eval(Bext_, BdotN);
    for (sctl::Long k = 0; k < COORD_DIM; k++) {
      for (sctl::Long i = 0; i < trg_Nt_; i++) {
        for (sctl::Long j = 0; j < trg_Np_; j++) {
          Bext[(k*trg_Nt_+i)*trg_Np_+j] += Bext_[(k*trg_Nt_+i)*trg_Np_+j] + 0.5 * B_[(k*NFP_*trg_Nt_+i)*trg_Np_+j];
        }
      }
    }
  }
  return Bext;
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

template <class Real> std::vector<Real> VirtualCasingTestData<Real>::SurfaceCoordinates(const sctl::Integer NFP, const sctl::Long Nt, const sctl::Long Np, const biest::SurfType surf_type) {
  biest::Surface<Real> S(NFP*Nt,Np, surf_type);
  const auto& X_ = S.Coord();
  std::vector<Real> X(COORD_DIM*Nt*Np);
  for (sctl::Long k = 0; k < COORD_DIM; k++) {
    for (sctl::Long i = 0; i < Nt*Np; i++) {
      X[k*Nt*Np+i] = X_[k*NFP*Nt*Np+i];
    }
  }
  return X;
}

template <class Real> std::tuple<std::vector<Real>, std::vector<Real>> VirtualCasingTestData<Real>::BFieldData(const sctl::Integer NFP, const sctl::Long surf_Nt, const sctl::Long surf_Np, const std::vector<Real>& X, const sctl::Long trg_Nt, const sctl::Long trg_Np) {
  constexpr static sctl::Integer COORD_DIM = 3;
  auto WriteVTK_ = [](const std::string& fname, const sctl::Vector<sctl::Vector<Real>>& coords, const sctl::Vector<sctl::Vector<Real>>& values) {
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
  auto add_source_loop = [](sctl::Vector<sctl::Vector<Real>>& source, sctl::Vector<sctl::Vector<Real>>& density, const std::initializer_list<Real> coord, const std::initializer_list<Real> normal, const Real radius) {
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
  sctl::Vector<Real> X_surf, X_trg;
  { // Set X_surf, X_trg
    CompleteVecField(X_surf, NFP, surf_Nt, surf_Np, sctl::Vector<Real>(X));
    Resample<Real>(X_trg, NFP*trg_Nt, trg_Np, X_surf, NFP*surf_Nt, surf_Np);
  }

  sctl::Vector<sctl::Vector<Real>> source0, density0;
  sctl::Vector<sctl::Vector<Real>> source1, density1;
  { // Set inside sources (source0, density0)
    sctl::Long N = 20000;
    sctl::Vector<Real> X(COORD_DIM*N), dX(COORD_DIM*N);
    { // Set X, dX
      sctl::Long Nt = 100, Np = 100;
      sctl::Vector<Real> coord(COORD_DIM*Nt);
      { // Set coord
        auto S = biest::Surface<Real>(Nt,Np);
        biest::SurfaceOp<Real>::Upsample(X_surf, NFP*surf_Nt, surf_Np, S.Coord(), Nt, Np);

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
  for (int i = 0; i < NFP; i++) { // Set outside sources (source1, density1)
    const sctl::Long N = source0[0].Dim()/COORD_DIM;
    const sctl::Long Nskip = i * N / NFP;

    const sctl::StaticArray<Real,COORD_DIM> X{source0[0][0*N+Nskip],source0[0][1*N+Nskip],source0[0][2*N+Nskip]};
    const sctl::StaticArray<Real,COORD_DIM> Xn{density0[0][0*N+Nskip],density0[0][1*N+Nskip],density0[0][2*N+Nskip]};
    const Real R = sctl::sqrt<Real>(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
    add_source_loop(source1, density1, {X[0],X[1],X[2]}, {Xn[0],Xn[1],Xn[2]}, R);
  }

  const auto Bint_ = eval_BiotSavart(X_trg, source0, density0);
  const auto Bext_ = eval_BiotSavart(X_trg, source1, density1);

  std::vector<Real> Bint(COORD_DIM*trg_Nt*trg_Np);
  std::vector<Real> Bext(COORD_DIM*trg_Nt*trg_Np);
  for (sctl::Integer k = 0; k < COORD_DIM; k++) { // Set Bint, Bext
    for (sctl::Long i = 0; i < trg_Nt*trg_Np; i++) {
      Bint[k*trg_Nt*trg_Np+i] = Bint_[k*NFP*trg_Nt*trg_Np+i];
      Bext[k*trg_Nt*trg_Np+i] = Bext_[k*NFP*trg_Nt*trg_Np+i];
    }
  }

  if (0) { // Visualization
    sctl::Vector<biest::Surface<Real>> Svec(1);
    Svec[0] = biest::Surface<Real>(NFP*trg_Nt,trg_Np);
    Svec[0].Coord() = X_trg;
    biest::WriteVTK("S", Svec, Bext_+Bint_);
    WriteVTK_("loop0", source0, density0);
    WriteVTK_("loop1", source1, density1);
  }

  SCTL_UNUSED(WriteVTK_);
  SCTL_UNUSED(DotProd);

  return std::make_tuple(std::move(Bext), std::move(Bint));
}
