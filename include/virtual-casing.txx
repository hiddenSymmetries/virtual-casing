
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
      if (cond > 5) {
        SCTL_WARN("The surface mesh is highly anisotropic! Quadrature generation will be very slow. Consider using a better surface discretization.");
        std::cout<<"Mesh aspect-ratio = "<<cond<<'\n';
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
    template <sctl::Integer UPSAMPLE, sctl::Integer PDIM> void SetupSingular1(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, sctl::Integer RDIM_) {
      if (RDIM_ >= 100) {
        constexpr sctl::Integer RDIM = 100;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      } else if (RDIM_ >= 90) {
        constexpr sctl::Integer RDIM = 90;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      } else if (RDIM_ >= 80) {
        constexpr sctl::Integer RDIM = 80;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      } else if (RDIM_ >= 70) {
        constexpr sctl::Integer RDIM = 70;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      } else if (RDIM_ >= 60) {
        constexpr sctl::Integer RDIM = 60;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      } else if (RDIM_ >= 50) {
        constexpr sctl::Integer RDIM = 50;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      } else if (RDIM_ >= 40) {
        constexpr sctl::Integer RDIM = 40;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      } else if (RDIM_ >= 30) {
        constexpr sctl::Integer RDIM = 30;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      } else if (RDIM_ >= 20) {
        constexpr sctl::Integer RDIM = 20;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      } else if (RDIM_ >= 10) {
        constexpr sctl::Integer RDIM = 10;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      } else {
        constexpr sctl::Integer RDIM = 1;
        SetupSingular0<UPSAMPLE, PDIM, RDIM>(Svec, ker);
      }
    }
    template <sctl::Integer UPSAMPLE> void SetupSingular2(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, sctl::Integer PDIM_, sctl::Integer RDIM_) {
      sctl::Integer SurfDim = std::min(Svec[0].NPol(), Svec[0].NTor());
      if (PDIM_ >= 64 && SurfDim > 64*2+1) {
        static constexpr sctl::Integer PDIM = 64;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 60 && SurfDim > 60*2+1) {
        static constexpr sctl::Integer PDIM = 60;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 56 && SurfDim > 56*2+1) {
        static constexpr sctl::Integer PDIM = 56;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 52 && SurfDim > 52*2+1) {
        static constexpr sctl::Integer PDIM = 52;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 48 && SurfDim > 48*2+1) {
        static constexpr sctl::Integer PDIM = 48;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 44 && SurfDim > 44*2+1) {
        static constexpr sctl::Integer PDIM = 44;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 40 && SurfDim > 40*2+1) {
        static constexpr sctl::Integer PDIM = 40;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 36 && SurfDim > 36*2+1) {
        static constexpr sctl::Integer PDIM = 36;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 32 && SurfDim > 32*2+1) {
        static constexpr sctl::Integer PDIM = 32;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 28 && SurfDim > 28*2+1) {
        static constexpr sctl::Integer PDIM = 28;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 24 && SurfDim > 24*2+1) {
        static constexpr sctl::Integer PDIM = 24;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 20 && SurfDim > 20*2+1) {
        static constexpr sctl::Integer PDIM = 20;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 16 && SurfDim > 16*2+1) {
        static constexpr sctl::Integer PDIM = 16;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >= 12 && SurfDim > 12*2+1) {
        static constexpr sctl::Integer PDIM = 12;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else if (PDIM_ >=  8 && SurfDim >  8*2+1) {
        static constexpr sctl::Integer PDIM =  8;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      } else {
        static constexpr sctl::Integer PDIM =  6;
        SetupSingular1<UPSAMPLE,PDIM>(Svec, ker, RDIM_);
      }
    }
    void SetupSingular3(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, sctl::Integer UPSAMPLE_, sctl::Integer PDIM_, sctl::Integer RDIM_) {
      if (RDIM_ >= 10) {
        SetupSingular2<10>(Svec, ker, PDIM_, RDIM_);
      } else if (RDIM_ >= 9) {
        SetupSingular2< 9>(Svec, ker, PDIM_, RDIM_);
      } else if (RDIM_ >= 8) {
        SetupSingular2< 8>(Svec, ker, PDIM_, RDIM_);
      } else if (RDIM_ >= 7) {
        SetupSingular2< 7>(Svec, ker, PDIM_, RDIM_);
      } else if (RDIM_ >= 6) {
        SetupSingular2< 6>(Svec, ker, PDIM_, RDIM_);
      } else if (RDIM_ >= 5) {
        SetupSingular2< 5>(Svec, ker, PDIM_, RDIM_);
      } else if (RDIM_ >= 4) {
        SetupSingular2< 4>(Svec, ker, PDIM_, RDIM_);
      } else if (RDIM_ >= 3) {
        SetupSingular2< 3>(Svec, ker, PDIM_, RDIM_);
      } else if (RDIM_ >= 2) {
        SetupSingular2< 2>(Svec, ker, PDIM_, RDIM_);
      } else {
        SetupSingular2< 1>(Svec, ker, PDIM_, RDIM_);
      }
    }

    template <sctl::Integer UPSAMPLE> void SetupSingular_(const sctl::Vector<biest::Surface<Real>>& Svec, const biest::KernelFunction<Real,COORD_DIM,KDIM0,KDIM1>& ker, sctl::Integer PDIM_) {
      sctl::Integer SurfDim = std::min(Svec[0].NPol(), Svec[0].NTor());
      if (PDIM_ >= 64 && SurfDim > 64*2+1) {
        static constexpr sctl::Integer PDIM = 64, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 60 && SurfDim > 60*2+1) {
        static constexpr sctl::Integer PDIM = 60, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 56 && SurfDim > 56*2+1) {
        static constexpr sctl::Integer PDIM = 56, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 52 && SurfDim > 52*2+1) {
        static constexpr sctl::Integer PDIM = 52, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 48 && SurfDim > 48*2+1) {
        static constexpr sctl::Integer PDIM = 48, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 44 && SurfDim > 44*2+1) {
        static constexpr sctl::Integer PDIM = 44, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 40 && SurfDim > 40*2+1) {
        static constexpr sctl::Integer PDIM = 40, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 36 && SurfDim > 36*2+1) {
        static constexpr sctl::Integer PDIM = 36, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 32 && SurfDim > 32*2+1) {
        static constexpr sctl::Integer PDIM = 32, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 28 && SurfDim > 28*2+1) {
        static constexpr sctl::Integer PDIM = 28, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 24 && SurfDim > 24*2+1) {
        static constexpr sctl::Integer PDIM = 24, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 20 && SurfDim > 20*2+1) {
        static constexpr sctl::Integer PDIM = 20, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 16 && SurfDim > 16*2+1) {
        static constexpr sctl::Integer PDIM = 16, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >= 12 && SurfDim > 12*2+1) {
        static constexpr sctl::Integer PDIM = 12, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else if (PDIM_ >=  8 && SurfDim >  8*2+1) {
        static constexpr sctl::Integer PDIM =  8, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
      } else {
        static constexpr sctl::Integer PDIM =  6, RDIM = PDIM*1.5;
        SetupSingular0<UPSAMPLE,PDIM,RDIM>(Svec, ker);
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

template <class Real> VirtualCasing<Real>::VirtualCasing() : comm_(sctl::Comm::Self()), BiotSavartFxU(comm_), LaplaceFxdU(comm_), Svec(1), digits_(10), dosetup(true) {
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
    BiotSavartFxU.SetupSingular(Svec, biest::BiotSavart3D<Real>::FxU(), digits_);
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
  BiotSavartFxU.Eval(Bext, J);
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

