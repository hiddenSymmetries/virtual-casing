#ifndef _BIEST_VIRTUAL_CASING_TEST_DATA_HPP_
#define _BIEST_VIRTUAL_CASING_TEST_DATA_HPP_

#include <biest.hpp>

template <class Real> void GenerateSurfaceCoordinates(sctl::Vector<Real>& X, long Nt, long Np, biest::SurfType surf_type = biest::SurfType::AxisymNarrow) {
  biest::Surface<Real> S(Nt,Np, surf_type);
  X = S.Coord();
}

template <class Real> void GenerateVirtualCasingTestData(sctl::Vector<Real>& Bext, sctl::Vector<Real>& Bint, long Nt, long Np, const sctl::Vector<Real>& X) {
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
  WriteVTK_("loop0", source0, density0);
  WriteVTK_("loop1", source1, density1);
  biest::WriteVTK("S", Svec, Bext+Bint);
  SCTL_UNUSED(WriteVTK_);
  SCTL_UNUSED(DotProd);
}

#endif
