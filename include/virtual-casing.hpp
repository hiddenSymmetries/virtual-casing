#include <biest.hpp>

#ifndef _BIEST_VIRTUAL_CASING_HPP_
#define _BIEST_VIRTUAL_CASING_HPP_

template <class Real, sctl::Integer COORD_DIM, sctl::Integer KDIM0, sctl::Integer KDIM1> class BIOpWrapper;

template <class Real> class VirtualCasing {
  static constexpr sctl::Integer COORD_DIM = 3;

  public:

    /**
     * Constructor
     */
    VirtualCasing();

    /**
     * Set surface coordinates with uniform discretization in toroidal and
     * poloidal directions.
     *
     * @param[in] Nt discretization order in toroidal direction.
     *
     * @param[in] Np discretization order in poloidal direction.
     *
     * @param[in] X the surface coordinates in the order {x11, x12, ..., x1Np,
     * x21, x22, ... , xNtNp, y11, ... , z11, ...}.
     */
    void SetSurface(sctl::Long Nt, sctl::Integer Np, const sctl::Vector<Real>& X);

    /**
     * Set approximate accuracy required.
     * @param[in] digits number of decimal digits of accuracy.
     */
    void SetAccuracy(sctl::Integer digits);

    /**
     * Recover the Bext component from the total field B = Bext + Bint by
     * applying the virtual-casing principle:
     * Bext = B/2 + gradG[B . n] + BiotSavart[n x B]
     *
     * @param[in] B the total magnetic field on the surface due to all currents.
     * B = {Bx11, Bx12, ..., Bx1Np, Bx21, Bx22, ... , BxNtNp, By11, ... , Bz11, ...},
     * where Nt and Np are the number of discretizations in toroidal and
     * poloidal directions.
     *
     * @param[out] Bext the component of magnetic field on the surface due to
     * currents in the exterior of the surface, computed using the
     * virtual-casing principle
     */
    void ComputeBext(sctl::Vector<Real>& Bext, const sctl::Vector<Real>& B) const;

  private:

    static void DotProd(sctl::Vector<Real>& AdotB, const sctl::Vector<Real>& A, const sctl::Vector<Real>& B);

    static void CrossProd(sctl::Vector<Real>& AcrossB, const sctl::Vector<Real>& A, const sctl::Vector<Real>& B);

    sctl::Comm comm_;
    mutable BIOpWrapper<Real,COORD_DIM,3,3> BiotSavartFxU;
    mutable BIOpWrapper<Real,COORD_DIM,1,3> LaplaceFxdU;
    sctl::Vector<biest::Surface<Real>> Svec;
    sctl::Integer digits_;
    mutable sctl::Vector<Real> dX, normal;
    mutable bool dosetup;
};

#include <virtual-casing.txx>

#endif
