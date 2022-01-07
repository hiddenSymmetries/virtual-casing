

      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: novers(:),ixyzso(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      real *8, allocatable :: hvec(:,:)
      real *8, allocatable :: rhs(:)
      real *8, allocatable :: hvecdiv(:)

      real *8, allocatable :: ynm(:),unm(:,:),xnm(:,:)
      complex *16, allocatable :: zynm(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8, allocatable :: rjvec(:,:),rho(:)
      real *8, allocatable :: curlj(:,:),gradrho(:,:)

      real *8, allocatable :: wnear(:)
      real *8, allocatable :: w(:,:)
      real *8 vtmp1(3),vtmp2(3)
      real *8 xyz_out(3),xyz_in(3)
      complex *16 zpars


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      igeomtype = 1

      ipars(1) = 2

      npatches = 12*(4**ipars(1)) 
      norder = 7 
      npols = (norder+1)*(norder+2)/2
      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 0


      call setup_geom(igeomtype,norder,npatches,ipars, &
            srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      ra = 0
      do i=1,npts
        ra = ra + wts(i)
      enddo
      call prin2('surface area of sphere=*',ra,1)

      nn = 3
      mm = 2
      nmax = nn
      allocate(w(0:nmax,0:nmax))
      allocate(zynm(npts),ynm(npts),unm(3,npts),xnm(3,npts))
      call l3getsph(nmax,mm,nn,12,srcvals,zynm,npts,w)

      do i=1,npts
        ynm(i) = real(zynm(i))
      enddo


      call surf_grad(npatches,norders,ixyzs,iptype,npts,srccoefs, &
       srcvals,ynm,unm)
      
      do i=1,npts
        unm(1:3,i) = unm(1:3,i)/sqrt(nn*(nn+1.0d0))
        call cross_prod3d(srcvals(10,i),unm(1,i),xnm(1,i))
      enddo
      allocate(hvecdiv(npts))

      call surf_div(npatches,norders,ixyzs,iptype,npts,srccoefs, &
       srcvals,xnm,hvecdiv)
      
      ra = 0
      do i=1,npts
        ra = ra + hvecdiv(i)**2*wts(i)
      enddo
      ra = sqrt(ra)
      call prin2('error in surface divergence of xnm=*',ra,1)

      allocate(rjvec(3,npts),rho(npts),curlj(3,npts),gradrho(3,npts))


      rcu = 0.0d0
      rcx = 1.1d0
      do i=1,npts
        rjvec(1:3,i) = rcu*unm(1:3,i) + rcx*xnm(1:3,i)
        rho(i) = ynm(i)
      enddo

      eps = 1.0d-5

      allocate(ipatch_id(npts),uvs_targ(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, &
       ipatch_id,uvs_targ)

      call lpcomp_virtualcasing(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,12,npts,srcvals,ipatch_id,uvs_targ,&
       eps,rjvec,rho,curlj,gradrho)


      errnc = 0
      errnd = 0
      rnc = 0
      rnd = 0

      runc = (nn+1.0d0)/(2*nn+1.0d0)*rcu
      rxnc = (nn+0.0d0)/(2*nn+1.0d0)*rcx

      call prin2('runc=*',runc,1)
      call prin2('rxnc=*',rxnc,1)
      call prin2('rcu=*',rcu,1)
      call prin2('rcx=*',rcx,1)
      rund = 0
      rxnd = -sqrt((nn+0.0d0)*(nn+1.0d0))/(2*nn+1.0d0)*rcx
      do i=1,npts
        call cross_prod3d(srcvals(10,i),curlj(1,i),vtmp1)
        vtmp1(1:3) = vtmp1(1:3) + rjvec(1:3,i)/2

        vtmp2(1) = runc*unm(1,i) + rxnc*xnm(1,i) 
        vtmp2(2) = runc*unm(2,i) + rxnc*xnm(2,i) 
        vtmp2(3) = runc*unm(3,i) + rxnc*xnm(3,i) 
        

        rnc = rnc + (vtmp2(1)**2 + vtmp2(2)**2 + vtmp2(3)**2)*wts(i)
        errnc = errnc + (vtmp1(1)-vtmp2(1))**2*wts(i)
        errnc = errnc + (vtmp1(2)-vtmp2(2))**2*wts(i)
        errnc = errnc + (vtmp1(3)-vtmp2(3))**2*wts(i)
        wtmp = 0
        call dot_prod3d(srcvals(10,i),curlj(1,i),wtmp)

        rnd = rnd + ynm(i)**2*wts(i)
        errnd = errnd + (wtmp - rxnd*ynm(i))**2*wts(i)
        if(i.lt.3) print *, vtmp1(1),vtmp2(1),vtmp1(1)/vtmp2(1)
        if(i.lt.3) print *, wtmp,rxnd*ynm(i),rxnd*ynm(i)/wtmp
      enddo

      errnd = sqrt(errnd/rnd)
      errnc = sqrt(errnc/rnc)


      call prin2('rnc=*',rnc,1)
      call prin2('rnd=*',rnd,1)

      call prin2('error in n times curl s0 = *',errnc,1)
      call prin2('error in n dot curl s0 =*',errnd,1)


      runc = 0
      rxnc = sqrt((nn+0.0d0)*(nn+1.0d0))/(2*nn+1.0d0)

      rynd = -(nn+1.0d0)/(2*nn+1.0d0)

      errnc = 0
      errnd = 0
      rnc = 0
      rnd = 0
      do i=1,npts
        call cross_prod3d(srcvals(10,i),gradrho(1,i),vtmp1)
        vtmp2(1:3) = rxnc*xnm(1:3,i)
        rnc = rnc + (vtmp2(1)**2+vtmp2(2)**2+vtmp2(3)**2)*wts(i)
        errnc = errnc + (vtmp2(1)-vtmp1(1))**2*wts(i)
        errnc = errnc + (vtmp2(2)-vtmp1(2))**2*wts(i)
        errnc = errnc + (vtmp2(3)-vtmp1(3))**2*wts(i)
        
        wtmp = 0
        call dot_prod3d(srcvals(10,i),gradrho(1,i),wtmp)
        wtmp = wtmp - rho(i)/2
        errnd = errnd + (wtmp-rynd*ynm(i))**2*wts(i)
        rnd = rnd + ynm(i)**2*wts(i)
      enddo

      errnc = sqrt(errnc/rnc)
      errnd = sqrt(errnd/rnd)
      call prin2('error in n times grad s0 = *',errnc,1)
      call prin2('error in n dot grad s0 =*',errnd,1)




      return
      end





      subroutine setup_geom(igeomtype,norder,npatches,ipars, & 
         srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval,xtri_wtorus_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, & 
           triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, & 
              ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4, &
          npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2), &
          nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,  &
              iptr3,iptr4, norder, &
              'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4, &
          npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.3) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 0
        vmax = 2*pi
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2), &
          nover,npatches,npatches,triaskel)
        call prinf('npatches=*',npatches,1)
         
        p1(1) = 1
        p1(2) = 2
        p1(3) = 0.25d0

        p2(1) = 1.2d0
        p2(2) = 1.0d0
        p2(3) = 1.7d0

!
!         numberof oscillations
!
        p4(1) = 5.0d0


        ptr1 => triaskel(1,1,1)
        ptr2 => p1(1)
        ptr3 => p2(1)
        ptr4 => p4(1)
        xtri_geometry => xtri_wtorus_eval
        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, &
              ptr3,ptr4, norder,&
              'Triangulated surface of the wtorus')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4, &
          npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      if(igeomtype.eq.4) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 0
        vmax = 2*pi
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2), &
          nover,npatches,npatches,triaskel)
        call prinf('npatches=*',npatches,1)
         
        p1(1) = 1.0d0
        p1(2) = 3.0d0
        p1(3) = 0.25d0

        p2(1) = 1.0d0
        p2(2) = 1.0d0
        p2(3) = 1.0d0

!
!         number of oscillations
!
        p4(1) = 0.0d0


        ptr1 => triaskel(1,1,1)
        ptr2 => p1(1)
        ptr3 => p2(1)
        ptr4 => p4(1)
        xtri_geometry => xtri_wtorus_eval
        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, &
              ptr3,ptr4, norder, &
              'Triangulated surface of the torus')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4, &
          npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end





      subroutine test_exterior_pt(npatches,norder,npts,srcvals, &
        srccoefs,wts,xyzout,isout)
!
!
!  this subroutine tests whether the pt xyzin, is
!  in the exterior of a surface, and also estimates the error
!  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
!  centered at the interior point. Whether a point 
!  is in the interior or not is tested using Gauss' 
!  identity for the flux due to a point charge
!
!
!  input:
!    npatches - integer
!       number of patches
!    norder - integer
!       order of discretization
!    npts - integer
!       total number of discretization points on the surface
!    srccoefs - real *8 (9,npts)
!       koornwinder expansion coefficients of geometry info
!    xyzout -  real *8 (3)
!       point to be tested
!
!  output: 
!    isout - boolean
!      whether the target is in the interior or not
!

      implicit none
      integer npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
      real *8 ra,ds
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2


      zk = 0

      ra = 0



      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          call h3d_sprime(xyzout,12,srcvals(1,i),0,dpars,1,zk,0,ipars, &
            val)

          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-1) isout = .false.
      if(abs(ra).le.1.0d-1) isout = .true.

      return
      end

   






      subroutine l3getsph(nmax,mm,nn,ndx,xyzs,ynms,npts,ynm)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      complex *16 ynms(npts),ima
      real *8 rat1(10000),rat2(10000)
      real *8 ynm(0:nmax,0:nmax)
      data ima/(0.0d0,1.0d0)/
  
      call ylgndrini(nmax, rat1, rat2)
  
      do i=1,npts
        x=xyzs(1,i)
        y=xyzs(2,i)
        z=xyzs(3,i)
        r=sqrt(x**2+y**2+z**2)
        call cart2polar(xyzs(1,i),r,theta,phi)
        ctheta = cos(theta)
        call ylgndrf(nmax, ctheta, ynm, rat1, rat2)
        ynms(i) = ynm(nn,abs(mm))*exp(ima*mm*phi)        
      enddo
       
      return
      end



