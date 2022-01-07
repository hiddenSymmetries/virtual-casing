      subroutine getnearquad_magnetostatics(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad, &
        wnear)
!
!  This subroutine generates the near field quadrature
!  for the representations:
!
!  wnear: 
!   \nabla_{1} S_{0}, \nabla_{2} S_{0}, \nabla_{3} S_{0},
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a patch centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!  
!  NOTES:
!    - wnear must be of size nquad,3 as 3 different layer
!      potentials are returned
!      * the first kernel is \nabla_{x} S_{0}
!      * the second kernel is \nabla_{y} S_{0}
!      * the third kernel is \nabla_{z} S_{0}
! 
!  Input arguments:
! 
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - ndtarg: integer
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: double precision(ndtarg,ntarg)
!        target info, the first three components
!        must be xyz coordinates
!    - ipatch_id: integer(ntarg)
!        id of patch of target i, id=-1, if target is off-surface
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates of patch if target on surface,
!        otherwise not used
!    - eps: real *8
!        precision requested
!    - dpars: real *8 (1)
!        the parameter k for the yukawa kernels
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(ntarg+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: integer
!        radius parameter for near field
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 4*nquad since there are 4 kernels
!        per source target pair
!
!  Output arguments
!    - wnear: complex *16(3*nquad)
!        The desired near field quadrature
!               
!
      implicit none
      integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      integer, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg),uvs_targ(2,ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: eps
      integer, intent(in) :: iquadtype,nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer, intent(in) :: nquad
      real *8, intent(out) :: wnear(nquad,3)

      integer ndz,ndd,ndi
      integer ipv
      real *8 dpars
      integer ipars
      complex *16 zpars

      procedure (), pointer :: fker

      external l3d_sgradx,l3d_sgrady,l3d_sgradz

      ndd = 0
      ndz = 0
      ndi = 0
      ipv = 1
      dpars = 0
      zpars = 0
      ipars = 0

      fker => l3d_sgradx
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,1))

      fker => l3d_sgrady
      ipv = 1
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,2))

      fker => l3d_sgradz
      ipv = 1
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,3))



      return
      end subroutine getnearquad_magnetostatics 
!
!
!
!
!
!
      subroutine lpcomp_virtualcasing(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,rjvec,rho,curlj,gradrho)

!
!  This subroutine evaluates the sets of potentials to the following
!  representations where j is assumed to be a tangential vector field, and rho
!  a scalar function defined on the surface:
!  1. curl (S_{0}[j]) 
!  2. grad S_{0}[rho]
!
!
!  The identity term is not included in the representation if targets are on 
!  surface
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Input arguments:
! 
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - ndtarg: integer
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: double precision(ndtarg,ntarg)
!        target info, the first three components
!        must be xyz coordinates
!    - ipatch_id: integer(ntarg)
!        id of patch of target i, id=-1, if target is off-surface
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates of patch if target on surface,
!        otherwise not used
!    - eps: real *8
!        precision requested
!    - rjvec: real *8(3,npts)
!        The current density j
!    - rho: real *8 (npts)
!        The charge density rho
!
!  Output arguments:
!    - curlj: real *8 (3,npts)
!         returns the potential curl S_{0}[j]
!    - gradrho: real *8 (3,npts)
!         returns the potential grad S_{0}[\rho]
!
!
!
!
      implicit none
      integer, intent(in) :: npatches,norders(npatches)
      integer, intent(in) :: iptype(npatches),ixyzs(npatches+1)
      integer, intent(in) :: npts,ndtarg,ntarg
      real *8, intent(in) :: srcvals(12,npts),srccoefs(9,npts)
      real *8, intent(in) :: targs(ndtarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      real *8, intent(in) :: rjvec(3,npts),rho(npts)
      real *8, intent(out) :: curlj(3,ntarg),gradrho(3,ntarg)

      integer nptso,nnz,nquad
      integer nover,npolso,npols,norder
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:,:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer iptype_avg,norder_avg,iquadtype,npts_over,ikerorder
      integer i
      real *8 rfac,rfac0
      complex *16 zpars

      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)
! Get centroid and bounding sphere info

      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
          srccoefs,cms,rads)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO

!
!  Find near quadrature corrections
!
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,ntarg,nnz)

      allocate(row_ptr(ntarg+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,row_ptr, & 
             col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind, &
              iquad)

!
! Estimate oversampling required for far-field, and oversample geometry
!
      ikerorder = 0

      allocate(novers(npatches),ixyzso(npatches+1))

      zpars = 0.0d0

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms, &
         rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,zpars, &
         nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over, &
             srcover,wover)
!
!  Compute near quadrature corrections
!

      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad,3))
      wnear = 0
      iquadtype = 1
      call getnearquad_magnetostatics(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad, &
        wnear)
!
!
!  compute layer potential
!
      call lpcomp_virtualcasing_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs, &
        eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,rjvec,rho, &
        novers,npts_over,ixyzso,srcover,wover,curlj,gradrho)
        

      return
      end subroutine lpcomp_virtualcasing
!
!
!
!
!
!
      subroutine lpcomp_virtualcasing_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs, &
        eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,rjvec,rho, &
        novers,nptso,ixyzso,srcover,whtsover,curlj,gradrho)

!
!  This subroutine evaluates the sets of potentials to the following
!  representations where j is assumed to be a tangential vector field, and rho
!  a scalar function defined on the surface:
!  1. curl (S_{0}[j]) 
!  2. grad S_{0}[rho]
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  The identity term is not included in teh representation
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
!
!  Input arguments:
! 
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - ndtarg: integer
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: double precision(ndtarg,ntarg)
!        target info, the first three components
!        must be xyz coordinates
!    - eps: real *8
!        precision requested
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 4*nquad since there are 4 kernels
!        per source target pair
!    - wnear: real *8(3*nquad)
!        Precomputed near field quadrature
!          * the first kernel is \nabla_{x} S_{0}
!          * the second kernel is \nabla_{y} S_{0}
!          * the third kernel is \nabla_{z} S_{0}
!    - rjvec: real *8(3,npts)
!        The current density j
!    - rho: real *8 (npts)
!        The charge density rho
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!
!  Output arguments:
!    - curlj: real *8 (3,ntarg)
!         returns the potential curl S_{0}[j]
!    - gradrho: real *8 (3,ntarg)
!         returns the potential grad S_{0}[\rho]
!

      implicit none
      integer npatches,norder,npols,npts
      integer ndtarg,ntarg
      integer norders(npatches),ixyzs(npatches+1)
      integer ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg) 
      integer nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      real *8 rjvec(3,npts),rho(npts)
      real *8 wnear(nquad,3)

      integer novers(npatches)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      real *8 curlj(3,ntarg),gradrho(3,ntarg)
      real *8, allocatable :: wts(:)

      real *8 rhom,rhop,rmum,uf,vf,wtmp
      real *8 u1,u2,u3,u4,w1,w2,w3,w4,w5

      real *8, allocatable :: sources(:,:),targtmp(:,:)
      real *8, allocatable :: charges0(:,:)
      real *8, allocatable :: sigmaover(:,:),abc0(:,:)
      real *8, allocatable :: pot_aux(:,:),grad_aux(:,:,:)
      real *8, allocatable :: pcurltmp(:,:)
      real *8, allocatable :: dpottmp(:),dgradtmp(:,:)
      real *8 vtmp1(3),vtmp2(3),vtmp3(3),rncj,errncj

      integer ns,nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(4)

      integer i,j,jpatch,jquadstart,jstart,count1,count2
      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8, allocatable :: ctmp0(:,:)
      real *8 thresh,ra,erra
      real *8 rr,rmin
      real *8 over4pi
      real *8 rbl(3),rbm(3)
      integer nss,ii,l,npover
      complex *16 ima,ztmp

      integer nd,ntarg0,nmax
      integer ndd,ndz,ndi,ier

      real *8 ttot,done,pi
      data ima/(0.0d0,1.0d0)/
      data over4pi/0.07957747154594767d0/

      parameter (ntarg0=1)

      ns = nptso
      ntarg = npts
      done = 1
      pi = atan(done)*4

      ifpgh = 0
      ifpghtarg = 2
      allocate(sources(3,ns),targtmp(3,npts))

      nmax = 0

!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,ixyzso,&
        nmax)

!
!  Allocate various densities
!

      allocate(sigmaover(4,ns),abc0(4,npts))

!
!  extract source and target info

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
!$OMP END PARALLEL DO      


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        targtmp(1,i) = targs(1,i)
        targtmp(2,i) = targs(2,i)
        targtmp(3,i) = targs(3,i)
      enddo
!$OMP END PARALLEL DO


      nd = 4
      allocate(charges0(nd,ns))

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npts
        abc0(1,i) = rjvec(1,i)
        abc0(2,i) = rjvec(2,i) 
        abc0(3,i) = rjvec(3,i)
        abc0(4,i) = rho(i)
      enddo
!$OMP END PARALLEL DO

      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
          npts,abc0,novers,ixyzso,ns,sigmaover)
        
!
!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,ns
        charges0(1:4,i) = sigmaover(1:4,i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      

      allocate(pot_aux(nd,ntarg),grad_aux(nd,3,ntarg))

!      print *, "before fmm"

      call lfmm3d_t_c_g_vec(nd,eps,ns,sources,charges0,npts,targtmp, &
        pot_aux,grad_aux,ier)

!$OMP PARALLEL DO DEFAULT(SHARED)         
      do i=1,ntarg
        curlj(1,i) = grad_aux(3,2,i) - grad_aux(2,3,i)
        curlj(2,i) = grad_aux(1,3,i) - grad_aux(3,1,i) 
        curlj(3,i) = grad_aux(2,1,i) - grad_aux(1,2,i)
        gradrho(1:3,i) = grad_aux(4,1:3,i)
      enddo
!$OMP END PARALLEL DO      

!      print *, "after fmm"

!
!  Add near quadrature correction
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,w1,w2,w3)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            w1 = wnear(jquadstart+l-1,1)
            w2 = wnear(jquadstart+l-1,2)
            w3 = wnear(jquadstart+l-1,3)
            curlj(1,i) = curlj(1,i) + w2*abc0(3,jstart+l-1) - &
              w3*abc0(2,jstart+l-1)
            curlj(2,i) = curlj(2,i) + w3*abc0(1,jstart+l-1) - &
              w1*abc0(3,jstart+l-1)
            curlj(3,i) = curlj(3,i) + w1*abc0(2,jstart+l-1) - &
              w2*abc0(1,jstart+l-1)
            gradrho(1,i) = gradrho(1,i) + w1*abc0(4,jstart+l-1)
            gradrho(2,i) = gradrho(2,i) + w2*abc0(4,jstart+l-1)
            gradrho(3,i) = gradrho(3,i) + w3*abc0(4,jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO     
      
!      print *, "after Laplace near correction"
!      print *, "nmax=",nmax
!      print *, "nd=",nd

      call get_fmm_thresh(12,ns,srcover,12,npts,srcvals,thresh)

!      print *, "Thresh=",thresh


!
! Subtract near contributions computed via fmm
!
      allocate(dpottmp(nd),dgradtmp(nd,3))
      allocate(ctmp0(nd,nmax),srctmp2(3,nmax))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,l,jstart,nss,dpottmp,dgradtmp)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            ctmp0(1:4,nss)=charges0(1:4,l)
          enddo
        enddo

        dpottmp = 0
        dgradtmp = 0

        call l3ddirectcg(nd,srctmp2,ctmp0,nss,targtmp(1,i), &
          ntarg0,dpottmp,dgradtmp,thresh)
        curlj(1,i) = curlj(1,i) - (dgradtmp(3,2)-dgradtmp(2,3))
        curlj(2,i) = curlj(2,i) - (dgradtmp(1,3)-dgradtmp(3,1))
        curlj(3,i) = curlj(3,i) - (dgradtmp(2,1)-dgradtmp(1,2))
        gradrho(1:3,i) = gradrho(1:3,i) - dgradtmp(4,1:3)
      enddo
!$OMP END PARALLEL DO      

!      print *, "finished pot eval"

      return
      end subroutine lpcomp_virtualcasing_addsub
!
!
!
!
