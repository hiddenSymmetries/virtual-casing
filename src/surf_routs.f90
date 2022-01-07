!
!  This file contains the following user callable routines:
!    
!
!
!
subroutine get_second_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sfform)
!
!------------------------
!  This subroutine computes the second fundamental form at
!  the discretization nodes on the surface.
!
!  The second fundamental form is
!  
!  .. math::
!    
!    \begin{bmatrix} x_{uu} \cdot n & x_{uv} \cdot n \\
!    x_{uv} \cdot n & x_{vv} \cdot n \end{bmatrix}
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - sfform: double precision(2,2,npts)
!        second fundamental form at the discretization nodes
!--------------------------
!
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: sfform(2,2,npts)
  integer i
  real *8 L,M,N
  real *8, allocatable :: dxuv(:,:)
  real *8, allocatable :: dxuv2(:,:,:)


  allocate(dxuv(6,npts))
! Calculating x_{uu}, x_{uv}, x_{uv} stored in dxuv2 
 
  do i=1,npts
    dxuv(1,i) = srcvals(4,i)  
    dxuv(2,i) = srcvals(5,i)  
    dxuv(3,i) = srcvals(6,i)  
    dxuv(4,i) = srcvals(7,i)  
    dxuv(5,i) = srcvals(8,i)  
    dxuv(6,i) = srcvals(9,i)   
  enddo

  allocate(dxuv2(6,2,npts))

  call get_surf_uv_grad(6,npatches,norders,ixyzs,iptype,npts,dxuv,dxuv2)



  do i=1,npts
    L = dxuv2(1,1,i)*srcvals(10,i) + dxuv2(2,1,i)*srcvals(11,i) + dxuv2(3,1,i)*srcvals(12,i) ! Calculation of L, M, N. L = x_uu \cdot n 
    M = dxuv2(1,2,i)*srcvals(10,i) + dxuv2(2,2,i)*srcvals(11,i) + dxuv2(3,2,i)*srcvals(12,i)  ! M = x_uv \cdot n
    N = dxuv2(4,2,i)*srcvals(10,i) + dxuv2(5,2,i)*srcvals(11,i) + dxuv2(6,2,i)*srcvals(12,i)  ! N = x_vv \cdot n
    sfform(1,1,i) = L
    sfform(2,1,i) = M
    sfform(1,2,i) = M
    sfform(2,2,i) = N
  enddo

  return
end subroutine get_second_fundamental_form
!
!
!
subroutine get_mean_curvature(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,mean_curv)
!
!------------------------
!  This subroutine computes the mean curvature at
!  the discretization nodes on the surface.
!
!  The mean curvature is
!  
!  .. math::
!    
!    0.5*Trace(II \cdot I^{-1}) \\
!    
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - mean_curv: double precision(npts)
!        mean curvature at the discretization nodes
!--------------------------
!
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: mean_curv(npts)
  integer i
  real *8, allocatable :: ffform(:,:,:),ffforminv(:,:,:)
  real *8, allocatable :: sfform(:,:,:)

  allocate(ffform(2,2,npts))

  call get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)



  allocate(ffforminv(2,2,npts))

  call get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)

  allocate(sfform(2,2,npts))

  call get_second_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sfform)
  
!  print *,"Point on surface:", srcvals(1,3),srcvals(2,3), srcvals(3,3) 
!  print *,"First fundamental form=", ffform(:,:, 3) 
!  print *,"Inverse first fundamental form=", ffforminv(:,:, 3)
!  print *,"Second fundamental form=", sfform(:,:, 3)
 


! Calculating mean curvature 
 
  do i=1,npts
    mean_curv(i) = -0.5*(sfform(1,1,i)*ffforminv(1,1,i) + &
                     sfform(1,2,i)*ffforminv(2,1,i) + &
                     sfform(2,1,i)*ffforminv(1,2,i) + &
                     sfform(2,2,i)*ffforminv(2,2,i))
  enddo
!  print *,"Mean=", mean_curv(3)
 
  return
end subroutine get_mean_curvature

!
!
!
!
!
!
subroutine surf_grad(npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,fin,gradf)
!
!-----------------------------
!  Compute the surface gradient of scalar function fin  
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!    - fin: double precision (npts)
!         vector function on surface
!  Output arguments:
!
!    - gradf: double precision(3,npts)
!        surface gradient 
!        
!-----------------------------
!
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),fin(npts)
  real *8, intent(out) :: gradf(3,npts)
  real *8, allocatable :: ffforminv(:,:,:)
  real *8, allocatable :: dfuv(:,:,:)
  real *8 E,F1,F2,G,W_sq,a,b,fdu,fdv
  integer i,istart,npols,j,l

  allocate(ffforminv(2,2,npts))



  call get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)

  allocate(dfuv(1,2,npts))
! Calculating f_{u}, f_{v}} stored in dfuv 
 
  call get_surf_uv_grad(1,npatches,norders,ixyzs,iptype,npts,fin,dfuv)



  do i=1,npts
    E = ffforminv(1,1,i)
    F1 = ffforminv(1,2,i)
    F2 = ffforminv(2,1,i)
    G = ffforminv(2,2,i)
    fdu = dfuv(1,1,i)
    fdv = dfuv(1,2,i)
    a = E*fdu+F1*fdv
    b = F2*fdu+G*fdv

    gradf(1,i) = a*srcvals(4,i) + b*srcvals(7,i)
    gradf(2,i) = a*srcvals(5,i) + b*srcvals(8,i)
    gradf(3,i) = a*srcvals(6,i) + b*srcvals(9,i)


                 
  enddo 

  return
end subroutine surf_grad
!




subroutine surf_grad2(npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,fin,gradf)
!
!-----------------------------
!  Compute the surface gradient of scalar function fin  
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!    - fin: double precision (npts)
!         vector function on surface
!  Output arguments:
!
!    - gradf: double precision(3,npts)
!        surface gradient 
!        
!-----------------------------
!
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),fin(3,npts)
  real *8, intent(out) :: gradf(3,npts)
  real *8, allocatable :: ffform(:,:,:)
  real *8, allocatable :: dfuv(:,:,:)
  real *8 E,F,G,W_sq,a(3),b(3),fdu,fdv
  integer i,istart,npols,j,l

  allocate(ffform(2,2,npts))



  call get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)

  allocate(dfuv(1,2,npts))
! Calculating f_{u}, f_{v}} stored in dfuv 
 
  call get_surf_uv_grad(1,npatches,norders,ixyzs,iptype,npts,fin,dfuv)



  do i=1,npts
    E = ffform(1,1,i)
    F = ffform(1,2,i)
    G = ffform(2,2,i)
    W_sq = E*G - F**2
    fdu = dfuv(1,1,i)
    fdv = dfuv(1,2,i)
    a(1) = (G*srcvals(4,i)-F*srcvals(7,i))/W_sq 
    a(2) = (G*srcvals(5,i)-F*srcvals(8,i))/W_sq 
    a(3) = (G*srcvals(6,i)-F*srcvals(9,i))/W_sq 
    b(1) = (E*srcvals(7,i)-F*srcvals(4,i))/W_sq 
    b(2) = (E*srcvals(8,i)-F*srcvals(5,i))/W_sq 
    b(3) = (E*srcvals(9,i)-F*srcvals(6,i))/W_sq 


    gradf(1,i) = a(1)*fdu+b(1)*fdv
    gradf(2,i) = a(2)*fdu+b(2)*fdv
    gradf(3,i) = a(3)*fdu+b(3)*fdv



                 
  enddo 

  return
end subroutine surf_grad2
!

!
!
!
!
!
subroutine surf_div(npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,fin,divf)
!
!-----------------------------
!  Compute the surface divergence of vector function fin  
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!    - fin: double precision (3,npts)
!         vector function on surface
!  Output arguments:
!
!    - divf: double precision(npts)
!        surface divergence 
!        
!-----------------------------
!
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),fin(3,npts)
  real *8, intent(out) :: divf(npts)
  real *8, allocatable :: ffform(:,:,:)
  real *8, allocatable :: dfuv(:,:,:)
  real *8 E,F,G,W_sq,a(3),b(3),fdu(3),fdv(3)
  integer i,istart,npols,j,l

  allocate(ffform(2,2,npts))



  call get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)

  allocate(dfuv(3,2,npts))
! Calculating f_{u}, f_{v}} stored in dfuv 
 
  call get_surf_uv_grad(3,npatches,norders,ixyzs,iptype,npts,fin,dfuv)



  do i=1,npts
    E = ffform(1,1,i)
    F = ffform(1,2,i)
    G = ffform(2,2,i)
    W_sq = E*G - F**2
    fdu(1) = dfuv(1,1,i)
    fdu(2) = dfuv(2,1,i)
    fdu(3) = dfuv(3,1,i)
    fdv(1) = dfuv(1,2,i)
    fdv(2) = dfuv(2,2,i)
    fdv(3) = dfuv(3,2,i)
    a(1) = (G*fdu(1)-F*fdv(1))/W_sq 
    a(2) = (G*fdu(2)-F*fdv(2))/W_sq 
    a(3) = (G*fdu(3)-F*fdv(3))/W_sq 
    b(1) = (E*fdv(1)-F*fdu(1))/W_sq 
    b(2) = (E*fdv(2)-F*fdu(2))/W_sq 
    b(3) = (E*fdv(3)-F*fdu(3))/W_sq 


    divf(i) = a(1)*srcvals(4,i)+a(2)*srcvals(5,i)+a(3)*srcvals(6,i)+ &
              b(1)*srcvals(7,i)+b(2)*srcvals(8,i)+b(3)*srcvals(9,i)
                 
  enddo 

  return
end subroutine surf_div
!
!
!
!
!
!

subroutine get_ab_cycles_torusparam(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ipars,m,na,avals,awts,apatches,auv,nb,bvals, &
  bwts,bpatches,buv)

  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
  integer npts,m,ipars(2)
  real *8 srcvals(12,npts),srccoefs(9,npts)
  integer apatches(na),bpatches(nb)
  real *8 avals(9,na),awts(na),auv(2,na)
  real *8 bvals(9,nb),bwts(nb),buv(2,nb)
  integer, intent(in) :: na,nb
  integer nu,nv,i,j,k,l
  real *8, allocatable :: xnodes(:),wts(:)
  real *8, allocatable :: polsu(:,:),polsv(:,:)
  real *8 uv(2)
  real *8 umat,vmat
  real *8 wtmp(12),rr
  integer ipt,itri,lpt,nmax,norder,npmax,npols
  integer itype

  itype = 1
  allocate(xnodes(m),wts(m))
  call legeexps(itype,m,xnodes,umat,vmat,wts)

  nu = ipars(1)
  nv = ipars(2)

  nmax = 20
  npmax = (nmax+1)*(nmax+2)/2
  allocate(polsu(npmax,m))
  allocate(polsv(npmax,m))
  do i=1,m
    uv(1) = (xnodes(i)+1)/2
    uv(2) = 0
    call koorn_pols(uv,nmax,npmax,polsu(1,i))

    uv(1) = 0
    uv(2) = (xnodes(i)+1)/2
    call koorn_pols(uv,nmax,npmax,polsv(1,i))
  enddo


  do i=1,nv
    itri = 2*i-1
    norder = norders(itri)
    npols = (norder+1)*(norder+2)/2
    do j=1,m
      ipt = (i-1)*m + j
      apatches(ipt) = itri
      auv(1,ipt) = 0
      auv(2,ipt) = (xnodes(j)+1)/2
      awts(ipt) = wts(j)/2
      wtmp(1:12) = 0
      do l=1,npols
        lpt = ixyzs(itri)+l-1
        wtmp(1:9) = wtmp(1:9) + srccoefs(1:9,lpt)*polsv(l,j)
      enddo
      call cross_prod3d(wtmp(4),wtmp(7),wtmp(10))
      rr = sqrt(wtmp(10)**2 + wtmp(11)**2 + wtmp(12)**2)
      avals(1:3,ipt) = wtmp(1:3)
      avals(4:6,ipt) = wtmp(7:9)
      avals(7:9,ipt) = wtmp(10:12)/rr
    enddo
  enddo

  do i=1,nu
    itri = 2*nv*(i-1) + 1
    norder = norders(itri)
    npols = (norder+1)*(norder+2)/2
    do j=1,m
      ipt = (i-1)*m + j
      bpatches(ipt) = itri
      buv(1,ipt) = (xnodes(j)+1)/2
      buv(2,ipt) = 0
      bwts(ipt) = wts(j)/2
      wtmp(1:12) = 0
      do l=1,npols
        lpt = ixyzs(itri)+l-1
        wtmp(1:9) = wtmp(1:9) + srccoefs(1:9,lpt)*polsu(l,j)
      enddo
      call cross_prod3d(wtmp(4),wtmp(7),wtmp(10))
      rr = sqrt(wtmp(10)**2 + wtmp(11)**2 + wtmp(12)**2)
      bvals(1:6,ipt) = wtmp(1:6)
      bvals(7:9,ipt) = wtmp(10:12)/rr
    enddo
  enddo


end subroutine get_ab_cycles_torusparam
!
!
!
!
!
!
!
!

subroutine get_cycle_readfile(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ngenus,m,nc,fname,cvals,cwts,cpatches,cuv, &
  icxyzs)

  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
  integer npts,m,ngenus
  real *8 srcvals(12,npts),srccoefs(9,npts)
  integer cpatches(nc)
  integer icxyzs(ngenus+1)
  real *8 cvals(9,nc),cwts(nc),cuv(2,nc)
  integer, intent(in) :: nc
  integer nu,nv,i,j,k,l,iuv(2,-3:3)
  real *8, allocatable :: xnodes(:),wts(:)
  real *8, allocatable :: pols(:,:,:),uvnodes(:,:,:),wtsall(:,:)
  real *8 uv(2)
  real *8 umat,vmat
  real *8 wtmp(12),rr
  integer ipt,itri,lpt,nmax,norder,npmax,npols
  integer itype,ngenus0,ne,iedge,netot
  character (len=*) fname

  itype = 1
  allocate(xnodes(m),wts(m))
  call legeexps(itype,m,xnodes,umat,vmat,wts)


  nmax = 20
  npmax = (nmax+1)*(nmax+2)/2
  allocate(pols(npmax,m,-3:3))
  allocate(uvnodes(2,m,-3:3))

  iuv(1,1) = 1
  iuv(2,1) = 0
  
  iuv(1,-1) = -1
  iuv(2,-1) = 0

  iuv(1,3) = 0
  iuv(2,3) = 1
  
  iuv(1,-3) = 0
  iuv(2,-3) = -1

  iuv(1,2) = -1
  iuv(2,2) = 1

  iuv(1,-2) = 1
  iuv(2,-2) = -1


  do i=1,m
    uv(1) = (xnodes(i)+1)/2
    uv(2) = 0
    uvnodes(1,i,1) = uv(1) 
    uvnodes(2,i,1) = uv(2) 
    call koorn_pols(uv,nmax,npmax,pols(1,i,1))
    
    uv(1) = 1.0d0-(xnodes(i)+1)/2
    uv(2) = 0.0d0
    uvnodes(1,i,-1) = uv(1) 
    uvnodes(2,i,-1) = uv(2) 
    call koorn_pols(uv,nmax,npmax,pols(1,i,-1))

    uv(1) = 0
    uv(2) = (xnodes(i)+1)/2
    uvnodes(1,i,3) = uv(1)
    uvnodes(2,i,3) = uv(2)
    call koorn_pols(uv,nmax,npmax,pols(1,i,3))
    
    uv(1) = 0
    uv(2) = 1.0d0-(xnodes(i)+1)/2
    uvnodes(1,i,-3) = uv(1)
    uvnodes(2,i,-3) = uv(2)
    call koorn_pols(uv,nmax,npmax,pols(1,i,-3))
    
    uv(2) = (xnodes(i)+1)/2
    uv(1) = 1.0d0-uv(2) 
    uvnodes(1,i,2) = uv(1)
    uvnodes(2,i,2) = uv(2)
    call koorn_pols(uv,nmax,npmax,pols(1,i,2))

    uv(1) = (xnodes(i)+1)/2
    uv(2) = 1.0d0-uv(1) 
    uvnodes(1,i,-2) = uv(1)
    uvnodes(2,i,-2) = uv(2)
    call koorn_pols(uv,nmax,npmax,pols(1,i,-2))
  enddo

  open(unit=33,file=fname)
  read(33,*) ngenus0
  icxyzs(1) = 1
  
  netot = 0
  do i=1,ngenus
    read(33,*) ne
    icxyzs(i+1) = icxyzs(i) + ne*m
    netot = netot+ne
  enddo
  call prinf('ngenus=*',ngenus,1)

  do i=1,netot
    read(33,*) itri,iedge
    norder = norders(itri)
    npols = (norder+1)*(norder+2)/2
    do j=1,m
      ipt = (i-1)*m + j
      cpatches(ipt) = itri
      cuv(1,ipt) = uvnodes(1,j,iedge)
      cuv(2,ipt) = uvnodes(2,j,iedge) 
      cwts(ipt) = wts(j)/2
      wtmp(1:12) = 0
      do l=1,npols
        lpt = ixyzs(itri)+l-1
        wtmp(1:9) = wtmp(1:9) + srccoefs(1:9,lpt)*pols(l,j,iedge)
      enddo
      call cross_prod3d(wtmp(4),wtmp(7),wtmp(10))
      rr = sqrt(wtmp(10)**2 + wtmp(11)**2 + wtmp(12)**2)
      cvals(1:3,ipt) = wtmp(1:3)
      cvals(4:6,ipt) = wtmp(4:6)*iuv(1,iedge) + wtmp(7:9)*iuv(2,iedge)
      cvals(7:9,ipt) = wtmp(10:12)/rr
    enddo
  enddo
  close(33)



end subroutine get_cycle_readfile
!
!
!
!
!
!

subroutine vtk_curv_plot(n,nda,avals,fname,title)
!
! This subroutine writes a vtk file to plot a given a curve 
! as a collection of line segments
!
!
!  Input arguments:
!    - n: integer
!        number of points
!    - nda: integer
!        leading dimension of data array
!    - avals: real *8 (nda,n)
!        curve to be plotted, the first three components
!        must be xyz coordinates
!    - fname: character (len=*)
!        file name where vtk output should be written
!
  implicit none
  integer, intent(in) :: nda,n
  real *8, intent(in) :: avals(nda,n)
  character (len=*), intent(in) :: fname,title


  integer i,j,k,l,n0,npuv,ipatch,ipt,i1,m,norder,npols,iunit1

  
  iunit1 = 877
  open(unit = iunit1, file=trim(fname), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i9,a)') "POINTS ", n, " float"

  do i = 1,n
    write(iunit1,"(E11.5,2(2x,e11.5))") avals(1,i), avals(2,i), avals(3,i)
  end do

  write(iunit1,'(a,i9,i9)') "CELLS ", n, n*3

  do i=1,n
    if(i.ne.n)  write(iunit1,'(a,i9,i9)') "2 ", i-1, i
    if(i.eq.n)  write(iunit1,'(a,i9,i9)') "2 ", i-1, 0
  enddo

  write(iunit1,'(a,i9)') "CELL_TYPES ", n
  do ipatch = 1,n
    write(iunit1,'(a)') "3"
  end do

  close(iunit1)

end subroutine vtk_curv_plot
!
!
!
!
!


subroutine fun_surf_interp(nd,npatches,norders,ixyzs,iptype,npts, &
  f,na,apatches,auv,finterp)

  implicit real *8 (a-h,o-z)
  integer nd,npatches,npts
  integer norders(npatches),ixyzs(npatches+1),iptype(npatches)
  real *8 f(nd,npts)
  integer na,apatches(na)
  real *8 auv(2,na),finterp(nd,na)

  real *8, allocatable :: fcoefs(:,:)
  real *8, allocatable :: pols(:)

  allocate(pols(1000))
  

  allocate(fcoefs(nd,npts))
  call surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts,&
    f, fcoefs)


  do i=1,na
    ip = apatches(i)
    norder = norders(ip)
    npols = (norder+1)*(norder+2)/2

    call koorn_pols(auv(1,i),norder,npols,pols)
    do idim=1,nd
      finterp(idim,i) = 0
    enddo
    do l=1,npols
      lpt = ixyzs(ip)+l-1
      do idim=1,nd
        finterp(idim,i) = finterp(idim,i) + fcoefs(idim,lpt)*pols(l)
      enddo
    enddo
  enddo

end subroutine fun_surf_interp
!
!
!
!
!
subroutine geom_coefs_interp(npatches,norders,ixyzs,iptype,npts, &
  srccoefs,na,apatches, auv,srcinterp)

  implicit real *8 (a-h,o-z)
  integer nd,npatches,npts
  integer norders(npatches),ixyzs(npatches+1),iptype(npatches)
  real *8 srccoefs(9,npts)
  integer na,apatches(na)
  real *8 auv(2,na),srcinterp(12,na)
  real *8, allocatable :: pols(:)


  allocate(pols(1000))


  do i=1,na
    ip = apatches(i)
    norder = norders(ip)
    npols = (norder+1)*(norder+2)/2

    call koorn_pols(auv(1,i),norder,npols,pols)
    do idim=1,9
      srcinterp(idim,i) = 0
    enddo
    do l=1,npols
      lpt = ixyzs(ip)+l-1
      do idim=1,9
        srcinterp(idim,i) = srcinterp(idim,i) + srccoefs(idim,lpt)*pols(l)
      enddo
    enddo

    call cross_prod3d(srcinterp(4,i),srcinterp(7,i),srcinterp(10,i))
    rr = sqrt(srcinterp(10,i)**2 + srcinterp(11,i)**2 +  &
      srcinterp(12,i)**2)
    srcinterp(10,i) = srcinterp(10,i)/rr
    srcinterp(11,i) = srcinterp(11,i)/rr
    srcinterp(12,i) = srcinterp(12,i)/rr
  enddo

end subroutine geom_coefs_interp
