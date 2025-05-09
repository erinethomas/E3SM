#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module global_norms_mod

  use kinds, only : iulog
  use edgetype_mod, only : EdgeBuffer_t

  implicit none
  private
  save

  public :: l1_snorm
  public :: l2_snorm
  public :: linf_snorm

  public :: l1_vnorm
  public :: l2_vnorm
  public :: linf_vnorm

  public :: print_cfl
  public :: dss_hvtensor
  public :: test_global_integral
  public :: global_integral
  public :: wrap_repro_sum

  private :: global_maximum

  ! EdgeBuffer_t variables are shared by all thread,
  ! therefore we need this to be module static
  type (EdgeBuffer_t), private :: edgebuf

contains


  ! ================================
  ! global_integral:
  !
  ! eq 81 in Williamson, et. al. p 218
  ! for spectral elements
  !
  ! ================================
  ! --------------------------
  function global_integral(elem, h,hybrid,npts,nets,nete) result(I_sphere)
    use kinds,       only : real_kind
    use hybrid_mod,  only : hybrid_t
    use element_mod, only : element_t
    use dimensions_mod, only : np, nelemd
    use physical_constants, only : dd_pi, domain_size
    use parallel_mod, only: global_shared_buf, global_shared_sum

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind) :: I_sphere

    real (kind=real_kind) :: I_priv
    real (kind=real_kind) :: I_shared
    common /gblintcom/I_shared

    ! Local variables

    integer :: ie,j,i
    real(kind=real_kind) :: I_tmp(1)

    real (kind=real_kind) :: da
    real (kind=real_kind) :: J_tmp(nets:nete)
!
! This algorythm is independent of thread count and task count.
! This is a requirement of consistancy checking in cam.
!
    J_tmp = 0.0D0

!JMD    print *,'global_integral: before loop'
       do ie=nets,nete
          do j=1,np
             do i=1,np
                da = elem(ie)%mp(i,j)*elem(ie)%metdet(i,j)
                J_tmp(ie) = J_tmp(ie) + da*h(i,j,ie)
             end do
          end do
       end do       
    do ie=nets,nete
      global_shared_buf(ie,1) = J_tmp(ie)
    enddo
!JMD    print *,'global_integral: before wrap_repro_sum'
    call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
!JMD    print *,'global_integral: after wrap_repro_sum'
    I_tmp = global_shared_sum(1)
!JMD    print *,'global_integral: after global_shared_sum'

    I_sphere = I_tmp(1)/domain_size


  end function global_integral

  ! ================================
  ! test_global_integral:
  !
  ! test that the global integral of 
  ! the area of the sphere is 1.
  !
  ! ================================

  subroutine test_global_integral(elem,hybrid,nets,nete,mindxout)
    use kinds,       only : real_kind
    use hybrid_mod,  only : hybrid_t
    use element_mod, only : element_t
    use dimensions_mod, only : np,ne, nelem, nelemd
    use mesh_mod,     only : MeshUseMeshFile          

    use reduction_mod, only : ParallelMin,ParallelMax
    use physical_constants, only : scale_factor,dd_pi
    use parallel_mod, only : abortmp, global_shared_buf, global_shared_sum
    use bndry_mod, only : bndry_exchangeV
    use control_mod, only : geometry

    type(element_t)      , intent(inout) :: elem(:)
    integer              , intent(in) :: nets,nete
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind),intent(out), optional :: mindxout

    real (kind=real_kind)             :: I_sphere

    ! Local variables
    real (kind=real_kind) :: h(np,np,nets:nete)
    ! Element statisics
    real (kind=real_kind) :: min_area,max_area,avg_area, max_ratio
    real (kind=real_kind) :: min_min_dx, max_min_dx, avg_min_dx
    real (kind=real_kind) :: min_normDinv, max_normDinv
    real (kind=real_kind) :: min_len
    integer :: ie,corner, i, j,nlon


    h(:,:,nets:nete)=1.0D0

    ! Calculate surface area by integrating 1.0d0 over domain
    ! (Should be 1 for unit sphere and Lx * Ly for plane)
    I_sphere = global_integral(elem, h(:,:,nets:nete),hybrid,np,nets,nete)

    min_area=1d99
    max_area=0
    avg_area=0_real_kind

    max_ratio = 0

    min_normDinv=1d99
    max_normDinv=0

    min_min_dx=1d99
    max_min_dx=0
    avg_min_dx=0_real_kind

    do ie=nets,nete
       
       elem(ie)%area = sum(elem(ie)%spheremp(:,:))
       min_area=min(min_area,elem(ie)%area)
       max_area=max(max_area,elem(ie)%area)

       min_normDinv = min(min_normDinv,elem(ie)%normDinv)
       max_normDinv = max(max_normDinv,elem(ie)%normDinv)

       max_ratio   = max(max_ratio,elem(ie)%dx_long/elem(ie)%dx_short)


       min_min_dx = min(min_min_dx,elem(ie)%dx_short)
       max_min_dx = max(max_min_dx,elem(ie)%dx_short)


       global_shared_buf(ie,1) = elem(ie)%area
       global_shared_buf(ie,2) = elem(ie)%dx_short

    enddo

    min_area=ParallelMin(min_area,hybrid)
    max_area=ParallelMax(max_area,hybrid)

    min_normDinv=ParallelMin(min_normDinv,hybrid)
    max_normDinv=ParallelMax(max_normDinv,hybrid)

    max_ratio=ParallelMax(max_ratio,hybrid)

    min_min_dx=ParallelMin(min_min_dx,hybrid)
    max_min_dx=ParallelMax(max_min_dx,hybrid)

    call wrap_repro_sum(nvars=2, comm=hybrid%par%comm)

    avg_area = global_shared_sum(1)/dble(nelem)
    avg_min_dx = global_shared_sum(2)/dble(nelem)

    ! Physical units for area
    min_area = min_area*scale_factor*scale_factor/1000000_real_kind
    max_area = max_area*scale_factor*scale_factor/1000000_real_kind
    avg_area = avg_area*scale_factor*scale_factor/1000000_real_kind


    ! for an equation du/dt = i c u, leapfrog is stable for |c u dt| < 1
    ! Consider a gravity wave at the equator, c=340m/s  
    ! u = exp(i kmax x/ a ) with x = longitude,  and kmax =  pi a / dx, 
    ! u = exp(i pi x / dx ),   so du/dt = c du/dx becomes du/dt = i c pi/dx u
    ! stable for dt < dx/(c*pi)
    ! CAM 26 level AMIP simulation: max gravity wave speed 341.75 m/s
    if (hybrid%masterthread) then
       write(iulog,* )""
       write(iulog,* )"Running Global Integral Diagnostic..."
         write(iulog,*)"Area of manifold is",I_sphere
       write(iulog,*)"Should be 1.0 to round off..."
       write(iulog,'(a,f9.3)') 'Element area:  max/min',(max_area/min_area)
       if (.not.MeshUseMeshFile .and. geometry == "sphere") then
           write(iulog,'(a,f6.3,f8.2)') "Average equatorial node spacing (deg, km) = ", &
                dble(90)/dble(ne*(np-1)), DD_PI*scale_factor/(2000.0d0*dble(ne*(np-1)))
       end if
       write(iulog,'(a,2f9.3)') 'norm of Dinv (min, max): ', min_normDinv, max_normDinv
       write(iulog,'(a,1f8.2)') 'Max Dinv-based element distortion: ', max_ratio
       write(iulog,'(a,3f8.2)') 'dx based on Dinv svd:          ave,min,max = ', avg_min_dx, min_min_dx, max_min_dx
       write(iulog,'(a,3f8.2)') "dx based on sqrt element area: ave,min,max = ", &
                sqrt(avg_area)/(np-1),sqrt(min_area)/(np-1),sqrt(max_area)/(np-1)
    end if

    if(present(mindxout)) then
        ! min_len now based on norm(Dinv)
        min_len = 0.002d0*scale_factor/(dble(np-1)*max_normDinv)
        mindxout=1000_real_kind*min_len
    end if

  end subroutine test_global_integral


!------------------------------------------------------------------------------------

  ! ================================
  ! print_cfl:
  !
  ! Calculate / output CFL info
  ! (both advective and based on
  ! viscosity or hyperviscosity)
  !
  ! ================================

  subroutine print_cfl(elem,hybrid,nets,nete)
!
!   estimate various CFL limits
!
    use kinds,       only : real_kind
    use hybrid_mod,  only : hybrid_t
    use element_mod, only : element_t
#ifdef MODEL_THETA_L
    use element_state, only : nu_scale_top
#endif
    use dimensions_mod, only : np
    use quadrature_mod, only : gausslobatto, quadrature_t

    use reduction_mod, only : ParallelMin,ParallelMax
    use physical_constants, only : scale_factor_inv
    use control_mod, only : nu, nu_q, nu_div, hypervis_order, nu_top,  &
                            hypervis_scaling, dcmip16_mu,dcmip16_mu_s
    use control_mod, only : tstep_type

    type(element_t)      , intent(inout) :: elem(:)
    integer              , intent(in) :: nets,nete
    type (hybrid_t)      , intent(in) :: hybrid

    ! Element statisics
    real (kind=real_kind) :: min_max_dx,max_unif_dx   ! used for normalizing scalar HV
    real (kind=real_kind) :: max_normDinv  ! used for CFL
    real (kind=real_kind) :: normDinv_hypervis
    real (kind=real_kind) :: lambda_max, lambda_vis, min_gw, lambda, nu_div_actual, nu_top_actual
    integer :: ie
    type (quadrature_t)    :: gp


    ! Eigenvalues calculated by folks at UMich (Paul U & Jared W)
    select case (np)
        case (2)
            lambda_max = 0.5d0
            lambda_vis = 0.0d0  ! need to compute this
        case (3)
            lambda_max = 1.5d0
            lambda_vis = 12.0d0
        case (4)
            lambda_max = 2.74d0
            lambda_vis = 30.0d0
        case (5)
            lambda_max = 4.18d0
            lambda_vis = 91.6742d0
        case (6)
            lambda_max = 5.86d0
            lambda_vis = 190.1176d0
        case (7)
            lambda_max = 7.79d0
            lambda_vis = 374.7788d0
        case (8)
            lambda_max = 10.0d0
            lambda_vis = 652.3015d0
        case DEFAULT
            lambda_max = 0.0d0
            lambda_vis = 0.0d0
    end select

    if ((lambda_max.eq.0d0).and.(hybrid%masterthread)) then
        print*, "lambda_max not calculated for NP = ",np
        print*, "Estimate of gravity wave timestep will be incorrect"
    end if
    if ((lambda_vis.eq.0d0).and.(hybrid%masterthread)) then
        print*, "lambda_vis not calculated for NP = ",np
        print*, "Estimate of viscous CFLs will be incorrect"
    end if

    do ie=nets,nete
      elem(ie)%variable_hyperviscosity = 1.0
    end do

    gp=gausslobatto(np)
    min_gw = minval(gp%weights)

    deallocate(gp%weights)

    max_normDinv=0
    min_max_dx=1d99
    do ie=nets,nete
        max_normDinv = max(max_normDinv,elem(ie)%normDinv)
        min_max_dx = min(min_max_dx,elem(ie)%dx_long)
    enddo
    max_normDinv=ParallelMax(max_normDinv,hybrid)
    min_max_dx=ParallelMin(min_max_dx,hybrid)


    if (hypervis_scaling/=0) then
       ! tensorHV.  New eigenvalues are the eigenvalues of the tensor V
       ! formulas here must match what is in cube_mod.F90
       ! for tensorHV, we scale out the rearth dependency
       lambda = max_normDinv**2
       normDinv_hypervis = (lambda_vis**2) * (max_normDinv**4) * &
            (lambda**(-hypervis_scaling/2) )
    else
       ! constant coefficient formula:
       normDinv_hypervis = (lambda_vis**2) * (scale_factor_inv*max_normDinv)**4
    endif

     if (hybrid%masterthread) then
       write(iulog,'(a,f10.2)') 'CFL estimates in terms of S=time step stability region'
       write(iulog,'(a,f10.2)') '(i.e. advection w/leapfrog: S=1, viscosity w/forward Euler: S=2)'
       write(iulog,'(a,f10.2,a)') 'SSP preservation (120m/s) RKSSP euler step dt  < S *', &
            min_gw/(120.0d0*max_normDinv*scale_factor_inv),'s'
       write(iulog,'(a,f10.2,a)') 'Stability: advective (120m/s)   dt_tracer < S *',&
            1/(120.0d0*max_normDinv*lambda_max*scale_factor_inv),'s'
       write(iulog,'(a,f10.2,a)') 'Stability: advective (120m/s)   dt_tracer < S *', &
                                   1/(120.0d0*max_normDinv*lambda_max*scale_factor_inv),'s'
       write(iulog,'(a,f10.2,a)') 'Stability: gravity wave(342m/s)   dt_dyn  < S *', &
                                   1/(342.0d0*max_normDinv*lambda_max*scale_factor_inv),'s'
       if (nu>0) then
          if (hypervis_order==1) then
              write(iulog,'(a,f10.2,a)') 'Stability: viscosity dt < S *',&
                   1/(nu*((scale_factor_inv*max_normDinv)**2)*lambda_vis),'s'
          endif
          if (hypervis_order==2) then
             !  dt < S  1/nu*normDinv
             write(iulog,'(a,f10.2,a)') "Stability: nu_q   hyperviscosity dt < S *", 1/(nu_q*normDinv_hypervis),'s'
             write(iulog,'(a,f10.2,a)') "Stability: nu_vor hyperviscosity dt < S *", 1/(nu*normDinv_hypervis),'s'
#ifdef MODEL_THETA_L
             nu_div_actual = nu_div
#else             
             ! bug in preqx nu_div implimentation:
             ! we apply nu_ration=(nu_div/nu) in laplace, so it is applied 2x
             ! making the effective nu_div = nu * (nu_div/nu)**2 
             ! should be fixed - but need all CAM defaults adjusted, 
             nu_div_actual = nu_div**2/nu
#endif
             write(iulog,'(a,f10.2,a)') "Stability: nu_div hyperviscosity dt < S *",&
                  1/(nu_div_actual*normDinv_hypervis),'s'
          endif
       endif
       if(nu_top>0) then
#ifdef MODEL_THETA_L
          nu_top_actual=maxval(nu_scale_top)*nu_top
          write(iulog,'(a,f10.2,a)') 'scaled nu_top viscosity CFL: dt < S*', &
               1.0d0/(nu_top_actual*((scale_factor_inv*max_normDinv)**2)*lambda_vis),'s'
#else
          nu_top_actual=4*nu_top
          write(iulog,'(a,f10.2,a)') '4*nu_top viscosity CFL: dt < S*', &
               1.0d0/(nu_top_actual*((scale_factor_inv*max_normDinv)**2)*lambda_vis),'s'
#endif
       end if

      if(dcmip16_mu>0)  write(iulog,'(a,f10.2,a)') 'dcmip16_mu   viscosity CFL: dt < S*', &
           1.0d0/(dcmip16_mu*  ((scale_factor_inv*max_normDinv)**2)*lambda_vis),'s'
      if(dcmip16_mu_s>0)write(iulog,'(a,f10.2,a)') 'dcmip16_mu_s viscosity CFL: dt < S*', &
           1.0d0/(dcmip16_mu_s*((scale_factor_inv*max_normDinv)**2)*lambda_vis),'s'

      write(iulog,*) 'tstep_type = ',tstep_type
    end if

  end subroutine print_cfl

  ! ================================
  ! dss_hvtensor:
  !
  ! apply dss and bilinear projection
  ! to tensor coefficients
  !
  ! ================================

  subroutine dss_hvtensor(elem,hybrid,nets,nete)
    use kinds,       only : real_kind
    use hybrid_mod,  only : hybrid_t
    use element_mod, only : element_t

    use dimensions_mod, only : np
    use quadrature_mod, only : gausslobatto, quadrature_t

    use control_mod, only : hypervis_scaling
    use edge_mod, only : initedgebuffer, FreeEdgeBuffer, edgeVpack, edgeVunpack
    use bndry_mod, only : bndry_exchangeV

    type(element_t)      , intent(inout) :: elem(:)
    integer              , intent(in) :: nets,nete
    type (hybrid_t)      , intent(in) :: hybrid

    ! Element statisics
    real (kind=real_kind) :: x, y, noreast, nw, se, sw
    real (kind=real_kind), dimension(np,np,nets:nete) :: zeta
    integer :: ie, i, j, rowind, colind
    type (quadrature_t)    :: gp

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  TENSOR, RESOLUTION-AWARE HYPERVISCOSITY
    !  The tensorVisc() array is computed in cube_mod.F90
    !  this block of code will DSS it so the tensor if C0
    !  and also make it bilinear in each element.
    !  Oksana Guba
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (hypervis_scaling /= 0) then
      call initEdgeBuffer(hybrid%par,edgebuf,elem,1)
      do rowind=1,2
        do colind=1,2
          do ie=nets,nete
            zeta(:,:,ie) = elem(ie)%tensorVisc(:,:,rowind,colind)*elem(ie)%spheremp(:,:)
            call edgeVpack(edgebuf,zeta(1,1,ie),1,0,ie)
          end do

          call bndry_exchangeV(hybrid,edgebuf)
          do ie=nets,nete
            call edgeVunpack(edgebuf,zeta(1,1,ie),1,0,ie)
            elem(ie)%tensorVisc(:,:,rowind,colind) = zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
          end do
        enddo !rowind
      enddo !colind
      call FreeEdgeBuffer(edgebuf)

      gp=gausslobatto(np)

      !IF BILINEAR MAP OF V NEEDED
      do rowind=1,2
        do colind=1,2
          ! replace hypervis w/ bilinear based on continuous corner values
          do ie=nets,nete
            noreast = elem(ie)%tensorVisc(np,np,rowind,colind)
            nw = elem(ie)%tensorVisc(1,np,rowind,colind)
            se = elem(ie)%tensorVisc(np,1,rowind,colind)
            sw = elem(ie)%tensorVisc(1,1,rowind,colind)
            do i=1,np
              x = gp%points(i)
              do j=1,np
                y = gp%points(j)
                elem(ie)%tensorVisc(i,j,rowind,colind) = 0.25d0*( &
                (1.0d0-x)*(1.0d0-y)*sw + &
                (1.0d0-x)*(y+1.0d0)*nw + &
                (x+1.0d0)*(1.0d0-y)*se + &
                (x+1.0d0)*(y+1.0d0)*noreast)
              end do
            end do
          end do
        enddo !rowind
      enddo !colind

      deallocate(gp%points)
      deallocate(gp%weights)
    endif

   end subroutine dss_hvtensor

  ! ================================
  ! global_maximum:
  !
  ! Find global maximum on sphere
  !
  ! ================================

  function global_maximum(h,hybrid,npts,nets,nete) result(Max_sphere)

    use kinds, only : real_kind
    use hybrid_mod, only : hybrid_t
    use reduction_mod, only : red_max, pmax_mt

    integer              , intent(in) :: npts,nets,nete     
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind) :: Max_sphere

    ! Local variables

    real (kind=real_kind) :: redp(1)

    Max_sphere = MAXVAL(h(:,:,nets:nete))

    redp(1) = Max_sphere
    call pmax_mt(red_max,redp,1,hybrid)
    Max_sphere = red_max%buf(1)

  end function global_maximum

  ! ==========================================================
  ! l1_snorm:
  !
  ! computes the l1 norm per Williamson et al, p. 218 eq(8)
  ! for a scalar quantity
  ! ===========================================================

  function l1_snorm(elem, h,ht,hybrid,npts,nets,nete) result(l1)

    use kinds, only : real_kind
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: ht(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: l1     

    ! Local variables

    real (kind=real_kind) :: dhabs(npts,npts,nets:nete)
    real (kind=real_kind) :: htabs(npts,npts,nets:nete)
    real (kind=real_kind) :: dhabs_int
    real (kind=real_kind) :: htabs_int
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dhabs(i,j,ie) = ABS(h(i,j,ie)-ht(i,j,ie))
             htabs(i,j,ie) = ABS(ht(i,j,ie))
          end do
       end do
    end do

    dhabs_int = global_integral(elem, dhabs(:,:,nets:nete),hybrid,npts,nets,nete)
    htabs_int = global_integral(elem, htabs(:,:,nets:nete),hybrid,npts,nets,nete)

    l1 = dhabs_int/htabs_int

  end function l1_snorm

  ! ===========================================================
  ! l1_vnorm:
  !
  ! computes the l1 norm per Williamson et al, p. 218 eq(97),
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function l1_vnorm(elem, v,vt,hybrid,npts,nets,nete) result(l1)
    use kinds, only : real_kind
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in), target :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: l1     

    ! Local variables

    real (kind=real_kind), dimension(:,:,:,:), pointer :: met
    real (kind=real_kind) :: dvsq(npts,npts,nets:nete)
    real (kind=real_kind) :: vtsq(npts,npts,nets:nete)
    real (kind=real_kind) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: dv1,dv2
    real (kind=real_kind) :: vt1,vt2
    real (kind=real_kind) :: dvsq_int
    real (kind=real_kind) :: vtsq_int

    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met
       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(i,j,1,1)*dv1 + met(i,j,1,2)*dv2
             dvco(i,j,2) = met(i,j,2,1)*dv1 + met(i,j,2,2)*dv2

             vtco(i,j,1) = met(i,j,1,1)*vt1 + met(i,j,1,2)*vt2
             vtco(i,j,2) = met(i,j,2,1)*vt1 + met(i,j,2,2)*vt2

             dvsq(i,j,ie) = SQRT(dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2)
             vtsq(i,j,ie) = SQRT(vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2)

          end do
       end do
    end do

    dvsq_int = global_integral(elem, dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_int = global_integral(elem, vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    l1 = dvsq_int/vtsq_int

  end function l1_vnorm

  ! ==========================================================
  ! l2_snorm:
  !
  ! computes the l2 norm per Williamson et al, p. 218 eq(83)
  ! for a scalar quantity on the pressure grid.
  !
  ! ===========================================================

  function l2_snorm(elem, h,ht,hybrid,npts,nets,nete) result(l2)
    use kinds, only : real_kind
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t), intent(in) :: elem(:)	
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: ht(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: l2   

    ! Local variables

    real (kind=real_kind) :: dh2(npts,npts,nets:nete)
    real (kind=real_kind) :: ht2(npts,npts,nets:nete)
    real (kind=real_kind) :: dh2_int
    real (kind=real_kind) :: ht2_int
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dh2(i,j,ie)=(h(i,j,ie)-ht(i,j,ie))**2
             ht2(i,j,ie)=ht(i,j,ie)**2
          end do
       end do
    end do

    dh2_int = global_integral(elem,dh2(:,:,nets:nete),hybrid,npts,nets,nete)
    ht2_int = global_integral(elem,ht2(:,:,nets:nete),hybrid,npts,nets,nete)

    l2 = SQRT(dh2_int)/SQRT(ht2_int)

  end function l2_snorm

  ! ==========================================================
  ! l2_vnorm:
  !
  ! computes the l2 norm per Williamson et al, p. 219 eq(98)
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function l2_vnorm(elem, v,vt,hybrid,npts,nets,nete) result(l2)
    use kinds, only : real_kind
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in), target :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: l2

    ! Local variables

    real (kind=real_kind), dimension(:,:,:,:), pointer :: met
    real (kind=real_kind) :: dvsq(npts,npts,nets:nete)
    real (kind=real_kind) :: vtsq(npts,npts,nets:nete)
    real (kind=real_kind) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: dv1,dv2
    real (kind=real_kind) :: vt1,vt2
    real (kind=real_kind) :: dvsq_int
    real (kind=real_kind) :: vtsq_int
    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met
       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(i,j,1,1)*dv1 + met(i,j,1,2)*dv2
             dvco(i,j,2) = met(i,j,2,1)*dv1 + met(i,j,2,2)*dv2

             vtco(i,j,1) = met(i,j,1,1)*vt1 + met(i,j,1,2)*vt2
             vtco(i,j,2) = met(i,j,2,1)*vt1 + met(i,j,2,2)*vt2

             dvsq(i,j,ie) = dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2
             vtsq(i,j,ie) = vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2

          end do
       end do
    end do

    dvsq_int = global_integral(elem, dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_int = global_integral(elem, vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    l2 = SQRT(dvsq_int)/SQRT(vtsq_int)

  end function l2_vnorm

  ! ==========================================================
  ! linf_snorm:
  !
  ! computes the l infinity norm per Williamson et al, p. 218 eq(84)
  ! for a scalar quantity on the pressure grid...
  !
  ! ===========================================================

  function linf_snorm(h,ht,hybrid,npts,nets,nete) result(linf)
    use kinds, only : real_kind
    use hybrid_mod, only : hybrid_t
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: ht(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: linf    

    ! Local variables

    real (kind=real_kind) :: dhabs(npts,npts,nets:nete)
    real (kind=real_kind) :: htabs(npts,npts,nets:nete)
    real (kind=real_kind) :: dhabs_max
    real (kind=real_kind) :: htabs_max
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dhabs(i,j,ie)=ABS(h(i,j,ie)-ht(i,j,ie))
             htabs(i,j,ie)=ABS(ht(i,j,ie))
          end do
       end do
    end do

    dhabs_max = global_maximum(dhabs(:,:,nets:nete),hybrid,npts,nets,nete)
    htabs_max = global_maximum(htabs(:,:,nets:nete),hybrid,npts,nets,nete)

    linf = dhabs_max/htabs_max

  end function linf_snorm


  ! ==========================================================
  ! linf_vnorm:
  !
  ! computes the linf norm per Williamson et al, p. 218 eq(99),
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function linf_vnorm(elem,v,vt,hybrid,npts,nets,nete) result(linf)
    use kinds, only : real_kind
    use hybrid_mod, only : hybrid_t
    use element_mod, only : element_t

    type(element_t)      , intent(in), target :: elem(:) 
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: linf     

    ! Local variables

    real (kind=real_kind), dimension(:,:,:,:), pointer :: met
    real (kind=real_kind) :: dvsq(npts,npts,nets:nete)
    real (kind=real_kind) :: vtsq(npts,npts,nets:nete)
    real (kind=real_kind) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: dv1,dv2
    real (kind=real_kind) :: vt1,vt2
    real (kind=real_kind) :: dvsq_max
    real (kind=real_kind) :: vtsq_max
    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met

       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(i,j,1,1)*dv1 + met(i,j,1,2)*dv2
             dvco(i,j,2) = met(i,j,2,1)*dv1 + met(i,j,2,2)*dv2

             vtco(i,j,1) = met(i,j,1,1)*vt1 + met(i,j,1,2)*vt2
             vtco(i,j,2) = met(i,j,2,1)*vt1 + met(i,j,2,2)*vt2

             dvsq(i,j,ie) = SQRT(dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2)
             vtsq(i,j,ie) = SQRT(vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2)

          end do
       end do
    end do

    dvsq_max = global_maximum(dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_max = global_maximum(vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    linf = dvsq_max/vtsq_max

  end function linf_vnorm


  subroutine wrap_repro_sum (nvars, comm, nsize)
    use dimensions_mod, only: nelemd
#ifdef CAM
    use shr_reprosum_mod, only: repro_sum => shr_reprosum_calc
#else
    use repro_sum_mod, only: repro_sum
#endif
    use parallel_mod, only: global_shared_buf, global_shared_sum, nrepro_vars, abortmp

    implicit none

    integer :: nvars            !  number of variables to be summed (cannot exceed nrepro_vars)
    integer :: comm             !  mpi communicator
    integer, optional :: nsize  !  local buffer size (defaults to nelemd - number of elements in mpi task)

    integer nsize_use,n,i

    if (present(nsize)) then
       nsize_use = nsize
    else
       nsize_use = nelemd
    endif
    if (nvars .gt. nrepro_vars) call abortmp('repro_sum_buffer_size exceeded')

#if (defined HORIZ_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif

#ifndef CAM
    ! CAM already does this, no need to do it twice
    do n=1,nvars
       do i=1,nsize_use
          if (global_shared_buf(i,n) /= global_shared_buf(i,n) ) then
               print *, "var,nvars:",n,nvars
               call abortmp('NaNs detected in repro sum input')
          endif
       enddo
    enddo
#endif    

! Repro_sum contains its own OpenMP, so only one thread should call it (AAM)
    call repro_sum(global_shared_buf, global_shared_sum, nsize_use, nelemd, nvars, commid=comm)

#if (defined HORIZ_OPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif

    end subroutine wrap_repro_sum

end module global_norms_mod
