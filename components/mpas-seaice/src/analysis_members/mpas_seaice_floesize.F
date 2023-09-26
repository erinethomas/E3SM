!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_floesize
!
!> \brief MPAS-Seaice analysis mode member: floesize
!> \author Erin E Thomas
!> \date   15 September 2023
!> \details
!>  MPAS-Seaice analysis mode member: floesize
!>  Calculates diagnositcs for Floe Size Distribution
!>
!-----------------------------------------------------------------------

module seaice_floesize

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: seaice_bootstrap_floesize, &
             seaice_init_floesize, &
             seaice_precompute_floesize, &
             seaice_compute_floesize, &
             seaice_restart_floesize, &
             seaice_finalize_floesize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_floesize
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  Erin E Thomas
!> \date    26 September 2023
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_floesize(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_bootstrap_floesize!}}}

!***********************************************************************
!
!  routine seaice_init_floesize
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  Erin E Thomas
!> \date    26 September 2023
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_floesize(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_init_floesize!}}}

!***********************************************************************
!
!  routine seaice_precompute_floesize
!
!> \brief   Precompute MPAS-Seaice analysis member
!> \author  Erin E Thomas
!> \date    26 September 2023
!> \details
!>
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_floesize(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_precompute_floesize!}}}

!***********************************************************************
!
!  routine seaice_compute_floesize
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Erin E Thomas
!> \date    26 September 2023
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member. It computes ice and snow layer
!>  floesize from ice area and volume
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_floesize(domain, instance, timeLevel, err)!{{{

     use seaice_constants, only: &
          seaicePuny

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (block_type), pointer :: &
           blockPtr

      type (mpas_pool_type), pointer :: &
           tracersPool, &
           tracersAggregatePool, &
           floesizeAMPool

      real(kind=RKIND), dimension(:), pointer :: &
           iceAreaCell, &
           iceVolumeCell, &
           snowVolumeCell, &
           iceThicknessCell, &
           snowThicknessCell

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           iceAreaCategory, &
           iceVolumeCategory, &
           snowVolumeCategory

      real(kind=RKIND), dimension(:,:), pointer :: &
           iceThicknessCategory, &
           snowThicknessCategory

      integer, pointer :: &
           nCellsSolve, &
           nCategories

      integer :: &
           iCell, &
           iCategory

      err = 0

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_subpool(blockPtr % structs, 'floesizeAM', floesizeAMPool)
         call MPAS_pool_get_subpool(blockPtr % structs, 'tracers', tracersPool)
         call MPAS_pool_get_subpool(blockPtr % structs, 'tracers_aggregate', tracersAggregatePool)

         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nCellsSolve', nCellsSolve)
         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nCategories', nCategories)

         call MPAS_pool_get_array(floesizeAMPool, 'iceThicknessCell',      iceThicknessCell)
         call MPAS_pool_get_array(floesizeAMPool, 'iceThicknessCategory',  iceThicknessCategory)
         call MPAS_pool_get_array(floesizeAMPool, 'snowThicknessCell',     snowThicknessCell)
         call MPAS_pool_get_array(floesizeAMPool, 'snowThicknessCategory', snowThicknessCategory)

         call MPAS_pool_get_array(tracersPool, 'iceAreaCategory', iceAreaCategory, 1)
         call MPAS_pool_get_array(tracersPool, 'iceVolumeCategory', iceVolumeCategory, 1)
         call MPAS_pool_get_array(tracersPool, 'snowVolumeCategory', snowVolumeCategory, 1)

         call MPAS_pool_get_array(tracersAggregatePool, 'iceAreaCell', iceAreaCell)
         call MPAS_pool_get_array(tracersAggregatePool, 'iceVolumeCell', iceVolumeCell)
         call MPAS_pool_get_array(tracersAggregatePool, 'snowVolumeCell', snowVolumeCell)

         ! initialize
         iceThicknessCell(:) = 0.0_RKIND
         iceThicknessCategory(:,:) = 0.0_RKIND
         snowThicknessCell(:) = 0.0_RKIND
         snowThicknessCategory(:,:) = 0.0_RKIND

         ! compute floesize
         do iCell = 1, nCellsSolve

            ! aggregate quantities
            if (iceAreaCell(iCell) > seaicePuny) then

               iceThicknessCell(iCell)  = iceVolumeCell(iCell)  / iceAreaCell(iCell)
               snowThicknessCell(iCell) = snowVolumeCell(iCell) / iceAreaCell(iCell)

            endif ! iceAreaCell

            ! category quantities
            do iCategory = 1, nCategories

               ! check if ice present
               if (iceAreaCategory(1,iCategory,iCell) > seaicePuny) then

                  iceThicknessCategory(iCategory,iCell)  = iceVolumeCategory(1,iCategory,iCell)  / iceAreaCategory(1,iCategory,iCell)
                  snowThicknessCategory(iCategory,iCell) = snowVolumeCategory(1,iCategory,iCell) / iceAreaCategory(1,iCategory,iCell)

               endif ! iceAreaCategory

            enddo ! iCategory
         enddo ! iCell

         blockPtr => blockPtr % next
      enddo

   end subroutine seaice_compute_floesize!}}}

!***********************************************************************
!
!  routine seaice_restart_floesize
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  Erin E Thomas
!> \date    26 September 2023
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_floesize(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_restart_floesize!}}}

!***********************************************************************
!
!  routine seaice_finalize_floesize
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  Erin E Thomas
!> \date    26 September 2023
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_floesize(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_finalize_floesize!}}}

!-----------------------------------------------------------------------

end module seaice_floesize
