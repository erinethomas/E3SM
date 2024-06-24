module ww3_cpl_indices
  
  use seq_flds_mod
  use mct_mod

  implicit none

  SAVE
  public                               ! By default make data private

  integer :: index_x2w_Sa_u     
  integer :: index_x2w_Sa_v     
  integer :: index_x2w_Sa_tbot  
  integer :: index_x2w_Si_ifrac
  integer :: index_x2w_Si_ithick
  integer :: index_x2w_Si_ifloe 
  integer :: index_x2w_So_t     
  integer :: index_x2w_So_u     
  integer :: index_x2w_So_v     
  integer :: index_x2w_So_bldepth     
  integer :: index_x2w_So_ssh

  integer :: index_w2x_Sw_ustokes_wavenumber_1
  integer :: index_w2x_Sw_vstokes_wavenumber_1
  integer :: index_w2x_Sw_ustokes_wavenumber_2
  integer :: index_w2x_Sw_vstokes_wavenumber_2
  integer :: index_w2x_Sw_ustokes_wavenumber_3
  integer :: index_w2x_Sw_vstokes_wavenumber_3
  integer :: index_w2x_Sw_ustokes_wavenumber_4
  integer :: index_w2x_Sw_vstokes_wavenumber_4
  integer :: index_w2x_Sw_ustokes_wavenumber_5
  integer :: index_w2x_Sw_vstokes_wavenumber_5
  integer :: index_w2x_Sw_ustokes_wavenumber_6
  integer :: index_w2x_Sw_vstokes_wavenumber_6

  integer :: index_w2x_Sw_Hs
  integer :: index_w2x_Sw_Fp
  integer :: index_w2x_Sw_Dp

  integer,dimension(:),allocatable :: index_w2x_Sw_wavespec

contains

  subroutine ww3_cpl_indices_set( )

    use seq_flds_mod, only : wav_nfreq,wav_ocn_coup,wav_ice_coup 
    type(mct_aVect) :: w2x      ! temporary
    type(mct_aVect) :: x2w      ! temporary

    integer :: i
    character(len= 2) :: freqnum
    character(len=64) :: name

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2w, rList=seq_flds_x2w_fields, lsize=1)
    call mct_aVect_init(w2x, rList=seq_flds_w2x_fields, lsize=1)

    index_x2w_Sa_u       = mct_avect_indexra(x2w,'Sa_u')       ! Zonal wind at lowest level (this should probably be at 10m)
    index_x2w_Sa_v       = mct_avect_indexra(x2w,'Sa_v')       ! Meridional wind at lowest level (see above)
    index_x2w_Sa_tbot    = mct_avect_indexra(x2w,'Sa_tbot')    ! Temperature at lowest level
    index_x2w_Si_ifrac   = mct_avect_indexra(x2w,'Si_ifrac')   ! Fractional sea ice coverage 
    index_x2w_Si_ithick  = mct_avect_indexra(x2w,'Si_ithick')  ! Sea ice thickness
    index_x2w_Si_ifloe   = mct_avect_indexra(x2w,'Si_ifloe')   ! Sea ice floe size
    index_x2w_So_t       = mct_avect_indexra(x2w,'So_t')       ! Sea surface temperature
    index_x2w_So_u       = mct_avect_indexra(x2w,'So_u')       ! Zonal sea surface water velocity
    index_x2w_So_v       = mct_avect_indexra(x2w,'So_v')       ! Meridional sea surface water velocity
    index_x2w_So_bldepth = mct_avect_indexra(x2w,'So_bldepth') ! Boundary layer depth
    index_x2w_So_ssh     = mct_avect_indexra(x2w,'So_ssh')     ! Sea surface height 

    if (wav_ocn_coup .eq. 'two' .or. wav_atm_coup .eq. 'two') then
       index_w2x_Sw_Charn     = mct_avect_indexra(w2x,'Sw_Charn') ! Charnock coeff accounting for the wave stress (Janssen 1989, 1991)
    endif
    if (wav_ocn_coup .eq. 'two') then
       index_w2x_Sw_Hs = mct_avect_indexra(w2x,'Sw_Hs') ! Significant wave height
       index_w2x_Sw_Fp = mct_avect_indexra(w2x,'Sw_Fp') ! Peak wave freqency  
       index_w2x_Sw_Dp = mct_avect_indexra(w2x,'Sw_Dp') ! Peak wave direction
       index_w2x_Sw_ustokes_wavenumber_1 = mct_avect_indexra(w2x,'Sw_ustokes_wavenumber_1') ! partitioned Stokes drift u 1
       index_w2x_Sw_vstokes_wavenumber_1 = mct_avect_indexra(w2x,'Sw_vstokes_wavenumber_1') ! partitioned Stokes drift v 1
       index_w2x_Sw_ustokes_wavenumber_2 = mct_avect_indexra(w2x,'Sw_ustokes_wavenumber_2') ! partitioned Stokes drift u 2
       index_w2x_Sw_vstokes_wavenumber_2 = mct_avect_indexra(w2x,'Sw_vstokes_wavenumber_2') ! partitioned Stokes drift v 2
       index_w2x_Sw_ustokes_wavenumber_3 = mct_avect_indexra(w2x,'Sw_ustokes_wavenumber_3') ! partitioned Stokes drift u 3
       index_w2x_Sw_vstokes_wavenumber_3 = mct_avect_indexra(w2x,'Sw_vstokes_wavenumber_3') ! partitioned Stokes drift v 3
       index_w2x_Sw_ustokes_wavenumber_4 = mct_avect_indexra(w2x,'Sw_ustokes_wavenumber_4') ! partitioned Stokes drift u 4
       index_w2x_Sw_vstokes_wavenumber_4 = mct_avect_indexra(w2x,'Sw_vstokes_wavenumber_4') ! partitioned Stokes drift v 4
       index_w2x_Sw_ustokes_wavenumber_5 = mct_avect_indexra(w2x,'Sw_ustokes_wavenumber_5') ! partitioned Stokes drift u 5
       index_w2x_Sw_vstokes_wavenumber_5 = mct_avect_indexra(w2x,'Sw_vstokes_wavenumber_5') ! partitioned Stokes drift v 5
       index_w2x_Sw_ustokes_wavenumber_6 = mct_avect_indexra(w2x,'Sw_ustokes_wavenumber_6') ! partitioned Stokes drift u 6
       index_w2x_Sw_vstokes_wavenumber_6 = mct_avect_indexra(w2x,'Sw_vstokes_wavenumber_6') ! partitioned Stokes drift v 6
    endif
    if (wav_ice_coup .eq. 'two') then
       allocate(index_w2x_Sw_wavespec(1:wav_nfreq))
       do i = 1,wav_nfreq
          write(freqnum,'(i2.2)') i
          name = 'Sw_wavespec' // freqnum
          index_w2x_Sw_wavespec(i) = mct_avect_indexra(w2x,trim(name)) ! full wave spectrum (fcn of frq)
       enddo
    endif
    call mct_aVect_clean(x2w)
    call mct_aVect_clean(w2x)

  end subroutine ww3_cpl_indices_set

end module ww3_cpl_indices
