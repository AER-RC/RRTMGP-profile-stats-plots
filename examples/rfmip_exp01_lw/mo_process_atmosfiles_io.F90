! This code is part
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! This module reads profiles needed for RRTMGP and related calculations assuming a certain
!   netCDF file layout

module mo_process_atmosfiles_io
  use mo_rte_kind,           only: wp
  use mo_optical_props,      only: ty_optical_props, ty_optical_props_arry, &
                                   ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_source_functions,   only: ty_source_func_lw
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_rrtmgp_util_reorder
  use mo_simple_netcdf,      only: read_field, read_string, var_exists, get_dim_size, &
                                   write_field, create_dim, create_var, dim_exists
  use netcdf
  implicit none
  private

  public :: read_atmos, write_atmos, is_lw, is_sw, &
            read_lw_bc, read_sw_bc, read_lw_rt, &
            write_fluxes, write_dir_fluxes,     &
            write_heating_rates, read_atmos_nCase, &
            write_fluxes_nCase, write_dir_fluxes_nCase,     &
            write_heating_rates_nCase, &
            write_spectral_disc, read_spectral_disc, &
            write_sw_surface_albedo, write_solar_zenith_angle, &
            write_lw_surface_emissivity, read_sfc_test_file, &
            write_optical_prop_values, read_optical_prop_values, &
            write_direction,     read_direction,  &
            write_lw_Planck_sources, read_lw_Planck_sources, &
            write_sw_solar_sources,  read_sw_solar_sources, &
            write_two_stream,    read_two_stream, &
            read_sources,        write_sources,   &
            write_gpt_fluxes,    read_gpt_fluxes

contains

  function reorder(arr_in) result(arr_out)
     real(wp), dimension(:,:,:), intent(in) :: arr_in
     real(wp), dimension(size(arr_in, 3), size(arr_in, 1), size(arr_in, 2)) :: arr_out
     call reorder123x312(arr_in, arr_out)
  end function reorder
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read profiles for all columns  -- T, p, and gas concentrations
  !   Allocation occurs on assignments (says the F2003 standard)
  !
  subroutine read_atmos(fileName,                          &
                        p_lay, t_lay, p_lev, t_lev,        &
                        gas_concs, col_dry)
    character(len=*),   intent(in   ) :: fileName
    real(wp), dimension(:,:), allocatable,                 &
                        intent(inout) :: p_lay, t_lay, p_lev, t_lev, col_dry
    type(ty_gas_concs), intent(inout) :: gas_concs
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, nlev, igas, ngases
    character(len=32), dimension(19) :: gas_names = (/ &
      'h2o    ', 'co2    ', 'o3     ', 'n2o    ', &
      'co     ', 'ch4    ', 'o2     ', 'n2     ', 'no2    ', &
      'cf4    ', 'ccl4   ', 'cfc11  ', 'cfc12  ', 'cfc22  ', &
      'hfc143a', 'hfc125 ', 'hfc23  ', 'hfc32  ', 'hfc134a'  &
    /)
    character(len=32), dimension(19) :: found_gases

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_atmos: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    nlev = get_dim_size(ncid, 'lev')
    if(nlev /= nlay+1) call stop_on_err("read_atmos: nlev should be nlay+1")

    !
    ! These lines assume that compilers follow the Fortran 2003 standard for
    !   allocating on assignment. This may require explicit compiler support
    !   e.g. -assume realloc_lhs flag for Intel
    !
    p_lay = read_field(ncid, 'p_lay', ncol, nlay)
    t_lay = read_field(ncid, 't_lay', ncol, nlay)
    p_lev = read_field(ncid, 'p_lev', ncol, nlev)
    t_lev = read_field(ncid, 't_lev', ncol, nlev)
    !
    ! Need to initialize gas_concentrations with the gas names we expect to find
    ! First, find gases in file
    !
    ngases = 0
    do igas = 1,size(gas_names)
      if (var_exists(ncid, 'vmr_' // trim(gas_names(igas)))) then
        ngases = ngases + 1
        found_gases(ngases) = trim(gas_names(igas))
      end if
    end do
    ! Next, initialize gas concentrations
    call stop_on_err(gas_concs%init(found_gases(1:ngases)))
    ! Finally, populate gas concentrations
    do igas = 1,size(gas_names)
      if (var_exists(ncid, 'vmr_' // trim(gas_names(igas)))) then
        call stop_on_err(gas_concs%set_vmr( &
          gas_names(igas), read_field(ncid, 'vmr_' // trim(gas_names(igas)), ncol, nlay) &
        ))
      end if
    end do

    ! col_dry has unchanged allocation status on return if the variable isn't present in the netCDF file
    if(var_exists(ncid, 'col_dry')) col_dry = read_field(ncid, 'col_dry', ncol, nlay)

    ncid = nf90_close(ncid)

  end subroutine read_atmos
  subroutine read_atmos_nCase(fileName,                          &
                        p_lay, t_lay, p_lev, t_lev,        &
                        gas_concs, col_dry)
    character(len=*),   intent(in   ) :: fileName
    real(wp), dimension(:,:,:), allocatable,                 &
                        intent(inout) :: p_lay, t_lay, p_lev, t_lev, col_dry
    type(ty_gas_concs), intent(inout),allocatable,dimension(:) :: gas_concs
    ! -------------------
    integer :: ncid, iCase
    integer :: ncol, nlay, nlev, nCase, ngases, igas
    real(wp),dimension(:,:,:),allocatable :: vmr_tmp
    character(len=32), dimension(19) :: gas_names = (/ &
      'h2o    ', 'co2    ', 'o3     ', 'n2o    ', &
      'co     ', 'ch4    ', 'o2     ', 'n2     ', 'no2    ', &
      'cf4    ', 'ccl4   ', 'cfc11  ', 'cfc12  ', 'cfc22  ', &
      'hfc143a', 'hfc125 ', 'hfc23  ', 'hfc32  ', 'hfc134a'  &
    /)
    character(len=32), dimension(19) :: found_gases

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) then
      call stop_on_err("read_atmos: can't find file " // trim(fileName))
    end if
    !
    ! Get dimension sizes from file
    !
    ncol  = get_dim_size(ncid, 'col')
    nlay  = get_dim_size(ncid, 'lay')
    nlev  = get_dim_size(ncid, 'lev')
    nCase = get_dim_size(ncid, 'record')
    if(nlev /= nlay+1) call stop_on_err("read_atmos: nlev should be nlay+1")
    !
    ! Read profiles from file.
    ! These lines assume that compilers follow the Fortran 2003 standard for
    !   allocating on assignment. This may require explicit compiler support
    !   e.g. -assume realloc_lhs flag for Intel
    !
    p_lay = read_field(ncid, 'p_lay', ncol, nlay, nCase)
    t_lay = read_field(ncid, 't_lay', ncol, nlay, nCase)
    p_lev = read_field(ncid, 'p_lev', ncol, nlev, nCase)
    t_lev = read_field(ncid, 't_lev', ncol, nlev, nCase)
    !
    ! Need to initialize gas_concentrations with all gas names we will use
    ! before we can populate the gas_concentrations object. First, find the
    ! gases present in the file...
    !
    ngases = 0
    do igas = 1,size(gas_names)
      if (var_exists(ncid, 'vmr_' // trim(gas_names(igas)))) then
        ngases = ngases + 1
        found_gases(ngases) = trim(gas_names(igas))
      end if
    end do
    if (ngases <= 0) then
      call stop_on_err('No gases found in file.')
    end if
    !
    ! ...Next, initialize gas concentrations object...
    !
    allocate(gas_concs(nCase))
    do iCase = 1,nCase
       call stop_on_err(gas_concs(iCase)%init(found_gases(1:ngases)))
    end do
    !
    ! ...Finally, set gas concentrations
    !
    do igas = 1,size(gas_names)
      if (var_exists(ncid, 'vmr_' // trim(gas_names(igas)))) then
        vmr_tmp = read_field(ncid, 'vmr_' // trim(gas_names(igas)), ncol, nlay, nCase)
        do iCase = 1,nCase
          call stop_on_err(gas_concs(iCase)%set_vmr( &
            gas_names(igas), vmr_tmp(:,:,iCase) &
          ))
        end do
      end if
    end do
    !
    ! Read col_dry from file
    ! col_dry has unchanged allocation status on return if the variable isn't present in the netCDF file
    !
    if(var_exists(ncid, 'col_dry')) col_dry = read_field(ncid, 'col_dry', ncol, nlay, nCase)
    !
    ! Close netcdf file
    !
    ncid = nf90_close(ncid)

  end subroutine read_atmos_nCase


  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write the atmospheric conditions that might be computed on the fly
  !
  subroutine write_atmos(fileName, t_lev, col_dry)
    character(len=*),         intent(in) :: fileName
    real(wp), dimension(:,:), intent(in) :: t_lev, col_dry

    integer :: ncid, ncol, nlev, nlay

    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_atmos: can't open file " // trim(fileName))
    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = get_dim_size(ncid, 'col')
    nlay  = get_dim_size(ncid, 'lay')
    nlev  = get_dim_size(ncid, 'lev')
    call create_var(ncid, "col_dry", ["col",  "lay"], [ncol, nlay])
    call create_var(ncid, "t_lev",   ["col",  "lev"], [ncol, nlev])
    call stop_on_err(write_field(ncid, "col_dry",  col_dry ))
    call stop_on_err(write_field(ncid, "t_lev",     t_lev ))

    ncid = nf90_close(ncid)
  end subroutine write_atmos
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Does this file contain variables needed to do SW calculations ?
  !
  function is_sw(fileName)
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_sw

    integer :: ncid

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("is_sw: can't find file " // trim(fileName))

    is_sw = var_exists(ncid, 'solar_zenith_angle')
    ncid = nf90_close(ncid)
  end function is_sw
  ! ----------------------
  function is_lw(fileName)
    character(len=*), intent(in   ) :: fileName
    logical                         :: is_lw

    is_lw = .not. is_sw(fileName)
  end function is_lw
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read LW boundary conditions for all columns
  !
  subroutine read_lw_bc(fileName, t_sfc, emis_sfc)
    character(len=*),                      intent(in   ) :: fileName
    real(wp), dimension(:),   allocatable, intent(inout) :: t_sfc
    real(wp), dimension(:,:), allocatable, intent(inout) :: emis_sfc
    ! -------------------
    integer :: ncid
    integer :: ncol, nband
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_lw_bc: can't find file " // trim(fileName))

    ncol  = get_dim_size(ncid, 'col')
    t_sfc     =  read_field(ncid, 't_sfc',           ncol)

    !
    ! In case the spectral dimension has been removed.
    !
    if(dim_exists(ncid, 'band')) then
      nband = get_dim_size(ncid, 'band')
      emis_sfc =  read_field(ncid, 'emis_sfc', nband, ncol)
    else
      emis_sfc =  reshape(read_field(ncid, 'emis_sfc', ncol), shape = [1,ncol])
    end if

    ncid = nf90_close(ncid)
  end subroutine read_lw_bc
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read LW radiative transfer parameters
  !
  subroutine read_lw_rt(fileName, n_quad_angles)
    character(len=*),                      intent(in   ) :: fileName
    integer,                               intent(  out) :: n_quad_angles
    ! -------------------
    integer :: ncid
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_lw_bc: can't find file " // trim(fileName))
    n_quad_angles  = get_dim_size(ncid, 'angle')
    ncid = nf90_close(ncid)
  end subroutine read_lw_rt
 !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read SW boundary conditions for all columns
  !
  subroutine read_sw_bc(fileName, sza, tsi, tsi_scaling, sfc_alb_dir, sfc_alb_dif)
    character(len=*),                      intent(in   ) :: fileName
    real(wp), dimension(:),   allocatable, intent(inout) :: sza, tsi
    real(wp), dimension(:,:), allocatable, intent(inout) :: sfc_alb_dir, sfc_alb_dif
    real(wp),                              intent(inout) :: tsi_scaling
    ! -------------------
    integer :: ncid
    integer :: ncol, nband
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_sw_bc: can't find file " // trim(fileName))

    ncol  = get_dim_size(ncid, 'col')
    nband = get_dim_size(ncid, 'band')

    sza         =  read_field(ncid, 'solar_zenith_angle',        ncol)
    tsi         =  read_field(ncid, 'total_solar_irradiance',    ncol)
    if(var_exists(ncid, 'band')) then
      sfc_alb_dir =  read_field(ncid, 'sfc_alb_direct',  nband, ncol)
      sfc_alb_dif =  read_field(ncid, 'sfc_alb_diffuse', nband, ncol)
    else
      sfc_alb_dir = reshape(read_field(ncid, 'sfc_alb_direct',  ncol), shape = [1,ncol])
      sfc_alb_dif = reshape(read_field(ncid, 'sfc_alb_diffuse', ncol), shape = [1,ncol])
    end if

    ! read tsi_scaling only if variable is present in the netCDF file
    if(var_exists(ncid, 'tsi_scaling')) tsi_scaling = read_field(ncid, 'tsi_scaling' )

    ncid = nf90_close(ncid)
  end subroutine read_sw_bc
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write broadband and by-band fluxes
  !
  subroutine write_fluxes(fileName, flux_up, flux_dn, flux_net, bnd_flux_up, bnd_flux_dn, bnd_flux_net)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:  ), intent(in) :: flux_up,     flux_dn,     flux_net
    real(wp), dimension(:,:,:), intent(in), optional :: bnd_flux_up, bnd_flux_dn, bnd_flux_net
    ! -------------------
    integer :: ncid
    integer :: ncol, nlev, nband
    real(wp), dimension(:,:,:), allocatable :: bnd_flux_out
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))
    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = get_dim_size(ncid, 'col')
    nlev  = get_dim_size(ncid, 'lev')
    nband = get_dim_size(ncid, 'band')

    call create_var(ncid,      "flux_up",          ["col",  "lev"],         [ncol, nlev])
    call create_var(ncid,      "flux_dn",          ["col",  "lev"],         [ncol, nlev])
    call create_var(ncid,      "flux_net",         ["col",  "lev"],         [ncol, nlev])
    if(present(bnd_flux_up )) call create_var(ncid, "band_flux_up",  ["band", "col ", "lev "], [nband, ncol, nlev])
    if(present(bnd_flux_dn )) call create_var(ncid, "band_flux_dn",  ["band", "col ", "lev "], [nband, ncol, nlev])
    if(present(bnd_flux_net)) call create_var(ncid, "band_flux_net", ["band", "col ", "lev "], [nband, ncol, nlev])

    call stop_on_err(write_field(ncid, "flux_up",  flux_up ))
    call stop_on_err(write_field(ncid, "flux_dn",  flux_dn ))
    call stop_on_err(write_field(ncid, "flux_net", flux_net))
    ! col,lay,bnd -> bnd,col,lay
    allocate(bnd_flux_out(nband, ncol, nlev))
    if(present(bnd_flux_up )) call stop_on_err(write_field(ncid, "band_flux_up",  reorder(bnd_flux_up )))
    if(present(bnd_flux_dn )) call stop_on_err(write_field(ncid, "band_flux_dn",  reorder(bnd_flux_dn )))
    if(present(bnd_flux_net)) call stop_on_err(write_field(ncid, "band_flux_net", reorder(bnd_flux_net)))
    deallocate(bnd_flux_out)

    ncid = nf90_close(ncid)
  end subroutine write_fluxes

  subroutine write_fluxes_nCase(fileName, flux_up, flux_dn, flux_net, bnd_flux_up, bnd_flux_dn, bnd_flux_net)
    character(len=*),             intent(in) :: fileName
    real(wp), dimension(:,:  ,:), intent(in) ::     flux_up,     flux_dn,     flux_net
    real(wp), dimension(:,:,:,:), optional, &
                                  intent(inout) :: bnd_flux_up, bnd_flux_dn, bnd_flux_net
    ! -------------------
    integer :: ncid
    integer :: ncol, nlev, nband, nCase, iCase
    real(wp),dimension(:,:,:,:),allocatable :: &
        bnd_flux_up2, bnd_flux_dn2, bnd_flux_net2
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))
    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = get_dim_size(ncid, 'col')
    nlev  = get_dim_size(ncid, 'lev')
    nband = get_dim_size(ncid, 'band')

    nCase = size(flux_dn,3)

    call create_var(ncid,      "flux_up",          ["col   ",  "lev   ", "record"],         [ncol, nlev, ncase])
    call create_var(ncid,      "flux_dn",          ["col   ",  "lev   ", "record"],         [ncol, nlev, ncase])
    call create_var(ncid,      "flux_net",         ["col   ",  "lev   ", "record"],         [ncol, nlev, ncase])
    if(present(bnd_flux_up )) call create_var(ncid, "band_flux_up",  ["band  ", "col   ", "lev   ", "record"], &
       [nband, ncol, nlev, nCase])
    if(present(bnd_flux_dn )) call create_var(ncid, "band_flux_dn",  ["band  ", "col   ", "lev   ", "record"], &
       [nband, ncol, nlev, ncase])
    if(present(bnd_flux_net)) call create_var(ncid, "band_flux_net", ["band  ", "col   ", "lev   ", "record"], &
       [nband, ncol, nlev, ncase])

    call stop_on_err(write_field(ncid, "flux_up",  flux_up ))
    call stop_on_err(write_field(ncid, "flux_dn",  flux_dn ))
    call stop_on_err(write_field(ncid, "flux_net", flux_net))
    if(present(bnd_flux_up )) allocate(bnd_flux_up2(nband,  ncol, nlev, nCase))
    if(present(bnd_flux_dn )) allocate(bnd_flux_dn2(nband,  ncol, nlev, nCase))
    if(present(bnd_flux_net)) allocate(bnd_flux_net2(nband, ncol, nlev, nCase))
    ! col,lay,bnd -> bnd,col,lay
    do iCase=1,nCase
       if(present(bnd_flux_up )) bnd_flux_up2(:,:,:,iCase)  = reorder(bnd_flux_up(:,:,:,iCase))
       if(present(bnd_flux_dn )) bnd_flux_dn2(:,:,:,iCase)  = reorder(bnd_flux_dn(:,:,:,iCase))
       if(present(bnd_flux_net)) bnd_flux_net2(:,:,:,iCase) = reorder(bnd_flux_net(:,:,:,iCase))
    enddo
    if(present(bnd_flux_up )) call stop_on_err(write_field(ncid, "band_flux_up",  bnd_flux_up2))
    if(present(bnd_flux_dn )) call stop_on_err(write_field(ncid, "band_flux_dn",  bnd_flux_dn2))
    if(present(bnd_flux_net)) call stop_on_err(write_field(ncid, "band_flux_net", bnd_flux_net2))

    ncid = nf90_close(ncid)
  end subroutine write_fluxes_nCase
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write direct-beam fluxes
  !
  subroutine write_dir_fluxes(fileName, flux_dir, bnd_flux_dir)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:  ), intent(in) ::     flux_dir
    real(wp), dimension(:,:,:), optional, &
                                intent(in) :: bnd_flux_dir
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, nband
    real(wp), allocatable, dimension(:,:,:) :: bnd_flux_out
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_dir_fluxes: can't open file " // trim(fileName))

    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = get_dim_size(ncid, 'col')
    nlay  = get_dim_size(ncid, 'lay')
    nband = get_dim_size(ncid, 'band')

    call create_var(ncid,      "flux_dir_dn",         ["col",  "lev"],         [ncol, nlay+1])
    if(present(bnd_flux_dir)) call create_var(ncid, "band_flux_dir_dn", ["band", "col ", "lev "], [nband, ncol, nlay+1])

    call stop_on_err(write_field(ncid, "flux_dir_dn",  flux_dir))
    if(present(bnd_flux_dir)) call stop_on_err(write_field(ncid, "band_flux_dir_dn",  reorder(bnd_flux_dir)))

    ncid = nf90_close(ncid)
  end subroutine write_dir_fluxes
  subroutine write_dir_fluxes_nCase(fileName, flux_dir, bnd_flux_dir)
    character(len=*),             intent(in) :: fileName
    real(wp), dimension(:,:  ,:), intent(in) ::     flux_dir
    real(wp), dimension(:,:,:,:), optional, &
                                  intent(in) :: bnd_flux_dir
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, nband, nCase, icase
    real(wp), dimension(:,:,:), allocatable :: bnd_flux_out
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_dir_fluxes: can't open file " // trim(fileName))

    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = get_dim_size(ncid, 'col')
    nlay  = get_dim_size(ncid, 'lay')
    nband = get_dim_size(ncid, 'band')
    nCase = size(flux_dir,3)

    call create_var(ncid,      "flux_dir_dn",         ["col   ",  "lev   ", "record"],         [ncol, nlay+1, ncase])
    if(present(bnd_flux_dir)) call create_var(ncid, "band_flux_dir_dn", ["band  ", "col   ", "lev   ","record"], [nband, ncol, nlay+1,ncase])

    call stop_on_err(write_field(ncid, "flux_dir_dn",  flux_dir))
    do iCase=1,nCase
       if(present(bnd_flux_dir)) call stop_on_err(write_field(ncid, "band_flux_dir_dn", reorder(bnd_flux_dir(:,:,:,iCase))))
    enddo

    ncid = nf90_close(ncid)
  end subroutine write_dir_fluxes_nCase


  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write heating rates (broadband, by-band)
  !
  subroutine write_heating_rates(fileName, heating_rate, bnd_heating_rate)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:  ), intent(in) ::     heating_rate
    real(wp), dimension(:,:,:), intent(in) :: bnd_heating_rate
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, nband
    real(wp), dimension(:,:,:), allocatable :: bnd_heating_out
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_heating_rates: can't open file " // trim(fileName))

    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = get_dim_size(ncid, 'col')
    nlay  = get_dim_size(ncid, 'lay')
    nband = get_dim_size(ncid, 'band')

    call create_var(ncid,      "heating_rate",          ["col", "lay"],         [ncol, nlay])
    call create_var(ncid, "band_heating_rate", ["band", "col ", "lay "], [nband, ncol, nlay])

    call stop_on_err(write_field(ncid,     "heating_rate",                     heating_rate))
    call stop_on_err(write_field(ncid, "band_heating_rate", reorder(bnd_heating_rate)))

    ncid = nf90_close(ncid)
  end subroutine write_heating_rates
  subroutine write_heating_rates_nCase(fileName, heating_rate, bnd_heating_rate)
    character(len=*),             intent(in) :: fileName
    real(wp), dimension(:,:  ,:), intent(in) ::     heating_rate
    real(wp), dimension(:,:,:,:), optional, &
                                  intent(in) :: bnd_heating_rate
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, nband, nCase, icase
    real(wp), dimension(:,:,:), allocatable :: bnd_heating_out
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_heating_rates: can't open file " // trim(fileName))

    !
    ! At present these dimension sizes aren't used
    !   We could certainly check the array sizes against these dimension sizes
    !
    ncol  = get_dim_size(ncid, 'col')
    nlay  = get_dim_size(ncid, 'lay')
    nband = get_dim_size(ncid, 'band')
    nCase = size(heating_rate,3)

    call create_var(ncid,      "heating_rate",          ["col   ", "lay   ", "record"],         [ncol, nlay, ncase])
    call stop_on_err(write_field(ncid,     "heating_rate",                     heating_rate))

    if(present(bnd_heating_rate)) then
      call create_var(ncid, "band_heating_rate", ["band  ", "col   ", "lay   ","record"], [nband, ncol, nlay,ncase])
      do iCase=1,nCase
         call stop_on_err(write_field(ncid, "band_heating_rate", reorder(bnd_heating_rate(:,:,:,iCase))))
      enddo
    end if

    ncid = nf90_close(ncid)
  end subroutine write_heating_rates_nCase

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write spectral discretization
  !
  subroutine write_spectral_disc(fileName, spectral_disc)
    character(len=*),        intent(in) :: fileName
    class(ty_optical_props), intent(in) :: spectral_disc

    integer :: ncid
    integer :: nband

    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_spectral_disc: can't open file " // trim(fileName))

    nband = spectral_disc%get_nband()
    call create_dim(ncid, 'band', nband)
    call create_dim(ncid, "pair", 2)

    call create_var(ncid, "band_lims_wvn", ["pair", "band"], [2, nband])
    call stop_on_err(write_field(ncid, "band_lims_wvn", spectral_disc%get_band_lims_wavenumber()))
    call create_var(ncid, "band_lims_gpt", ["pair", "band"], [2, nband], NF90_INT)
    call stop_on_err(write_field(ncid, "band_lims_gpt", spectral_disc%get_band_lims_gpoint()))

    ncid = nf90_close(ncid)
  end subroutine write_spectral_disc
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Read spectral discretization
  !
  subroutine read_spectral_disc(fileName, spectral_disc)
    character(len=*),       intent(in   ) :: fileName
    class(ty_optical_props), intent(inout) :: spectral_disc

    integer :: ncid
    integer :: nband
    integer,  dimension(:,:), allocatable :: band_lims_gpt
    real(wp), dimension(:,:), allocatable :: band_lims_wvn

    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_spectral_disc: can't open file " // trim(fileName))

    nband = get_dim_size(ncid, 'band')
    if (get_dim_size(ncid, 'pair') /= 2) &
      call stop_on_err("read_spectral_disc: pair dimension not 2 in file "//trim(fileName) )

    band_lims_wvn = read_field(ncid, 'band_lims_wvn', 2, nband)
    band_lims_gpt = read_field(ncid, 'band_lims_gpt', 2, nband)
    call stop_on_err(spectral_disc%init(band_lims_wvn, band_lims_gpt, read_string(ncid, 'name', 32)))

    ncid = nf90_close(ncid)
  end subroutine read_spectral_disc

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write shortwave direct and diffuse albedo
  !
  subroutine write_sw_surface_albedo(fileName, sfc_alb_direct, sfc_alb_diffuse)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:), intent(in) :: sfc_alb_direct, sfc_alb_diffuse ! Dimensions (nband,ncol)

    integer :: ncid
    integer :: ncol, nband

    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_sw_surface_albedo: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nband = get_dim_size(ncid, 'band')

    call create_var(ncid, "sfc_alb_direct", ["band", "col "], [nband, ncol])
    call stop_on_err(write_field(ncid, "sfc_alb_direct", sfc_alb_direct))

    call create_var(ncid, "sfc_alb_diffuse", ["band", "col "], [nband, ncol])
    call stop_on_err(write_field(ncid, "sfc_alb_diffuse", sfc_alb_diffuse))

    ncid = nf90_close(ncid)
  end subroutine write_sw_surface_albedo
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write solar zenith angles
  !
  subroutine write_solar_zenith_angle(fileName, solar_zenith_angle)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:), intent(in) :: solar_zenith_angle ! Dimensions (ncol)

    integer :: ncid
    integer :: ncol

    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_solar_zenith_angle: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')

    call create_var(ncid, "solar_zenith_angle", ["col"], [ncol])
    call stop_on_err(write_field(ncid, "solar_zenith_angle", solar_zenith_angle))

    ncid = nf90_close(ncid)
  end subroutine write_solar_zenith_angle
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Write longwave emissivity
  !
  subroutine write_lw_surface_emissivity(fileName, emis_sfc)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:), intent(in) :: emis_sfc ! Dimensions (ncol,nband)

    integer :: ncid
    integer :: ncol, nband

    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_lw_surface_emissivity: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nband = get_dim_size(ncid, 'band')

    call create_var(ncid, "emis_sfc", ["band", "col "], [nband, ncol])
    call stop_on_err(write_field(ncid, "emis_sfc", emis_sfc))

    ncid = nf90_close(ncid)
  end subroutine write_lw_surface_emissivity
 !--------------------------------------------------------------------------------------------------------------------
 !
 ! Read surface SW albedo and LW emissivity spectra from the surface test file
 !
  subroutine read_sfc_test_file(fileName, sfc_alb, sfc_emis)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:), allocatable, intent(inout) :: sfc_alb  ! Dimensions (nband,nspectra)
    real(wp), dimension(:,:), allocatable, intent(inout) :: sfc_emis ! Dimensions (nband,nspectra)

    integer :: ncid
    integer :: nswband, nlwband, nspectra

    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_sfc_test_file: can't open file " // trim(fileName))

    if (.not. var_exists(ncid, 'SW_albedo')) &
      call stop_on_err("read_sfc_test_file: file " //trim(fileName) // " doesn't contain SW_albedo field.")
    nswband = get_dim_size(ncid, 'swband')
    nspectra = get_dim_size(ncid, 'spectra')
     ! allocate on assignment
    sfc_alb = read_field(ncid, 'SW_albedo', nswband, nspectra)

    if (.not. var_exists(ncid, 'LW_emissivity')) &
      call stop_on_err("read_sfc_test_file: file " //trim(fileName) // " doesn't contain LW_emissivity field.")
    nlwband = get_dim_size(ncid, 'lwband')
    sfc_emis = read_field(ncid, 'LW_emissivity', nlwband, nspectra)

    ncid = nf90_close(ncid)
  end subroutine read_sfc_test_file
 !--------------------------------------------------------------------------------------------------------------------
  !
  ! Paired procedures for reading and writing intermediate results used in unit testing
  !
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Optical properties (tau, ssa and g or p if provided)
  !
  subroutine write_optical_prop_values(fileName, opt_props)
    character(len=*),             intent(in) :: fileName
    class(ty_optical_props_arry), intent(in) :: opt_props
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, ngpt, nmom, nband
    integer :: dimLengths(3)
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_optical_prop_values: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    ngpt = opt_props%get_ngpt()
    nband = opt_props%get_nband()

    !
    ! Spectral discretization
    !
    call create_dim(ncid, "gpt", ngpt)
    call create_dim(ncid, 'band', nband)
    call create_dim(ncid, "pair", 2)

    call create_var(ncid, "band_lims_wvn", ["pair", "band"], [2, nband])
    call stop_on_err(write_field(ncid, "band_lims_wvn", opt_props%get_band_lims_wavenumber()))
    call create_var(ncid, "band_lims_gpt", ["pair", "band"], [2, nband], NF90_INT)
    call stop_on_err(write_field(ncid, "band_lims_gpt", opt_props%get_band_lims_gpoint()))

    !
    ! Values
    !
    call create_var(ncid, "tau", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
    call stop_on_err(write_field(ncid, "tau", opt_props%tau))

    select type (opt_props)
      class is (ty_optical_props_2str) ! two-stream calculation
        call create_var(ncid, "ssa", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
        call stop_on_err(write_field(ncid, "ssa", opt_props%ssa))

        call create_var(ncid, "g", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
        call stop_on_err(write_field(ncid, "g", opt_props%g))

      class is (ty_optical_props_nstr) ! n-stream calculation
        call create_var(ncid, "ssa", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
        call stop_on_err(write_field(ncid, "ssa", opt_props%ssa))

        nmom = size(opt_props%p, 1)
        call create_dim(ncid, "mom", nmom)
        call create_var(ncid, "p", ["mom", "col", "lay", "gpt"], [nmom, ncol, nlay, ngpt])
        call stop_on_err(write_field(ncid, "p", opt_props%p))
    end select

    ncid = nf90_close(ncid)

  end subroutine write_optical_prop_values
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_optical_prop_values(fileName, opt_props)
    character(len=*),                          intent(in ) :: fileName
    class(ty_optical_props_arry), allocatable, intent(out) :: opt_props
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, ngpt, nmom, nband
    real(wp), dimension(:,:), allocatable :: band_lims_wvn
    integer,  dimension(:,:), allocatable :: band_lims_gpt
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_optical_prop_values: can't open file " // trim(fileName))

    if(.not. var_exists(ncid, 'tau')) &
      call stop_on_err("read_optical_prop_values: file " //trim(fileName) // " doesn't contain tau field.")

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    ngpt = get_dim_size(ncid, 'gpt')
    nband = get_dim_size(ncid, 'band')

    if(var_exists(ncid, 'p')) then
      nmom = get_dim_size(ncid, 'mom')
      allocate(ty_optical_props_nstr::opt_props)
    else if (var_exists(ncid, 'g')) then
      allocate(ty_optical_props_2str::opt_props)
    else
      allocate(ty_optical_props_1scl::opt_props)
    end if

    !
    ! Spectral discretization
    !
    if (get_dim_size(ncid, 'pair') /= 2) &
      call stop_on_err("read_optical_prop_values: pair dimension not 2 in file " // trim(fileName) )
    band_lims_wvn = read_field(ncid, 'band_lims_wvn', 2, nband)
    band_lims_gpt = read_field(ncid, 'band_lims_gpt', 2, nband)

    call stop_on_err(opt_props%init(band_lims_wvn, band_lims_gpt, read_string(ncid, 'name', 32)))
    select type (opt_props)
      class is (ty_optical_props_1scl)      ! No scattering
        call stop_on_err(opt_props%alloc_1scl(ncol, nlay))
      class is (ty_optical_props_2str) ! two-stream calculation
        call stop_on_err(opt_props%alloc_2str(ncol, nlay))
        opt_props%ssa = read_field(ncid, "ssa", ncol, nlay, ngpt)
        opt_props%g   = read_field(ncid, "g",   ncol, nlay, ngpt)
      class is (ty_optical_props_nstr) ! n-stream calculation
        call stop_on_err(opt_props%alloc_nstr(nmom, ncol, nlay))
        opt_props%ssa = read_field(ncid, "ssa",       ncol, nlay, ngpt)
        opt_props%p   = read_field(ncid, "p",   nmom, ncol, nlay, ngpt)
    end select
    opt_props%tau = read_field(ncid, "tau", ncol, nlay, ngpt)

    ncid = nf90_close(ncid)
  end subroutine read_optical_prop_values
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Which direction is up? Stored as a global attribute.
  !
  subroutine write_direction(fileName, top_at_1)
    character(len=*),           intent(in) :: fileName
    logical,                    intent(in) :: top_at_1
    ! -------------------
    integer :: ncid, status
    !
    ! Vertical ordering, stored as a global attribute
    !
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_direction: can't open file " // trim(fileName))

    if(nf90_redef(ncid) /= NF90_NOERR) &
      call stop_on_err("write_direction: can't put file into redefine mode")
    if(nf90_put_att(ncid, NF90_GLOBAL, "top_at_1", merge(1, 0, top_at_1)) /= NF90_NOERR) &
      call stop_on_err("write_direction: can't write attribute top_at_1" )
    if(nf90_enddef(ncid) /= NF90_NOERR) &
      call stop_on_err("write_direction: can't end redefinition??")

    ncid = nf90_close(ncid)
  end subroutine write_direction
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_direction(fileName, top_at_1)
    character(len=*),           intent(in ) :: fileName
    logical,                    intent(out) :: top_at_1
     ! -------------------
    integer :: ncid, status
    integer :: top

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_direction: can't open file " // trim(fileName))

    status = nf90_get_att(ncid, NF90_GLOBAL, "top_at_1", top)
    top_at_1 = top == 1

    ncid = nf90_close(ncid)
  end subroutine read_direction
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Sources of upward and downward diffuse radiation, for eah layer and at the surface
  !
  subroutine write_sources(fileName, source_up, source_dn, source_sfc)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:,:), intent(in) :: source_up, source_dn
    real(wp), dimension(:,:  ), intent(in) :: source_sfc
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, ngpt
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_sources: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    ngpt = size(source_up, 3)
    call create_dim(ncid, "gpt", ngpt)

    call create_var(ncid, "source_up", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
    call stop_on_err(write_field(ncid, "source_up", source_up))
    call create_var(ncid, "source_dn", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
    call stop_on_err(write_field(ncid, "source_dn", source_dn))

    call create_var(ncid, "source_sfc", ["col",       "gpt"], [ncol,       ngpt])
    call stop_on_err(write_field(ncid, "source_sfc", source_sfc))

    ncid = nf90_close(ncid)
  end subroutine write_sources
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_sources(fileName, source_up, source_dn, source_sfc)
    character(len=*),           intent(in ) :: fileName
    real(wp), dimension(:,:,:), allocatable, &
                                intent(out) :: source_up, source_dn
    real(wp), dimension(:,:  ), allocatable, &
                                intent(out) :: source_sfc
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, ngpt, nmom
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_sources: can't open file " // trim(fileName))

    if(.not. var_exists(ncid, 'source_up')) &
      call stop_on_err("read_sources: file " //trim(fileName) // " doesn't contain source_up field.")

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    ngpt = get_dim_size(ncid, 'gpt')

    source_up  = read_field(ncid, 'source_up',  ncol, nlay, ngpt)
    source_dn  = read_field(ncid, 'source_dn',  ncol, nlay, ngpt)
    source_sfc = read_field(ncid, 'source_sfc', ncol,       ngpt)

    ncid = nf90_close(ncid)

  end subroutine read_sources

  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Longwave sources at layer centers; edges in two directions; surface
  !   Also directionality since this will be needed for solution
  !
  subroutine write_lw_Planck_sources(fileName, sources)
    character(len=*),        intent(in) :: fileName
    type(ty_source_func_lw), intent(in) :: sources

    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, ngpt
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_lw_Planck_sources: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    if(any([sources%get_ncol(), sources%get_nlay()] /= [ncol, nlay])) &
      call stop_on_err("write_lw_Planck_sources: inconsistent sizes in file, sources")

    ngpt = sources%get_ngpt()
    call create_dim(ncid, "gpt", ngpt)

    call create_var(ncid, "lay_src", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
    call stop_on_err(write_field(ncid, "lay_src", sources%lay_source))

    call create_var(ncid, "lev_src_inc", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
    call stop_on_err(write_field(ncid, "lev_src_inc", sources%lev_source_inc))
    call create_var(ncid, "lev_src_dec", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
    call stop_on_err(write_field(ncid, "lev_src_dec", sources%lev_source_dec))

    call create_var(ncid, "sfc_src", ["col", "gpt"], [ncol, ngpt])
    call stop_on_err(write_field(ncid, "sfc_src", sources%sfc_source))

    ncid = nf90_close(ncid)
  end subroutine write_lw_Planck_sources
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_lw_Planck_sources(fileName, sources)
    character(len=*),        intent(in   ) :: fileName
    type(ty_source_func_lw), intent(inout) :: sources
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, ngpt, nmom, nband
    integer,  dimension(:,:), allocatable :: band_lims_gpt
    real(wp), dimension(:,:), allocatable :: band_lims_wvn
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_lw_Planck_sources: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    ngpt = get_dim_size(ncid, 'gpt')
    nband = get_dim_size(ncid, 'band')
    !
    ! Error checking
    !
    if (get_dim_size(ncid, 'pair') /= 2) &
      call stop_on_err("read_spectral_disc: pair dimension not 2 in file "//trim(fileName) )
    if(.not. var_exists(ncid, 'lay_src')) &
      call stop_on_err("read_lw_Planck_sources: file " //trim(fileName) // " doesn't contain lay_src field.")

    !
    ! Spectral discretization
    !
    band_lims_wvn = read_field(ncid, 'band_lims_wvn', 2, nband)
    band_lims_gpt = read_field(ncid, 'band_lims_gpt', 2, nband)
    call stop_on_err(sources%init(band_lims_wvn, band_lims_gpt, read_string(ncid, 'name', 32)))
    call stop_on_err(sources%alloc(ncol, nlay))

    sources%lay_source     = read_field(ncid, 'lay_src',     ncol, nlay, ngpt)
    sources%lev_source_inc = read_field(ncid, 'lev_src_inc', ncol, nlay, ngpt)
    sources%lev_source_dec = read_field(ncid, 'lev_src_dec', ncol, nlay, ngpt)
    sources%sfc_source     = read_field(ncid, 'sfc_src',     ncol,       ngpt)

    ncid = nf90_close(ncid)

  end subroutine read_lw_Planck_sources
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Shortwave source at TOA
  !   Also directionality since this will be needed for solution
  !
  subroutine write_sw_solar_sources(fileName, toa_src)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:  ), intent(in) :: toa_src
    ! -------------------
    integer :: ncid
    integer :: ncol, ngpt
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_sw_solar_sources: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    ngpt = size(toa_src, 2)
    call create_dim(ncid, "gpt", ngpt)

    call create_var(ncid, "toa_src", ["col", "gpt"], [ncol, ngpt])
    call stop_on_err(write_field(ncid, "toa_src", toa_src))

    ncid = nf90_close(ncid)
  end subroutine write_sw_solar_sources
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_sw_solar_sources(fileName, toa_src)
    character(len=*),           intent(in ) :: fileName
    real(wp), dimension(:,:), allocatable, &
                                intent(out) :: toa_src
    ! -------------------
    integer :: ncid, status
    integer :: ncol, nlay, ngpt, nmom
    integer :: top
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_sw_solar_sources: can't open file " // trim(fileName))

    if(.not. var_exists(ncid, 'toa_src')) &
      call stop_on_err("read_sw_solar_sources: file " //trim(fileName) // " doesn't contain toa_src field.")

    ncol = get_dim_size(ncid, 'col')
    ngpt = get_dim_size(ncid, 'gpt')

    toa_src     = read_field(ncid, 'toa_src', ncol, ngpt)

    ncid = nf90_close(ncid)
  end subroutine read_sw_solar_sources
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Two-stream results: reflection and transmission for diffuse and direct radiation; also extinction
  !
  subroutine write_two_stream(fileName, Rdif, Tdif, Rdir, Tdir, Tnoscat)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:,:), intent(in) :: Rdif, Tdif
    real(wp), dimension(:,:,:), optional, &
                                intent(in) :: Rdir, Tdir, Tnoscat
    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, ngpt

    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_two_stream: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    ngpt = get_dim_size(ncid, 'gpt')

    call create_var(ncid, "Rdif", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
    call stop_on_err(write_field(ncid, "Rdif", Rdif))
    call create_var(ncid, "Tdif", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
    call stop_on_err(write_field(ncid, "Tdif", Tdif))
    if(present(Rdir)) then
      call create_var(ncid, "Rdir", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
      call stop_on_err(write_field(ncid, "Rdir", Rdir))
    end if
    if(present(Tdir)) then
      call create_var(ncid, "Tdir", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
      call stop_on_err(write_field(ncid, "Tdir", Tdir))
    end if
    if(present(Tnoscat)) then
      call create_var(ncid, "Tnoscat", ["col", "lay", "gpt"], [ncol, nlay, ngpt])
      call stop_on_err(write_field(ncid, "Tnoscat", Tnoscat))
    end if

    ncid = nf90_close(ncid)

  end subroutine write_two_stream
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_two_stream(fileName, Rdif, Tdif, Rdir, Tdir, Tnoscat)
    character(len=*),           intent(in ) :: fileName
    real(wp), dimension(:,:,:), allocatable, &
                                intent(out) :: Rdif, Tdif
    real(wp), dimension(:,:,:), allocatable, optional, &
                                intent(out) :: Rdir, Tdir, Tnoscat

    ! -------------------
    integer :: ncid
    integer :: ncol, nlay, ngpt
    ! -------------------

    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_two_stream: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'lay')
    ngpt = get_dim_size(ncid, 'gpt')

    Rdif    = read_field(ncid, 'Rdif',    ncol, nlay, ngpt)
    Tdif    = read_field(ncid, 'Tdif',    ncol, nlay, ngpt)
    if(present(Rdir)) Rdir       = read_field(ncid, 'Rdir',    ncol, nlay, ngpt)
    if(present(Tdir)) Tdir       = read_field(ncid, 'Tdir',    ncol, nlay, ngpt)
    if(present(Tnoscat)) Tnoscat = read_field(ncid, 'Tnoscat', ncol, nlay, ngpt)

    ncid = nf90_close(ncid)
  end subroutine read_two_stream
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! g-point fluxes
  !
  subroutine write_gpt_fluxes(fileName, gpt_flux_up, gpt_flux_dn, gpt_flux_dn_dir)
    character(len=*),           intent(in) :: fileName
    real(wp), dimension(:,:,:), intent(in) :: gpt_flux_up, gpt_flux_dn
    real(wp), dimension(:,:,:), optional, &
                                intent(in) :: gpt_flux_dn_dir
    ! -------------------
    integer :: ncid
    integer :: ncol, nlev, ngpt
    ! -------------------
    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_gpt_fluxes: can't open file " // trim(fileName))

    ncol = size(gpt_flux_up, 1)
    nlev = size(gpt_flux_up, 2)
    ngpt = size(gpt_flux_up, 3)
    call create_dim(ncid, "gpt", ngpt)

    call create_var(ncid, "gpt_flux_up", ["col", "lev", "gpt"], [ncol, nlev, ngpt])
    call stop_on_err(write_field(ncid, "gpt_flux_up", gpt_flux_up))
    call create_var(ncid, "gpt_flux_dn", ["col", "lev", "gpt"], [ncol, nlev, ngpt])
    call stop_on_err(write_field(ncid, "gpt_flux_dn", gpt_flux_dn))

    if(present(gpt_flux_dn_dir)) then
      call create_var(ncid, "gpt_flux_dn_dir", ["col", "lev", "gpt"], [ncol, nlev, ngpt])
      call stop_on_err(write_field(ncid, "gpt_flux_dn_dir", gpt_flux_dn_dir))
    end if
    ncid = nf90_close(ncid)
  end subroutine write_gpt_fluxes
  !--------------------------------------------------------------------------------------------------------------------
  subroutine read_gpt_fluxes(fileName, gpt_flux_up, gpt_flux_dn, gpt_flux_dn_dir)
    character(len=*),           intent(in ) :: fileName
    real(wp), dimension(:,:,:), allocatable, &
                                intent(out) :: gpt_flux_up, gpt_flux_dn
    real(wp), dimension(:,:,:), allocatable, optional, &
                                intent(out) :: gpt_flux_dn_dir
    ! -------------------
    integer :: ncid
    integer :: ncol, nlev, ngpt
    ! -------------------
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_gpt_fluxes: can't open file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlev = get_dim_size(ncid, 'lev')
    ngpt = get_dim_size(ncid, 'gpt')

    gpt_flux_up    = read_field(ncid, 'gpt_flux_up', ncol, nlev, ngpt)
    gpt_flux_dn    = read_field(ncid, 'gpt_flux_dn', ncol, nlev, ngpt)
    if(present(gpt_flux_dn_dir)) &
      gpt_flux_dn_dir = read_field(ncid, 'gpt_flux_dn_dir', ncol, nlev, ngpt)

    ncid = nf90_close(ncid)
  end subroutine read_gpt_fluxes
  !--------------------------------------------------------------------------------------------------------------------
  subroutine stop_on_err(msg)
    !
    ! Print error message and stop
    !
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg
    if(len_trim(msg) > 0) then
      write(error_unit,*) trim(msg)
      stop
    end if
  end subroutine
  !--------------------------------------------------------------------------------------------------------------------
end module mo_process_atmosfiles_io
