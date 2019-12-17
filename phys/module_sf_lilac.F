!-------------------------------------------------------------------------
module module_create_gindex
  implicit none
!========================================================================
contains
!========================================================================
  subroutine create_gindex(ide, jde, its, ite, jts, jte, gindex)
    ! Create a gindex array on each task. This gives the list of global indices owned by
    ! each processor, on the mass point grid.
    integer, intent(in) :: ide  ! domain end index, i
    integer, intent(in) :: jde  ! domain end index, j
    integer, intent(in) :: its  ! task start index, i
    integer, intent(in) :: ite  ! task end index, i
    integer, intent(in) :: jts  ! task start index, j
    integer, intent(in) :: jte  ! task end index, j
    integer, allocatable, intent(out) :: gindex(:)

    integer :: ite_limited  ! task end index on the mass point grid, i
    integer :: jte_limited  ! task end index on the mass point grid, j
    integer :: num_points
    integer :: i, j, n

    ! The very last index in both row & column space is just used on the momentum grid.
    ! Here we are just working with the mass point grid, so we need to ignore that last
    ! index.
    ite_limited = min(ite, ide-1)
    jte_limited = min(jte, jde-1)

    num_points = ((ite_limited - its + 1) * (jte_limited - jts + 1))

    print *, "ide :", ide
    print *, "jde :", jde
    print *, "its :", its
    print *, "ite :", ite
    print *, "jts :", jts
    print *, "jte :", jte
    print *, "ite_limited :", ite_limited
    print *, "jte_limited :", jte_limited
    print *, "num_points :, " , num_points

    allocate(gindex(num_points))

    print *, "here 3"
    n = 0
    do j = jts, jte_limited
       do i = its, ite_limited
          n = n + 1
          print *, "about to set ", n , "of", num_points
          ! In the following, note that we use ide-1 rather than ide for the same reason
          ! that we need ite_limited: ide gives the domain end index on the momentum grid,
          ! but here we're just dealing with the mass point grid, which has one less point
          ! in each direction.
          gindex(n) = (j-1)*(ide-1) + i

          ! This should be added in later:
          !INTEGER , DIMENSION(total_computational_MPI_ranks) :: count
          ! INTEGER, INTENT(IN) :: total_computational_MPI_ranks
          !DO m = 1, total_computational_MPI_ranks
          ! IF (m==1) THEN
          !WRITE(UNIT=10,FMT='("MPI RANK = ",I6," OF ",I6," Processes")') m-1,total_computational_MPI_ranks
          !DO c = 1, count(m)
          !    WRITE(UNIT=10,FMT='(I8i)' gindex
          ! END DO
          !ENDIF
          !END DO

       end do
    end do
  end subroutine create_gindex
end module module_create_gindex

!-------------------------------------------------------------------------

module module_lilac_reshape
  implicit none
!========================================================================
contains
!========================================================================

  subroutine lilac_reshape (ide, jde, its, ite, jts, jte, var_2d, var_1d)

    ! input/output variables
    integer, intent(in)               :: ide  ! domain end index, i
    integer, intent(in)               :: jde  ! domain end index, j
    integer, intent(in)               :: its  ! task start index, i
    integer, intent(in)               :: ite  ! task end index, i
    integer, intent(in)               :: jts  ! task start index, j
    integer, intent(in)               :: jte  ! task end index, j
    integer, allocatable, intent(in)  :: var_2d (:,:)
    integer, allocatable, intent(out) :: var_1d(:)

    ! local variables
    integer :: ite_limited  ! task end index on the mass point grid, i
    integer :: jte_limited  ! task end index on the mass point grid, j
    integer :: num_points
    integer :: i, j, n

    num_points  = ((ite_limited - its + 1) * (jte_limited - jts + 1))
    ite_limited = min(ite, ide-1)
    jte_limited = min(jte, jde-1)


    allocate (var_1d (num_points))

    n= 0

    do j = jts, jte_limited
       do i = its, ite_limited
          n = n + 1
          print *, "about to set var1d of ",n,"to var_2d of" i, j
          var_1d(n) = var_2d(i, j)

       end do
    end do


  end subroutine lilac_reshape
end module module_lilac_reshape


!========================================================================
! WRF Lilac subroutine
!========================================================================

module module_sf_lilac

  implicit none
!========================================================================
contains
!========================================================================
    use shr_kind_mod , only : r8 => shr_kind_r8
    use mpi          , only : MPI_COMM_WORLD    , MPI_COMM_NULL , MPI_Init , MPI_FINALIZE , MPI_SUCCESS

    ! parameters
    real, parameter :: rair   = 287.0423  !gas constant for dry air      [J/kg/K]
    real, parameter :: cpair  = 1.00464e3 !specific heat of dry air      [J/kg/K]
!-------------------------------------------------------------------------
    subroutine export_to_lilac (ide, jde, its, ite, jts, jte, var_name , var_2d)
        ! This subroutine reshapes the 2d array to 1d for Lilac
        ! and set it in the lilac_cap.

        ! in/out variables
        integer,              intent(in) :: ids,ide, jds,jde, kds,kde,  &  ! d -> domain
        integer,              intent(in) :: ims,ime, jms,jme, kms,kme,  &  ! m -> memory
        integer,              intent(in) :: its,ite, jts,jte, kts,kte      ! t -> tile
        character(len=*),     intent(in) :: var_name
        integer, allocatable, intent(in) :: var_2d (:,:)

        ! local variables
        integer, allocatable :: var_1d (:)

        ! -------------------------------------------------------------------
        ! 1) convert the 2d to 1d array
        call lilac_reshape (ide, jde, its, ite, jts, jte, var_name , var_2d)
       
        ! 2) setting the value in lilac atmos cap
        call lilac_atmcap_atm2lnd(var_name    , var_2d)

        print *, "Variable ", var_name , "is set in lilac atmospheric cap!"

    end subroutine
    !-------------------------------------------------------------------------

    subroutine import_from_lilac (ide, jde, its, ite, jts, jte, var_name , var_1d, var_2d)
        ! This subroutine reshapes the 1d array to 2d from Lilac
        ! and send back to atmosphere

        ! in/out variables
        integer,              intent(in) :: ids,ide, jds,jde, kds,kde,  &  ! d -> domain
        integer,              intent(in) :: ims,ime, jms,jme, kms,kme,  &  ! m -> memory
        integer,              intent(in) :: its,ite, jts,jte, kts,kte      ! t -> tile
        character(len=*),     intent(in) :: var_name
        integer, allocatable, intent(in) :: var_2d (:,:)
        integer, allocatable, intent(in) :: var_1d (:)

    end subroutine

        ! NOTE : This subroutine is NOT being used now! 
    subroutine wrf_lilac_init ()
          ! Here everything goes that is basically pre initalization in atm_driver
          ! This is mostly for description of what has been done in atm_driver and
          ! how it is translated for WRF

          ! -----------------------------------------------------------------------------
          ! Initiallize MPI
          ! -----------------------------------------------------------------------------
          ! initialization of MPI is already done earlier in WRF
          ! ==> no need for it here

          ! get the masterproc for printing out stuff...

          ! either with MPI_COMM_RANK or define a subroutine that does that for you.....


          ! -----------------------------------------------------------------------------
          ! Read in namelist file ...
          ! -----------------------------------------------------------------------------
          ! We should decide about this later on what we do and do NOT need from
          ! namelist

          ! The following will read this on all processors - might want to do a read just on the
          ! master processor and broadcast in the future

          ! no need to read the mesh from atmosphere anymore based on discussion
          ! with Mariana on 2019-12-06

          !open(newunit=fileunit, status="old", file="atm_driver_in")
          !read(fileunit, atm_driver_input, iostat=ierr)
          !if (ierr > 0) then
          !   print *, 'Error on reading atm_driver_in'
          !   call MPI_ABORT(MPI_COMM_WORLD, ierr)
          !end if
          !close(fileunit)


          ! -----------------------------------------------------------------------------
          ! Read mesh file to get number of points (n_points)
          ! -----------------------------------------------------------------------------

          ! no need to read the mesh from atmosphere anymore based on discussion
          ! with Mariana on 2019-12-06

          !print *, "DEBUG: atm_mesh_file = ",trim(atm_mesh_file)
          !call read_netcdf_mesh(atm_mesh_file, nglobal)
          !if (mytask == 0 ) then
          !   print *, " atm_driver mesh file ",trim(atm_mesh_file)
          !   print *, " number of global points in mesh is:", nglobal
          !end if


          !-----------------------------------------------------------------------------
          ! atmosphere domain decomposition
          !-----------------------------------------------------------------------------

          !here we are figuring out the atmosphere domain decomposition....
          ! we are using the create_gindex 
          ! not confirmed API but:

          ! integer, allocatable :: gindex(:) 
          ! create_gindex(ide, jde, its, ite, jts, jte, gindex)

          !------------------------------------------------------------------------
          ! Initialize lilac
          !------------------------------------------------------------------------

          if (mytask == 0 ) then
             print *, " initializing lilac "
          end if
          ! Since Mariana said she is going to change the interface for this....
          !call lilac_init(atm_mesh_file, atm_global_index, atm_lons, atm_lats, &

          atm_global_index = gindex_atm
          !call lilac_init( atm_global_index, atm_lons, atm_lats, &
          !     atm_calendar, atm_timestep, &
          !     atm_start_year, atm_start_mon, atm_start_day, atm_start_secs, &
          !     atm_stop_year, atm_stop_mon, atm_stop_day, atm_stop_secs)


    end subroutine

    !========================================================================

    subroutine wrf_lilac_run (                                    &
                  restart_alarm_is_ringing, stop_alarm_is_ringing, & ! in : Alarms for restart and stop
                  dz8w   , t3d , u_phys , v_phys,                  &
                  qv3d   , p8w3d ,glw   ,                               & !  in: forcing
                  swvisdir, swvisdif, swnirdir, swnirdif,          & ! shortwave radiation...
                  precip_in , dt,  sr,                            & ! total input precipitation  (mm; accumlated precipitation within dt)
                  dx , xlat , xlong , topo_hgt ,                                &
                                                                 ! Detailed precipitation
                  dzs , num_soil_layers , &                            ! in : model configurations
                  xland , xice, xice_thres, &                  !  in: land ice categories
                  tsk, hfx, qfx, tsk, lh, grdflx, &            ! in/out : lsm equivalant
                  smstav, smstot, sfcrunof , udrunoff, albedo, & !in/out: lsm equivalant
                  smois, snowc, snowh, snow,  &
                  tslb, emiss, 

                  snowc,qsfc,rainbl,maxpatch,                     &
                  dtbl,xtime, dt,dzs,             &
                  smois,tslb,snow,canwat,                         &
                  chs,chs2,snowh,                            &
                  shdmin,shdmax,                                  &
                  acsnom,acsnow,                                  &

                  ids,ide, jds,jde, kds,kde, &
                  ims,ime, jms,jme, kms,kme, &
                  its,ite, jts,jte, kts,kte, &
                  sf_urban_physics, nlcat,   &
                  num_timestep, year, month, day, julian_day )   ! in : Time-related 

                  !qgh,gsw,swdown,               &
                  !ra_sw_physics,history_interval,glw, &
                  !ivgtyp,isltyp,vegfra,        &
                  !znt,z0, tmn           &  ! 

                                                                                                       ! input/output variables
          logical ,                                         intent (in)    :: restart_alarm_is_ringing ! Logical flag for if the alarm rings this timestep for restart [-]
          logical ,                                         intent (in)    :: stop_alarm_is_ringing    ! Logical flag for if the alarm rings this timestep for ending the run []
          real , dimension( ims:ime, kms:kme, jms:jme),     intent (in)    :: dz8w                     ! thickness of first atmospheric layer [m] atmospheric bottom level height above surface [m]
          real ,    dimension( ims:ime, kms:kme, jms:jme ), intent (in)    :: t3d                      ! 3D atmospheric temperature [k]
          real ,    dimension( ims:ime, kms:kme, jms:jme ), intent (in)    :: u_phy                    ! 3D U wind component [m/s]
          real ,    dimension( ims:ime, kms:kme, jms:jme ), intent (in)    :: v_phy                    ! 3D U wind component [m/s]
          real ,    dimension( ims:ime, kms:kme, jms:jme ), intent (in)    :: qv3d                     ! 3D water vapor mixing ratio [kg/kg_dry]
          real ,    dimension( ims:ime, kms:kme, jms:jme ), intent (in)    :: p8w3d                    ! 3D pressure [Pa]
          real ,    dimension( ims:ime,jms:jme),             intent(in)    :: glw                     ! longwave down at surface [w m-2]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: tsk                      ! surface radiative temperature [k]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: hfx                      ! sensible heat flux [w m-2]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: qfx                      ! latent heat flux [kg s-1 m-2]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: lh                       ! latent heat flux [w m-2]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: grdflx                   ! ground/snow heat flux [w m-2]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: smstav                   ! soil moisture availailability for evatranspiration [not used in noahmp]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: smstot                   ! total soil water [mm][not used]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: sfcrunoff                ! accumulated surface runoff [m]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: udrunoff                 ! accumulated sub-surface runoff [m]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: albedo                   ! total grid albedo []
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: emiss !  surface bulk emissivity
          real ,    dimension( ims:ime, 1:nsoil, jms:jme ), intent (inout) :: smois                    ! volumetric soil moisture [m3/m3]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: snowc                    ! snow cover fraction [-]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: snowh                    ! physical snow depth [m]
          real ,    dimension( ims:ime,          jms:jme ), intent (inout) :: snow                     ! snow water equivalent [mm]
          real ,    dimension( ims:ime, 1:nsoil, jms:jme ), intent (inout) :: tslb                     ! soil temperature [k]
          real ,    dimension( ims:ime,          jms:jme ), intent (in   ) :: sr                       ! frozen precipitation ratio [-]
          real ,    dimension( ims:ime,          jms:jme ), intent (in   ) :: precip_in                ! total input precipitation [mm]
          real,    dimension(1:nsum_soil_layers),           intent (in   ) :: dzs                      ! thickness of soil layers [m]
          integer,                                          intent (in   ) :: num_soil_layers          ! number of soil layers
          real,    dimension( ims:ime,          jms:jme ),  intent (in)    :: prec_in                  ! total input precipitation [mm]
          real,                                             intent (in   ) :: dt                       ! timestep [s]
          real,                                             intent (in   ) :: dx                       ! horizontal grid spacing [m]
          real,    dimension( ims:ime,          jms:jme ),  intent (in)    :: xlat                     ! latitude [rad]
          real,    dimension( ims:ime,          jms:jme ),  intent (in)    :: xlong                    ! longitutde [rad]
          real,    dimension( ims:ime,          jms:jme ),  intent (in)    :: topo_hgt                 ! Terrain height [m]
          real,    dimension( ims:ime,          jms:jme ),  intent (in)    :: xland                    ! land mask (
          real,    dimension( ims:ime,          jms:jme ),  intent (in)    :: xice                     ! fraction of seaice in a grid
          real,    dimension( ims:ime,          jms:jme ),  intent (in)    :: xice_thres               ! fraction of a grid determining seaice

           real,    dimension( ims:ime,          jms:jme ), intent(inout) ::qsfc      ! bulk surface mixing ratio

           real , dimension (ims:ime,jms:jme),          intent(in) :: swvisdir  ! vis direct beam solar rad onto surface (W/m**2)
           real , dimension (ims:ime,jms:jme),          intent(in) :: swvisdif  ! vis diffuse solar rad onto surface (W/m**2)
           real , dimension (ims:ime,jms:jme),          intent(in) :: swnirdir  ! nir direct beam solar rad onto surface (W/m**2)
           real , dimension (ims:ime,jms:jme),          intent(in) :: swnirdif  ! nir diffuse solar rad onto surface(W/m**2)




           real, dimension(ims:ime,jms:jme),intent(in)  ::  swvisdir 
           real, dimension(ims:ime,jms:jme),intent(in)  ::  swvisdif
           real, dimension(ims:ime,jms:jme),intent(in)  ::  swnirdir
           real, dimension(ims:ime,jms:jme),intent(in)  ::  swnirdif

           real, dimension(ims:ime,jms:jme),intent(out)  ::  alswvisdir
           real, dimension(ims:ime,jms:jme),intent(out)  ::  alswvisdif
           real, dimension(ims:ime,jms:jme),intent(out)  ::  alswnirdir
           real, dimension(ims:ime,jms:jme),intent(out)  ::  alswnirdif
           real, dimensions ( ), intent(in)  :: forc_txy ! atmospheric bottom level temperature [k] 



          integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,  &  ! d -> domain
                                           ims,ime, jms,jme, kms,kme,  &  ! m -> memory
                                           its,ite, jts,jte, kts,kte      ! t -> tile
          integer , intent (in) :: num_timestep  ! timestep number
          integer , intent (in) :: year      ! year YYYY
          integer , intent (in) :: month      ! month MM
          integer , intent (in) :: day      ! day DD
          integer , intent (in) :: julian_day      ! Julian day 


          ! local variables 
          real , dimension( ims:ime, jms:jme) :: zlvl ! model forcing height [m]
          real , dimension( ims:ime, jms:jme) :: q_ratio  ! water vapor mixing ratio [kg/kg_dry]
          !-------------------------------------------------------------------------

          
          ! based on noah-mp first time-step sea handling....
          if(itimestep == 1)then
          jloop: do j=jts,jte
              do i=its,ite
                 if((xland(i,j)-1.5) >= 0.) then    ! open water case
                    if(xice(i,j) == 1. .and. iprint) print *,' sea-ice at water point, i=',i,'j=',j
                    smstav(i,j) = 1.0
                    smstot(i,j) = 1.0
                    do k = 1, nsoil
                       smois(i,k,j) = 1.0
                       tslb(i,k,j)  = 273.16
                    enddo
                 else
                    if(xice(i,j) == 1.) then        ! sea-ice case
                       smstav(i,j) = 1.0
                       smstot(i,j) = 1.0
                       do k = 1, nsoil
                          smois(i,k,j) = 1.0
                       enddo
                    endif
                 endif
              enddo
           endif                            ! end of initialization over ocean

          



          ! dz8w = thickness of full levels --> zlvl = dz8w/2

          do j=jts,jte
            do i=its,ite
              zlvl      (i,j) = 0.5 * dz8w (i,km,j) ! only first layer
              forc_uxy  (i,j) =  u_phy (i,km,j)
              forc_vxy  (i,j) =  v_phy (i,km,j)
              ! calculate potential temperature
              forc_th   (i,j) = forc_t (i,j) * (forc_psrf(i,j)/forc_pbot(i,j)) ** (rair/cpair)     ! pot temp = temp * (p0/p1) ^ (r/cp)
              forc_pbotxy (i,j) = p8w3d(i,km,j)   ! atmosphere bottom level pressure [pa]
              forc_t    (i,j) = t3d (i,km,j)    ! atmosphere bottom level temperature [k]
              q_ratio   (i,j) = qv3d(i,km,j)/(1.0 + qv3d(i,1,j)) ! converting to mixing ration [kg/kg] -- only first atmos layer
              prec_in   (i,j) = prec_in(i,km,j)/dt
              snow_in   (i,j) = prec_in (i,km,j) * sr(i,j)     /dt  ! devide by timestep for flux 
              rain_in   (i,j) = prec_in (i,km,j) * (1-sr(i,j)) /dt
              forc_psrf (i,j) = p8w3d (i,km,j)
              ! We send all the rain to rainc
              rainc_in  (i,j) = rain_in(i,j)
              rainl_in  (i,j) = 0.0
              snowc_in  (i,j) = snow_in(i,j)
              snowl_in  (i,j) = 0.0
            enddo
          enddo

          call export_to_lilac (ide, jde, its, ite, jts, jte, "Sa_z"    , zlvl(:,:))
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Sa_topo" , topo_hgt(:,:))
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Sa_u"    , forc_uxy(:,:))
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Sa_v"    , forc_vxy(:,:))
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Sa_ptem"  , forc_th(:,:))
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Sa_pbot"  , forc_pbotxy(:,:))
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Sa_tbot"  , forc_t(:,:) )
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Sa_shum"  , q_ratio(:,:) )

          call export_to_lilac (ide, jde, its, ite, jts, jte, "Faxa_lwdn"  , glw(:,:) )
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Faxa_rainc" , rainc_in (:,:)  )
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Faxa_rainl" , rainl_in (:,:) )
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Faxa_snowc" , snowc_in(:,:) )
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Faxa_snowl" , snowl_in(:,:) )

          call export_to_lilac (ide, jde, its, ite, jts, jte, "Faxa_swndn" , swnirdir  (:,:) )
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Faxa_swvdr" , swvisdir  (:,:) )
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Faxa_swndf" , swnirdif  (:,:) )
          call export_to_lilac (ide, jde, its, ite, jts, jte, "Faxa_swvdf" , swvisdif  (:,:) )


          ! -----------------------------------------------------------------------------
          !  Run lilac (lilac_run)
          ! -----------------------------------------------------------------------------

          call lilac_run (restart_alarm_is_ringing=.true.,stop_alarm_is_ringing=.true.)

          ! -----------------------------------------------------------------------------
          !  import from lilac 
          ! -----------------------------------------------------------------------------

          !call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_lfrin",  )
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_t",  )    !actual temperature of land
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_tref",  ) !2m reference temperature
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_qref",  )
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_avsdr",  )
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_anidr",  )
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_avsdf",  )
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_anidf",  )
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_snowh",  )
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_u10",  )
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_fv",  )
          call import_from_lilac (ide, jde, its, ite, jts, jte, "Sl_ram1",  )

    end subroutine wrf_lilac_run


    ! allocate and fill in atm2lnd
    subroutine wrf_to_lilac ( &
                    dz8w   , forc_txy , forc_uxy , forc_qxy, &   ! zgcmxy     , t_phy , u_phy, v_phy, &
                    ps  ,  flwdsxy, &                            ! p8w , glw
                    swvisdir, swvisdif, swnirdir, swnirdif, &    ! shortwave radiation...
                    prec_in  , dt &                              ! total input precipitation  (mm; accumlated precipitation within dt)
                                                                 ! Detailed precipitation
                   ,dzs , n_soil , dx, &                            ! in : model configurations
                    xland , xice, xice_thres, &                  !  in: land ice categories

                    ! t3d, qv3d, p8w3d , sr , &                  !  in: forcing
                    ,qfx        ,lhxy        ,soiflx      ,qgh       &
                   ,gsw, swdown,ra_sw_physics & ! let's pass sw_radiation option....

                   ,history_interval ,flwdsxy    ,smstav      ,smstot      ,qsfxy     &
                   ,qdnxy      ,ivgtyp     ,isltyp      ,vegfra      ,albxy     &
                   ,znt        ,z0         ,tmn         ,xland       ,xice      &
                   ,emiss      ,snowc      ,qsfc        ,prec        ,maxpatch  &
                   ,num_soil_layers        ,dt          ,xtime      ,dtwrf ,dzs &
                   ,smois      ,tslb       ,snow        ,canwat      ,chs       &
                   ,chs2                                                        &
                   ,sh2o       ,snowh      ,forc_uxy    ,forc_vxy    ,shdmin    &
                   ,shdmax     ,acsnom     ,acsnow      ,dx          ,xlat      &
                   ,xlong,ht                                                    &
                   ,ids,ide, jds,jde, kds,kde                    &
                   ,ims,ime, jms,jme, kms,kme                    &
                   ,its,ite, jts,jte, kts,kte                    &
                   ,inest, sf_urban_physics,nlcat               )!&


                  real,  dimension (ims:ime, kms:kme,jms:jme), intent(in) :: dz8w     !atm bottom level height above surface (m)
                  real,  dimension (ims:ime, kms:kme,jms:jme), intent(in) :: forc_txy   !atm bottom level temperature (Kelvin)
                  real,  dimension (ims:ime, kms:kme,jms:jme), intent(in) :: forc_uxy   !atm bottom level zonal wind (m/s)
                  real,  dimension (ims:ime, kms:kme,jms:jme), intent(in) :: forc_vxy   !atm bottom level meridional wind (m/s)
                  real,  dimension (ims:ime, kms:kme,jms:jme), intent(in) :: forc_qxy   !atm bottom level specific humidity (kg/kg)
                  real,  dimension (ims:ime, kms:kme,jms:jme), intent(in) :: ps
                  real,  dimension (ims:ime,jms:jme)           intent(in) :: flwdsxy   ! downward longwave rad onto surface (W/m**2)
                  real , dimension (ims:ime,jms:jme),          intent(in) :: swvisdir  ! vis direct beam solar rad onto srf (W/m**2)
                  real , dimension (ims:ime,jms:jme),          intent(in) :: swvisdif  ! vis diffuse solar rad onto srf (W/m**2)
                  real , dimension (ims:ime,jms:jme),          intent(in) :: swnirdir  ! nir direct beam solar rad onto srf (W/m**2)
                  real , dimension (ims:ime,jms:jme),          intent(in) :: swnirdif  ! nir diffuse solar rad onto srf(W/m**2)




                  real , dimension (ims:ime,jms:jme),          intent(in) :: qfx  ! nir diffuse solar rad onto srf(W/m**2)
                  
                  ! Optional Detailed Precipitation Partitioning Inputs

                  !1d arrays for filling in lilac
                  real, dimension :: Sa_z    (:)
                  real, dimension :: Sa_topo (:)
                  real, dimension :: Sa_u    (:)
                  real, dimension :: Sa_v    (:)
                  real, dimension :: Sa_ptem (:)
                  real, dimension :: Sa_shum (:)
                  real, dimension :: Sa_pbot (:)
                  real, dimension :: Sa_tbot (:)
                  real, dimension :: Faxa_lwdn (:)
                  real, dimension :: Faxa_rainc (:)
                  real, dimension :: Faxa_rainl (:)
                  real, dimension :: Faxa_snowc (:)
                  real, dimension :: Faxa_snowl (:)
                  real, dimension :: Faxa_swndr (:)
                  real, dimension :: Faxa_swvdr (:)
                  real, dimension :: Faxa_swndf (:)
                  real, dimension :: Faxa_swvdf (:)






    end subroutine wrf_to_lilac


end module module_sf_lilac









