!-------------------------------------------------------------------------
module module_create_gindex
  implicit none
contains

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

!========================================================================

module module_lilac_reshape
  implicit none
contains

! TODO (NS) : lilac_reshape should be one 2d --> 1d and 3d -->1d ...
!Either add to subroutines for each of them
!Or add an option for each of them to this subroutine

  subroutine lilac_reshape (ide, jde, its, ite, jts, jte, var_2d, var_1d)
    integer, intent(in) :: ide  ! domain end index, i
    integer, intent(in) :: jde  ! domain end index, j
    integer, intent(in) :: its  ! task start index, i
    integer, intent(in) :: ite  ! task end index, i
    integer, intent(in) :: jts  ! task start index, j
    integer, intent(in) :: jte  ! task end index, j
    integer, allocatable, intent(in)  :: var_2d (:,:)
    integer, allocatable, intent(out) :: var_1d(:)

    integer :: ite_limited  ! task end index on the mass point grid, i
    integer :: jte_limited  ! task end index on the mass point grid, j
    integer :: num_points
    integer :: i, j, n



    ! The very last index in both row & column space is just used on the momentum grid.
    ! Here we are just working with the mass point grid, so we need to ignore that last
    ! index.

    num_points = ((ite_limited - its + 1) * (jte_limited - jts + 1))
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




!========================================================================
! WRF Lilac subrouti
!========================================================================

module module_sf_lilac
!-------------------------------------------------------------------------
    use shr_kind_mod, only: r8 => shr_kind_r8
    use mpi           ,  only : MPI_COMM_WORLD, MPI_COMM_NULL, MPI_Init, MPI_FINALIZE, MPI_SUCCESS

    
    ! NOTE : This subroutine is NOT being used now! 
    subroutine lilac_wrf_init ()
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


    ! NOTE : This subroutine is NOT being used now! 
    subroutine lilac_clm_init ( gindex_atm, atm2lnd    &
                   ,zgcmxy     ,forc_qxy   ,ps   ,forc_txy    ,tsxy  &  
                   ,shxy       ,qfx        ,lhxy        ,soiflx      ,qgh       &
                   ,gsw, swdown,ra_sw_physics &   
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
                   ,inest, sf_urban_physics,nlcat)


          type (atm2lnd_type), intent(inout):: atm2lnd


          
 
          ! arguments
          integer , allocatable, intent(inout) :: gindex_atm(:)
          type (atm2lnd_type), intent(inout):: atm2lnd

          integer                              :: ntasks
          integer                              :: mytask
          character(len=128)                   :: filename
          !-----------------------------------------------------------------
          ! MPI initalization is done in WRF initializations

          ! decomposition in WRF
          ! allocate and initalize atm2lnd
          call init_atm2lnd_type(begc, endc,atm2lnd)


          atm2lnd%forc_t     (begc:endc ) = forc_txy
          atm2lnd%forc_u     (begc:endc ) = forc_uxy
          atm2lnd%forc_v     (begc:endc ) = forc_vxy
          atm2lnd%forc_wind  (begc:endc ) = sqrt(forc_uxy**2 + forc_vxy**2)
          atm2lnd%forc_q     (begc:endc ) = forc_qxy
          !atm2lnd%forc_rh    (begc:endc ) = 
          atm2lnd%forc_hgt   (begc:endc ) = zgcmxy
          atm2lnd%forc_hgt_u (begc:endc ) = zgcmxy !observational height of wind [m]
          atm2lnd%forc_hgt_t (begc:endc ) = zgcmxy !observational height of wind [m]
          atm2lnd%forc_hgt_q (begc:endc ) = zgcmxy !observational height of wind [m]

          atm2lnd%forc_pbot  (begc:endc ) = forc_pbotxy
          atm2lnd%forc_psrf  (begc:endc ) = forc_psrfxy
          !atm2lnd%forc_th    (begc:endc ) = 
          !atm2lnd%forc_vp    (begc:endc ) = 
          !atm2lnd%forc_rho   (begc:endc ) = 
          atm2lnd%forc_pco2  (begc:endc ) = pco2 * clm_a2l%forc_pbot(g)
          atm2lnd%forc_po2  (begc:endc )  = po2  * clm_a2l%forc_pbot(g)

          !Fluxes
          atm2lnd%forc_lwrad (begc:endc ) = flwdsxy
          clm_a2l%forc_solad (begc:endc ,1) = forc_solsxy
          clm_a2l%forc_solad (begc:endc ,2) = forc_sollxy
          atm2lnd%forc_solai (begc:endc ,1) = forc_solsdxy
          atm2lnd%forc_solai (begc:endc ,2) = forc_solldxy
          atm2lnd%forc_solar (begc:endc )   =   forc_solsxy + forc_sollxy &
                                                + forc_solsdxy + forc_solldxy


    end subroutine lilac_clm_init


          ! allocate and fill in atm2lnd

    subroutine lilac_ctsm_run (restart_alarm_is_ringing, stop_alarm_is_ringing, &
                  dz8w,qv_curr,p8w, t_phy,tsk,                   &
                  hfx,qfx,lh,grdflx,qgh,gsw,swdown,               &
                  ra_sw_physics,history_interval,glw,smstav,smstot, &
                  sfcrunoff,udrunoff,ivgtyp,isltyp,vegfra,        &
                  albedo,znt,z0, tmn,xland,xice, emiss,           &
                  snowc,qsfc,rainbl,maxpatch,                     &
                  num_soil_layers,dtbl,xtime, dt,dzs,             &
                  smois,tslb,snow,canwat,                         &
                  chs,chs2,sh2o,snowh,                            &
                  u_phy,v_phy,                                    &
                  shdmin,shdmax,                                  &
                  acsnom,acsnow,                                  &
                  dx,xlat,xlong,ht,                               &
                  ids,ide, jds,jde, kds,kde,                      &
                  ims,ime, jms,jme, kms,kme,                      &
                  i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte, &
                  inest,sf_urban_physics, nlcat )

          logical , intent (in) :: restart_alarm_is_ringing  ! Logical flag for if the alarm rings this timestep for restart
          logical , intent (in) :: stop_alarm_is_ringing     ! Logical flag for if the alarm rings this timestep for ending the run 


          
          ! This does the same thing as atm_driver_to_lilac 
          ! In LILAC, this is atm_driver_to_lilac
          call wrf_to_lilac (   &
                  dz8w,qv_curr,p8w, t_phy,tsk,                   &
                  hfx,qfx,lh,grdflx,qgh,gsw,swdown,               &
                  ra_sw_physics,history_interval,glw,smstav,smstot, &
                  sfcrunoff,udrunoff,ivgtyp,isltyp,vegfra,        &
                  albedo,znt,z0, tmn,xland,xice, emiss,           &
                  snowc,qsfc,rainbl,maxpatch,                     &
                  num_soil_layers,dtbl,xtime, dt,dzs,             &
                  smois,tslb,snow,canwat,                         &
                  chs,chs2,sh2o,snowh,                            &
                  u_phy,v_phy,                                    &
                  shdmin,shdmax,                                  &
                  acsnom,acsnow,                                  &
                  dx,xlat,xlong,ht,                               &
                  ids,ide, jds,jde, kds,kde,                      &
                  ims,ime, jms,jme, kms,kme,                      &
                  i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte, &
                  inest,sf_urban_physics, nlcat )


          call lilac_run (restart_alarm_is_ringing=.true.,stop_alarm_is_ringing=.true.)

    end subroutine lilac_ctsm_run


    subroutine wrf_to_lilac ( &
                    zgcmxy   , forc_txy , forc_uxy , forc_qxy, &     !dz8w     , t_phy , u_phy, v_phy, &
                    ps  ,  flwdsxy, &      ! p8w , glw 
                    swvisdir, swvisdif, swnirdir, swnirdif, &  ! & swvisdir, swvisdif, swnirdir, swnirdif
                    prec  , ! !total input precipitation  (mm; accumlated precipitation within DT)
                    !dz8w     ,forc_qxy   ,ps   ,forc_txy    ,tsxy  &
                    !dz8w     ,forc_qxy   ,ps   ,forc_txy    ,tsxy  &
                   ,shxy       ,qfx        ,lhxy        ,soiflx      ,qgh       &
                   ,gsw, swdown,ra_sw_physics &
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
                   ,inest, sf_urban_physics,nlcat               )


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


                 !integer :: im_size
                 !integer :: jm_size

                  im_size = ime - ims
                  jm_size = jme - jms

                  !Sa_topo = reshape (dz8w   (ims:ime,kms,jms:jme), (/im_size* jm_size/))

                  !flwdsxy = reshape (Faxa_lwdn(ims:ime,kms,jms:jme), (/jm_size* jm_size/))

                  ! convert the 2d to 1d array
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Sa_z   , zgcmxy)
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Sa_topo   , )
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Sa_u  , forc_uxy)
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Sa_v  , forc_vxy)
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Sa_ptem  , forc_th)
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Sa_pbot  , forc_pbotxy )
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Sa_tbot  , forc_txy)
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Sa_shum  , forc_qxy )

                  call lilac_reshape (ide, jde, its, ite, jts, jte, Faxa_lwdn , flwdsxy)
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Faxa_rainc , )
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Faxa_rainl , )
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Faxa_snowc , )
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Faxa_rainc , )

                  call lilac_reshape (ide, jde, its, ite, jts, jte, Faxa_swndn , forc_sollxy)
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Faxa_swvdr , forc_solsxy)
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Faxa_swndf , forc_solldxy)
                  call lilac_reshape (ide, jde, its, ite, jts, jte, Faxa_swvdf , forc_solsdxy)


                  call lilac_atmcap_atm2lnd('Sa_z'       , Sa_z)
                  call lilac_atmcap_atm2lnd('Sa_topo'    , Sa_topo)
                  call lilac_atmcap_atm2lnd('Sa_u'       , Sa_u)
                  call lilac_atmcap_atm2lnd('Sa_v'       , Sa_v)
                  call lilac_atmcap_atm2lnd('Sa_ptem'    , Sa_ptem)
                  call lilac_atmcap_atm2lnd('Sa_pbot'    , Sa_pbot)
                  call lilac_atmcap_atm2lnd('Sa_tbot'    , Sa_tbot)
                  call lilac_atmcap_atm2lnd('Sa_shum'    , Sa_shum)
                  call lilac_atmcap_atm2lnd('Faxa_lwdn'  , Faxa_lwdn)
                  call lilac_atmcap_atm2lnd('Faxa_rainc' , Faxa_rainc)
                  call lilac_atmcap_atm2lnd('Faxa_rainl' , Faxa_rainl)
                  call lilac_atmcap_atm2lnd('Faxa_snowc' , Faxa_snowc)
                  call lilac_atmcap_atm2lnd('Faxa_snowl' , Faxa_snowl)
                  call lilac_atmcap_atm2lnd('Faxa_swndr' , Faxa_swndr)
                  call lilac_atmcap_atm2lnd('Faxa_swvdr' , Faxa_swvdr)
                  call lilac_atmcap_atm2lnd('Faxa_swndf' , Faxa_swndf)
                  call lilac_atmcap_atm2lnd('Faxa_swvdf' , Faxa_swvdf)


    end subroutine wrf_to_lilac







    ! For now the below subroutine is not used and we are directly using reshape
    subroutine ptrcast (pointer1d, pointer2d)
        ! If we want to avoid a copy for wrf to lilac...
    ! this convert the 2d pointers to 1d pointer for passing to lilac
        real, pointer :: pointer1d(:)
        real, pointer :: pointer2d(:,:)
        integer       :: n
        integer       :: m


        allocate (pointer2d (n , m) )
        pointer1d(n*m) => pointer2d

    end subroutine


    subrouine lilac_reshape ()
    end subroutine lilac_reshape


    subroutine init_atm2lnd_type(begc, endc, atm2lnd)
        !
        ! !DESCRIPTION:
        ! Initialize atmospheric variables required by the land
        !
        ! !ARGUMENTS:
          implicit none
          integer, intent(in) :: begc, endc
          type (atm2lnd_type), intent(inout):: atm2lnd

        ! !LOCAL VARIABLES:
          real(r8) :: ival   ! initial value
          integer  :: begc_atm, end_atm
        !-----------------------------------------------------------------
        !allocating these values from atmosphere for now!

          allocate(atm2lnd%forc_t     (begc:endc) )
          allocate(atm2lnd%forc_u     (begc:endc) )
          allocate(atm2lnd%forc_v     (begc:endc) )
          allocate(atm2lnd%forc_wind  (begc:endc) )
          allocate(atm2lnd%forc_q     (begc:endc) )
          allocate(atm2lnd%forc_rh    (begc:endc) )
          allocate(atm2lnd%forc_hgt   (begc:endc) )
          allocate(atm2lnd%forc_hgt_u (begc:endc) )
          allocate(atm2lnd%forc_hgt_t (begc:endc) )
          allocate(atm2lnd%forc_hgt_q (begc:endc) )
          allocate(atm2lnd%forc_pbot  (begc:endc) )
          allocate(atm2lnd%forc_th    (begc:endc) )
          allocate(atm2lnd%forc_vp    (begc:endc) )
          allocate(atm2lnd%forc_rho   (begc:endc) )
          allocate(atm2lnd%forc_psrf  (begc:endc) )
          allocate(atm2lnd%forc_pco2  (begc:endc) )
          allocate(atm2lnd%forc_lwrad (begc:endc) )
          allocate(atm2lnd%forc_solad (begc:endc ,numrad))
          allocate(atm2lnd%forc_solai (begc:endc ,numrad))
          allocate(atm2lnd%forc_solar (begc:endc))
          allocate(atm2lnd%forc_rain  (begc:endc))
          allocate(atm2lnd%forc_snow  (begc:endc))
          allocate(atm2lnd%forc_ndep  (begc:endc))
          allocate(atm2lnd%rainf      (begc:endc))


          ! should we initialize to a special value????
          atm2lnd%forc_t     (begc:endc ) = ival
          atm2lnd%forc_u     (begc:endc ) = ival
          atm2lnd%forc_v     (begc:endc ) = ival
          atm2lnd%forc_wind  (begc:endc ) = ival
          atm2lnd%forc_q     (begc:endc ) = ival
          atm2lnd%forc_rh    (begc:endc ) = ival
          atm2lnd%forc_hgt   (begc:endc ) = ival
          atm2lnd%forc_hgt_u (begc:endc ) = ival
          atm2lnd%forc_hgt_t (begc:endc ) = ival
          atm2lnd%forc_hgt_q (begc:endc ) = ival
          atm2lnd%forc_pbot  (begc:endc ) = ival
          atm2lnd%forc_th    (begc:endc ) = ival
          atm2lnd%forc_vp    (begc:endc ) = ival
          atm2lnd%forc_rho   (begc:endc ) = ival
          atm2lnd%forc_psrf  (begc:endc ) = ival
          atm2lnd%forc_pco2  (begc:endc ) = ival
          atm2lnd%forc_lwrad (begc:endc ) = ival
          atm2lnd%forc_solad (begc:endc ,1:numrad) = ival
          atm2lnd%forc_solai (begc:endc ,1:numrad) = ival
          atm2lnd%forc_solar (begc:endc ) = ival
          atm2lnd%forc_rain  (begc:endc ) = ival
          atm2lnd%forc_snow  (begc:endc ) = ival
          atm2lnd%forc_ndep  (begc:endc ) = ival
          atm2lnd%rainf      (begc:endc ) = ival

    end subroutine init_atm2lnd_type


end module module_sf_lilac


!========================================================================

module lilac_type_mod

    type atm2lnd_type

      real(r8), pointer :: forc_t(:)       !atmospheric temperature (Kelvin)
      real(r8), pointer :: forc_u(:)       !atm wind speed, east direction (m/s)
      real(r8), pointer :: forc_v(:)       !atm wind speed, north direction (m/s)
      real(r8), pointer :: forc_wind(:)    !atmospheric wind speed
      real(r8), pointer :: forc_q(:)       !atmospheric specific humidity (kg/kg)
      real(r8), pointer :: forc_hgt(:)     !atmospheric reference height (m)
      real(r8), pointer :: forc_hgt_u(:)   !obs height of wind [m] (new)
      real(r8), pointer :: forc_hgt_t(:)   !obs height of temperature [m] (new)
      real(r8), pointer :: forc_hgt_q(:)   !obs height of humidity [m] (new)
      real(r8), pointer :: forc_pbot(:)    !atmospheric pressure (Pa)
      real(r8), pointer :: forc_th(:)      !atm potential temperature (Kelvin)
      real(r8), pointer :: forc_vp(:)      !atmospheric vapor pressure (Pa)
      real(r8), pointer :: forc_rho(:)     !density (kg/m**3)
      real(r8), pointer :: forc_rh(:)      !atmospheric relative humidity (%)
      real(r8), pointer :: forc_psrf(:)    !surface pressure (Pa)
      real(r8), pointer :: forc_pco2(:)    !CO2 partial pressure (Pa)
      real(r8), pointer :: forc_lwrad(:)   !downwrd IR longwave radiation (W/m**2)
      real(r8), pointer :: forc_solad(:,:) !direct beam radiation (numrad)
                                           !(vis=forc_sols , nir=forc_soll )

      real(r8), pointer :: forc_solai(:,:) !diffuse radiation (numrad)
                                           !(vis=forc_solsd, nir=forc_solld)
      real(r8), pointer :: forc_solar(:)   !incident solar radiation
      real(r8), pointer :: forc_rain(:)    !rain rate [mm/s]
      real(r8), pointer :: forc_snow(:)    !snow rate [mm/s]
      real(r8), pointer :: forc_ndep(:)    !nitrogen deposition rate (gN/m2/s)
      real(r8), pointer :: rainf(:)        !ALMA rain+snow [mm/s]

    end type atm2lnd_type

end module lilac_type_mod






