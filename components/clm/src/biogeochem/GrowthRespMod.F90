module GrowthRespMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for growth respiration fluxes,
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use pftvarcon        , only : grperc, grpnow, npcropmin
  use VegetationPropertiesType   , only : veg_vp, vegetation_properties_type
  use CNCarbonFluxType , only : carbonflux_type
  use VegetationType        , only : veg_pp, vegetation_physical_properties
  use VegetationDataType    , only : veg_cf, vegetation_carbon_flux

  !

  use shr_log_mod   , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use subgridAveMod   , only : p2c_1d_filter
  use ColumnDataType  , only : column_carbon_flux
  use ColumnType      , only : col_pp


  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: GrowthResp
  public :: GrowthResp_test
  public :: print_cf
  !-----------------------------------------------------------------------

contains

!        subroutine acc_initialization(mygpu,ngpus)
!                use openacc
!                use spmdMod,  only : iam
!
!                integer, intent(inout) :: mygpu
!                integer, intent(inout) :: ngpus
!
!                ngpus = acc_get_num_devices(acc_device_nvidia)
!                call acc_set_device_num(mod(iam,ngpus),acc_device_nvidia)
!
!                mygpu = acc_get_device_num(acc_device_nvidia)
!                print *, "iam, mygpu: ", iam, mygpu
!        end subroutine




        subroutine summary_rr_test(veg_cf, bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, col_cf)

            !
            ! !DESCRIPTION:
            ! summarize root respiration
           ! !$acc routine seq
            ! !USES:
            !
            ! !ARGUMENTS:
            type(vegetation_carbon_flux) :: veg_cf
            type(bounds_type), intent(in) :: bounds
            integer, intent(in) :: num_soilp
            integer, intent(in) :: filter_soilp(:)
            integer, intent(in) :: num_soilc
            integer, intent(inout) :: filter_soilc(:)

            type(column_carbon_flux) :: col_cf

            !
            ! !LOCAL VARIABLES
            integer :: fp, p, u_p, l_p, u_c, l_c
            !------------------------------------------------------------

            u_p = bounds%endp
            l_p = bounds%begp
            u_c = bounds%endc
            l_c = bounds%begc

            do fp = 1,num_soilp
              p = filter_soilp(fp)
              ! root respiration (RR)
              veg_cf%rr(p) = &
              veg_cf%froot_mr(p) + &
              veg_cf%cpool_froot_gr(p) + &
              veg_cf%cpool_livecroot_gr(p) + &
              veg_cf%cpool_deadcroot_gr(p) + &
              veg_cf%transfer_froot_gr(p) + &
              veg_cf%transfer_livecroot_gr(p) + &
              veg_cf%transfer_deadcroot_gr(p) + &
              veg_cf%cpool_froot_storage_gr(p) + &
              veg_cf%cpool_livecroot_storage_gr(p) + &
              veg_cf%cpool_deadcroot_storage_gr(p)
            enddo
              !call p2c_test(bounds, num_soilc, filter_soilc, &
              !     veg_cf%rr(bounds%begp:bounds%endp), &
              !     col_cf%rr(bounds%begc:bounds%endc))
              call p2c_test(bounds,num_soilc, filter_soilc,veg_cf%rr,col_cf%rr)

           !!!! !$acc update self(filter_soilc)

        end subroutine summary_rr_test


        subroutine p2c_test(bounds, numfc, filterc, pftarr,colarr)
                !(bounds, numfc, filterc,  pftarr, colarr)

            !
            ! !DESCRIPTION:
            ! perform pft to column averaging for single level pft arrays
            !
           ! !$acc routine seq
            ! !ARGUMENTS:
            type(bounds_type)  :: bounds
            integer  :: numfc
            integer, intent(inout) :: filterc(:)
            real(r8) :: pftarr( : )
            real(r8) :: colarr( : )
            !type(column_physical_properties) , target :: col_pp

            ! !LOCAL VARIABLES:
            integer :: fc, c ,p  ! indices
            !-----------------------------------------------------------------------

            ! Enforce expected array sizes
            !SHR_ASSERT_ALL((ubound(pftarr) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
            !SHR_ASSERT_ALL((ubound(colarr) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

            do fc = 1,numfc
               c = filterc(fc)
               colarr(c) = 0._r8
               do p = col_pp%pfti(c), col_pp%pftf(c)
                  if (veg_pp%active(p)) colarr(c) = colarr(c) + pftarr(p) * veg_pp%wtcol(p)
               end do
            end do


       end subroutine p2c_test
  !-----------------------------------------------------------------------
  subroutine GrowthResp_test(num_soilp, filter_soilp, veg_cf, veg_pp,veg_vp,&
          grperc,grpnow,npcropmin)
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables
    !
    ! !USES:
    use spmdMod, only : iam
    use openacc

    !use VegetationPropertiesType   , only : vegetation_properties_type

    !use VegetationType        , only : vegetation_physical_properties
    !use VegetationDataType    , only : vegetation_carbon_flux
    !
    ! !ARGUMENTS:
    integer, intent(in) :: num_soilp       ! number of soil patches in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil patches
    type(vegetation_carbon_flux)           , target :: veg_cf
    type(vegetation_physical_properties)   , target :: veg_pp
    type(vegetation_properties_type) :: veg_vp
    real(r8), allocatable :: grperc(:)      !growth respiration parameter
    real(r8), allocatable :: grpnow(:)
    integer :: npcropmin

    !
    ! !LOCAL VARIABLES:
    integer :: p                ! indices
    integer :: fp               ! lake filter pft index
    integer :: itype_p          ! temp index for veg_pp%itype(p)
    integer :: ngpus

    !-----------------------------------------------------------------------

      ! Loop through patches
      ! start pft loop -- fp loop seems independent

 !$acc data present(veg_cf,veg_pp,veg_vp) copyin(grperc,grpnow,npcropmin,filter_soilp, iam)
 !$acc parallel loop independent private(itype_p, p)
      do fp = 1,num_soilp

        p = filter_soilp(fp)
        itype_p = veg_pp%itype(p)
         if (itype_p >= npcropmin) then ! skip 2 generic crops

            veg_cf%cpool_livestem_gr(p)          = veg_cf%cpool_to_livestemc(p) * grperc(itype_p)

            veg_cf%cpool_livestem_storage_gr(p)  = veg_cf%cpool_to_livestemc_storage(p) * &
                 grperc(itype_p) * grpnow(itype_p)

            veg_cf%transfer_livestem_gr(p)       = veg_cf%livestemc_xfer_to_livestemc(p) * &
                 grperc(itype_p) * (1._r8 - grpnow(itype_p))

            veg_cf%cpool_grain_gr(p)             = veg_cf%cpool_to_grainc(p) * grperc(itype_p)

            veg_cf%cpool_grain_storage_gr(p)     = veg_cf%cpool_to_grainc_storage(p) * &
                 grperc(itype_p) * grpnow(itype_p)

            veg_cf%transfer_grain_gr(p)          = veg_cf%grainc_xfer_to_grainc(p) * grperc(itype_p) &
                 * (1._r8 - grpnow(itype_p))
         end if

         ! leaf and fine root growth respiration
         veg_cf%cpool_leaf_gr(p)          = veg_cf%cpool_to_leafc(p) * grperc(itype_p) + iam+3
         veg_cf%cpool_leaf_storage_gr(p)  = veg_cf%cpool_to_leafc_storage(p) * grperc(itype_p) * &
              grpnow(itype_p) + iam+3
         veg_cf%transfer_leaf_gr(p)       = veg_cf%leafc_xfer_to_leafc(p) * grperc(itype_p) * &
              (1._r8 - grpnow(itype_p))
         veg_cf%cpool_froot_gr(p)         = veg_cf%cpool_to_frootc(p) * grperc(itype_p) + iam+3
         veg_cf%cpool_froot_storage_gr(p) = veg_cf%cpool_to_frootc_storage(p) * grperc(itype_p) * &
              grpnow(itype_p) + iam+3
         veg_cf%transfer_froot_gr(p)      = veg_cf%frootc_xfer_to_frootc(p) * grperc(itype_p) * &
              (1._r8 - grpnow(itype_p))

         if (veg_vp%woody(itype_p) == 1._r8) then
            veg_cf%cpool_livestem_gr(p)          = veg_cf%cpool_to_livestemc(p) * grperc(itype_p) + iam+3
            veg_cf%cpool_livestem_storage_gr(p)  = veg_cf%cpool_to_livestemc_storage(p) * &
                 grperc(itype_p) * grpnow(itype_p) + iam+3
            veg_cf%transfer_livestem_gr(p)       = veg_cf%livestemc_xfer_to_livestemc(p) * &
                 grperc(itype_p) * (1._r8 - grpnow(itype_p))
            veg_cf%cpool_deadstem_gr(p)          = veg_cf%cpool_to_deadstemc(p) * grperc(itype_p) + iam+3
            veg_cf%cpool_deadstem_storage_gr(p)  = veg_cf%cpool_to_deadstemc_storage(p) * &
                 grperc(itype_p) * grpnow(itype_p) + iam+3
            veg_cf%transfer_deadstem_gr(p)       = veg_cf%deadstemc_xfer_to_deadstemc(p) * &
                 grperc(itype_p) * (1._r8 - grpnow(itype_p))
            veg_cf%cpool_livecroot_gr(p)         = veg_cf%cpool_to_livecrootc(p) * grperc(itype_p) + iam+3
            veg_cf%cpool_livecroot_storage_gr(p) = veg_cf%cpool_to_livecrootc_storage(p) * &
                 grperc(itype_p) * grpnow(itype_p) + iam+3
            veg_cf%transfer_livecroot_gr(p)      = veg_cf%livecrootc_xfer_to_livecrootc(p) * &
                 grperc(itype_p) * (1._r8 - grpnow(itype_p))
            veg_cf%cpool_deadcroot_gr(p)         = veg_cf%cpool_to_deadcrootc(p) * grperc(itype_p)
            veg_cf%cpool_deadcroot_storage_gr(p) = veg_cf%cpool_to_deadcrootc_storage(p) * &
                 grperc(itype_p) * grpnow(itype_p)
            veg_cf%transfer_deadcroot_gr(p)      = veg_cf%deadcrootc_xfer_to_deadcrootc(p) * &
                 grperc(itype_p) * (1._r8 - grpnow(itype_p))
         end if

      end do
 !$acc end data

 end subroutine GrowthResp_test

  subroutine print_cf(num_soilp, filter_soilp, veg_cf, veg_vp, veg_pp, npcropmin, grpnow, grperc)

    use spmdMod, only : iam

    integer, intent(in) :: num_soilp       ! number of soil patches in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil patches
    type(vegetation_carbon_flux)           , target :: veg_cf
    type(vegetation_physical_properties)   , target :: veg_pp
    type(vegetation_properties_type) :: veg_vp
    real(r8), allocatable :: grperc(:)      !growth respiration parameter
    real(r8), allocatable :: grpnow(:)
    integer :: npcropmin

    !
    ! !LOCAL VARIABLES:
    integer :: p                ! indices
    integer :: fp               ! lake filter pft index
    integer :: itype_p          ! temp index for veg_pp%itype(p)



        do fp = 1, num_soilp

                p = filter_soilp(fp)
                itype_p = veg_pp%itype(p)
             !   print *, "p:  ",p,"npcropmin : ",npcropmin
             !   print *, "p:  ",p,"itype_p", veg_pp%itype(p)
             !   print *, "p:   ",p,"woody: ",veg_vp%woody(itype_p)
             !   print *, "p:   ",p,"grperc: ", grperc(itype_p)
             !   print *, "p:   ",p,"grpnow:  ", grpnow(itype_p)

             !   if (itype_p >= npcropmin) then
             !           print *, "p:   ",p,"cpool_livestem_gr:",veg_cf%cpool_livestem_gr(p)
             !           print *, "p:   ",p, &
             !                   "cpool_livestem_storage_gr:",veg_cf%cpool_livestem_storage_gr(p)
             !           print *,"p:",p,"transfer_livestem_gr:",&
             !                   veg_cf%transfer_livestem_gr(p)
             !           print *,"p: ",p,"cpool_to_livestemc",&
             !                 veg_cf%cpool_to_livestemc(p)
             !
             !           print *,"p:   ",p,"cpool_grain_gr:",veg_cf%cpool_grain_gr(p)
             !           print *,"p:   ",p,"cpool_to_grainc:",veg_cf%cpool_to_grainc(p)
             !
             !           print *,"p:   ",p,"cpool_grain_storage_gr:",&
             !                   veg_cf%cpool_grain_storage_gr(p)
             !           print *,"p:   ",p,"cpool_to_grainc_storage:"&
             !                   ,veg_cf%cpool_to_grainc_storage(p)

             !
             !            print *,"p:   ",p,"transfer_grain_gr:",&
             !                   veg_cf%transfer_grain_gr(p)
             !           print *,"p:   ",p,"grainc_xfer_to_grainc:",&
             !                   veg_cf%grainc_xfer_to_grainc(p)


             !   end if

                print *,iam,"p:   ",p,"cpool_leaf_gr:",veg_cf%cpool_leaf_gr(p)
                print *,iam,"p:   ",p,"cpool_leaf_storage:",veg_cf%cpool_leaf_storage_gr(p)
             !   print *,iam,"p:   ",p,"cpool_to_livestemc_storage:",&
             !           veg_cf%cpool_to_livestemc_storage(p)
             !   print *,iam,"p:   ",p,"livestemc_xfer_to_livestemc:",&
             !           veg_cf%livestemc_xfer_to_livestemc(p)
             !   print *,iam,"p:   ",p,"cpool_to_leafc:",veg_cf%cpool_to_leafc(p)
             !   print *,iam,"p:   ",p,"cpool_to_leafc_storage:",&
             !           veg_cf%cpool_to_leafc_storage(p)
             !   print *,iam,"p:   ",p,"transfer_leaf_gr(p):",&
             !           veg_cf%transfer_leaf_gr(p)
             !   print *,iam,"p:   ",p,"leafc_xfer_to_leafc:",&
             !           veg_cf%leafc_xfer_to_leafc(p)
                print *,iam,"p:   ",p,"cpool_froot_gr:",veg_cf%cpool_froot_gr(p)
                print *,iam,"p:   ",p,"cpool_froot_storage_gr",&
                        veg_cf%cpool_froot_storage_gr(p)
             !   print *,iam,"p:   ",p,"cpool_to_frootc:",veg_cf%cpool_to_frootc(p)
             !   print *,iam,"p:   ",p,"transfer_froot_gr:",veg_cf%transfer_froot_gr(p)
             !   print *,iam,"p:   ",p,"frootc_xfer_to_frootc:",&
             !           veg_cf%frootc_xfer_to_frootc(p)
                print *,iam,"p:   ",p,"cpool_deadstem_gr:",veg_cf%cpool_deadstem_gr(p)
             !   print *,iam,"p:   ",p,"cpool_to_deadstemc:",&
             !           veg_cf%cpool_to_deadstemc(p)
                print *,iam,"p:   ",p,"cpool_deadstem_storage_gr:",&
                        veg_cf%cpool_deadstem_storage_gr(p)
             !   print *,iam,"p:   ",p,"cpool_to_deadstemc_storage:",&
             !           veg_cf%cpool_to_deadstemc_storage(p)
             !   print *,iam,"p:   ",p,"transfer_deadstem_gr:",veg_cf%transfer_deadstem_gr(p)
             !   print *,iam,"p:   ",p,"deadstemc_xfer_to_deadstemc:",&
             !           veg_cf%deadstemc_xfer_to_deadstemc(p)
                print *,iam,"p:   ",p,"cpool_livecroot_gr:",veg_cf%cpool_livecroot_gr(p)
             !   print *,iam,"p:   ",p,"cpool_to_livecrootc:",veg_cf%cpool_to_livecrootc(p)
             !   print *,iam,"p:   ",p,"cpool_livecroot_storage_gr:",&
             !           veg_cf%cpool_livecroot_storage_gr(p)
             !   print *,iam,"p:   ",p,"cpool_to_livecrootc_storage:",&
             !           veg_cf%cpool_to_livecrootc_storage(p)
             !   print *,iam, "p:   ",p,"transfer_livecroot_gr:",&
             !           veg_cf%transfer_livecroot_gr(p)
             !   print *,iam,"p:   ",p,"livecrootc_xfer_to_livecrootc:",&
             !           veg_cf%livecrootc_xfer_to_livecrootc(p)
             !   print *,iam,"p:   ",p,"cpool_deadcroot_gr:",veg_cf%cpool_deadcroot_gr(p)
             !   print *,iam,"p:   ",p,"cpool_to_deadcrootc:",veg_cf%cpool_to_deadcrootc(p)
             !   print *,iam,"p:   ",p,"cpool_deadcroot_storage_gr:",&
             !           veg_cf%cpool_deadcroot_storage_gr(p)
             !   print *,iam,"p:   ",p,"cpool_to_deadcrootc_storage:",&
             !           veg_cf%cpool_to_deadcrootc_storage(p)
             !   print *,iam,"p:   ",p,"transfer_deadcroot_gr:",&
             !           veg_cf%transfer_deadcroot_gr(p)
             !   print *,iam,"p:   ",p,"deadcrootc_xfer_to_deadcrootc:",&
             !           veg_cf%deadcrootc_xfer_to_deadcrootc(p)


        end do


  end subroutine

  subroutine GrowthResp(num_soilp, filter_soilp, carbonflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables
    !
    ! !USES:

    !use pftvarcon        , only : grperc, grpnow, npcropmin
    !use VegetationPropertiesType   , only : veg_vp
    !use CNCarbonFluxType , only : carbonflux_type
    !use VegetationType        , only : veg_pp
    !use VegetationDataType    , only : veg_cf
    !
    ! !ARGUMENTS:
    integer, intent(in) :: num_soilp       ! number of soil patches in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil patches
    type(carbonflux_type), intent(inout) :: carbonflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p                ! indices
    integer :: fp               ! lake filter pft index
    !-----------------------------------------------------------------------

    associate(                                                                      &
         ivt                           =>    veg_pp%itype                         , & ! Input:  [integer (:)]  pft vegetation type

         woody                         =>    veg_vp%woody                         , & ! Input:  [real(r8) (:)]  binary flag for woody lifeform (1=woody, 0=not woody)

         cpool_to_leafc                =>    veg_cf%cpool_to_leafc                , & ! Input:  [real(r8) (:)]
         cpool_to_leafc_storage        =>    veg_cf%cpool_to_leafc_storage        , & ! Input:  [real(r8) (:)]
         cpool_to_frootc               =>    veg_cf%cpool_to_frootc               , & ! Input:  [real(r8) (:)]
         cpool_to_frootc_storage       =>    veg_cf%cpool_to_frootc_storage       , & ! Input:  [real(r8) (:)]
         cpool_to_livestemc            =>    veg_cf%cpool_to_livestemc            , & ! Input:  [real(r8) (:)]
         cpool_to_livestemc_storage    =>    veg_cf%cpool_to_livestemc_storage    , & ! Input:  [real(r8) (:)]
         cpool_to_deadstemc            =>    veg_cf%cpool_to_deadstemc            , & ! Input:  [real(r8) (:)]
         cpool_to_deadstemc_storage    =>    veg_cf%cpool_to_deadstemc_storage    , & ! Input:  [real(r8) (:)]
         cpool_to_livecrootc           =>    veg_cf%cpool_to_livecrootc           , & ! Input:  [real(r8) (:)]
         cpool_to_livecrootc_storage   =>    veg_cf%cpool_to_livecrootc_storage   , & ! Input:  [real(r8) (:)]
         cpool_to_deadcrootc           =>    veg_cf%cpool_to_deadcrootc           , & ! Input:  [real(r8) (:)]  allocation to dead coarse root C (gC/m2/s)
         cpool_to_deadcrootc_storage   =>    veg_cf%cpool_to_deadcrootc_storage   , & ! Input:  [real(r8) (:)]  allocation to dead coarse root C storage (gC/m2/s)
         cpool_to_grainc               =>    veg_cf%cpool_to_grainc               , & ! Input:  [real(r8) (:)]  allocation to grain C (gC/m2/s)
         cpool_to_grainc_storage       =>    veg_cf%cpool_to_grainc_storage       , & ! Input:  [real(r8) (:)]  allocation to grain C storage (gC/m2/s)
         grainc_xfer_to_grainc         =>    veg_cf%grainc_xfer_to_grainc         , & ! Input:  [real(r8) (:)]  grain C growth from storage (gC/m2/s)
         leafc_xfer_to_leafc           =>    veg_cf%leafc_xfer_to_leafc           , & ! Input:  [real(r8) (:)]  leaf C growth from storage (gC/m2/s)
         frootc_xfer_to_frootc         =>    veg_cf%frootc_xfer_to_frootc         , & ! Input:  [real(r8) (:)]  fine root C growth from storage (gC/m2/s)
         livestemc_xfer_to_livestemc   =>    veg_cf%livestemc_xfer_to_livestemc   , & ! Input:  [real(r8) (:)]  live stem C growth from storage (gC/m2/s)
         deadstemc_xfer_to_deadstemc   =>    veg_cf%deadstemc_xfer_to_deadstemc   , & ! Input:  [real(r8) (:)]  dead stem C growth from storage (gC/m2/s)
         livecrootc_xfer_to_livecrootc =>    veg_cf%livecrootc_xfer_to_livecrootc , & ! Input:  [real(r8) (:)]  live coarse root C growth from storage (gC/m2/s)
         deadcrootc_xfer_to_deadcrootc =>    veg_cf%deadcrootc_xfer_to_deadcrootc , & ! Input:  [real(r8) (:)]  dead coarse root C growth from storage (gC/m2/s)
         cpool_grain_gr                =>    veg_cf%cpool_grain_gr                , & ! InOut:  [real(r8) (:)]
         cpool_grain_storage_gr        =>    veg_cf%cpool_grain_storage_gr        , & ! InOut:  [real(r8) (:)]
         transfer_grain_gr             =>    veg_cf%transfer_grain_gr             , & ! InOut:  [real(r8) (:)]
         cpool_leaf_gr                 =>    veg_cf%cpool_leaf_gr                 , & ! InOut:  [real(r8) (:)]
         cpool_leaf_storage_gr         =>    veg_cf%cpool_leaf_storage_gr         , & ! InOut:  [real(r8) (:)]
         transfer_leaf_gr              =>    veg_cf%transfer_leaf_gr              , & ! InOut:  [real(r8) (:)]
         cpool_froot_gr                =>    veg_cf%cpool_froot_gr                , & ! InOut:  [real(r8) (:)]
         cpool_froot_storage_gr        =>    veg_cf%cpool_froot_storage_gr        , & ! InOut:  [real(r8) (:)]
         transfer_froot_gr             =>    veg_cf%transfer_froot_gr             , & ! InOut:  [real(r8) (:)]
         cpool_livestem_gr             =>    veg_cf%cpool_livestem_gr             , & ! InOut:  [real(r8) (:)]
         cpool_livestem_storage_gr     =>    veg_cf%cpool_livestem_storage_gr     , & ! InOut:  [real(r8) (:)]
         transfer_livestem_gr          =>    veg_cf%transfer_livestem_gr          , & ! InOut:  [real(r8) (:)]
         cpool_deadstem_gr             =>    veg_cf%cpool_deadstem_gr             , & ! InOut:  [real(r8) (:)]
         cpool_deadstem_storage_gr     =>    veg_cf%cpool_deadstem_storage_gr     , & ! InOut:  [real(r8) (:)]
         transfer_deadstem_gr          =>    veg_cf%transfer_deadstem_gr          , & ! InOut:  [real(r8) (:)]
         cpool_livecroot_gr            =>    veg_cf%cpool_livecroot_gr            , & ! InOut:  [real(r8) (:)]
         cpool_livecroot_storage_gr    =>    veg_cf%cpool_livecroot_storage_gr    , & ! InOut:  [real(r8) (:)]
         transfer_livecroot_gr         =>    veg_cf%transfer_livecroot_gr         , & ! InOut:  [real(r8) (:)]
         cpool_deadcroot_gr            =>    veg_cf%cpool_deadcroot_gr            , & ! InOut:  [real(r8) (:)]
         cpool_deadcroot_storage_gr    =>    veg_cf%cpool_deadcroot_storage_gr    , & ! InOut:  [real(r8) (:)]
         transfer_deadcroot_gr         =>    veg_cf%transfer_deadcroot_gr           & ! InOut:  [real(r8) (:)]
         )

!!!acc data copy(veg_pp<grresp>, veg_vp<grresp>, veg_cf)

      ! Loop through patches
      ! start pft loop -- fp loop seems independent
      do fp = 1,num_soilp

        p = filter_soilp(fp)
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops

            cpool_livestem_gr(p)          = cpool_to_livestemc(p) * grperc(ivt(p))

            cpool_livestem_storage_gr(p)  = cpool_to_livestemc_storage(p) * &
                 grperc(ivt(p)) * grpnow(ivt(p))

            transfer_livestem_gr(p)       = livestemc_xfer_to_livestemc(p) * &
                 grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))

            cpool_grain_gr(p)             = cpool_to_grainc(p) * grperc(ivt(p))

            cpool_grain_storage_gr(p)     = cpool_to_grainc_storage(p) * &
                 grperc(ivt(p)) * grpnow(ivt(p))

            transfer_grain_gr(p)          = grainc_xfer_to_grainc(p) * grperc(ivt(p)) &
                 * (1._r8 - grpnow(ivt(p)))
         end if

         ! leaf and fine root growth respiration
         cpool_leaf_gr(p)          = cpool_to_leafc(p) * grperc(ivt(p))
         cpool_leaf_storage_gr(p)  = cpool_to_leafc_storage(p) * grperc(ivt(p)) * &
              grpnow(ivt(p))
         transfer_leaf_gr(p)       = leafc_xfer_to_leafc(p) * grperc(ivt(p)) * &
              (1._r8 - grpnow(ivt(p)))
         cpool_froot_gr(p)         = cpool_to_frootc(p) * grperc(ivt(p))
         cpool_froot_storage_gr(p) = cpool_to_frootc_storage(p) * grperc(ivt(p)) * &
              grpnow(ivt(p))
         transfer_froot_gr(p)      = frootc_xfer_to_frootc(p) * grperc(ivt(p)) * &
              (1._r8 - grpnow(ivt(p)))

         if (woody(ivt(p)) == 1._r8) then
            cpool_livestem_gr(p)          = cpool_to_livestemc(p) * grperc(ivt(p))
            cpool_livestem_storage_gr(p)  = cpool_to_livestemc_storage(p) * &
                 grperc(ivt(p)) * grpnow(ivt(p))
            transfer_livestem_gr(p)       = livestemc_xfer_to_livestemc(p) * &
                 grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))
            cpool_deadstem_gr(p)          = cpool_to_deadstemc(p) * grperc(ivt(p))
            cpool_deadstem_storage_gr(p)  = cpool_to_deadstemc_storage(p) * &
                 grperc(ivt(p)) * grpnow(ivt(p))
            transfer_deadstem_gr(p)       = deadstemc_xfer_to_deadstemc(p) * &
                 grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))
            cpool_livecroot_gr(p)         = cpool_to_livecrootc(p) * grperc(ivt(p))
            cpool_livecroot_storage_gr(p) = cpool_to_livecrootc_storage(p) * &
                 grperc(ivt(p)) * grpnow(ivt(p))
            transfer_livecroot_gr(p)      = livecrootc_xfer_to_livecrootc(p) * &
                 grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))
            cpool_deadcroot_gr(p)         = cpool_to_deadcrootc(p) * grperc(ivt(p))
            cpool_deadcroot_storage_gr(p) = cpool_to_deadcrootc_storage(p) * &
                 grperc(ivt(p)) * grpnow(ivt(p))
            transfer_deadcroot_gr(p)      = deadcrootc_xfer_to_deadcrootc(p) * &
                 grperc(ivt(p)) * (1._r8 - grpnow(ivt(p)))
         end if

      end do

!!! !acc end data

    end associate

  end subroutine GrowthResp
end module GrowthRespMod
