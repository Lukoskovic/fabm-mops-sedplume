#include "fabm_driver.h"

module mops_zooplankton

   use fabm_types
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_zooplankton
      type (type_state_variable_id) :: id_c, id_phy, id_po4, id_din, id_oxy, id_det, id_dop, id_dic, id_sed, id_alk
      type (type_diagnostic_variable_id) :: id_f2
      ! VS: introducing id_det_prod_zoo (see below)
      type (type_diagnostic_variable_id) :: id_det_prod_zoo
      ! VS: diagnostic to complete all carbon_c flux diagnostics
      type (type_diagnostic_variable_id) :: id_zooexu

      real(rk) :: ACmuzoo, ACkphy, AClambdaz, AComniz, ACeff, graztodop, zlambda
      real(rk) :: ACmuzoofac, ACefffac, AComnizfac, zlambdafac  ! PLUME
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_mops_zooplankton), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      call self%get_parameter(self%ro2ut, 'ro2ut', 'mol O2/mol P','redfield -O2:P ratio', default=151.13958_rk)
      call self%get_parameter(self%ACmuzoo, 'ACmuzoo', '1/d','max. grazing rate', default=1.893_rk)
      call self%get_parameter(self%ACkphy, 'ACkphy', 'mmol P/m3','half-saturation constant', default=SQRT(self%ACmuzoo/1.0_rk)/rnp)
      call self%get_parameter(self%ACeff, 'ACeff', '1','assimilation efficiency', default=0.75_rk)
      call self%get_parameter(self%graztodop, 'graztodop', '1','fraction of grazing that goes into DOP', default=0.0_rk)
      call self%get_parameter(self%AClambdaz, 'AClambdaz', '1/d','excretion', default=0.03_rk)
      call self%get_parameter(self%AComniz, 'AComniz', 'm3/(mmol P * day)','density dependent loss rate', default=4.548_rk)
      call self%get_parameter(self%zlambda, 'zlambda', '1/d','mortality', default=0.01_rk)

! VS remove minimum value parameter to avoid clipping in TMM implementation acc. to Jorns, October 16, 2024
      call self%register_state_variable(self%id_c, 'c', 'mmol P/m3', 'concentration', minimum=0.0_rk)
      call self%register_diagnostic_variable(self%id_f2, 'f2', 'mmol P/m3/d', 'grazing')
      ! VS: diagnostic to complete all carbon_c flux diagnostics
      call self%register_diagnostic_variable(self%id_zooexu, 'zooexu', 'mmol P/m3/d', 'exudation')
      ! VS introducing detritus production by zooplankton as diagnostic variable
      call self%register_diagnostic_variable(self%id_det_prod_zoo, 'det_prod_zoo', 'mmol P/m3/d', 'detritus produced by zooplankton')

      call self%register_state_dependency(self%id_phy, 'phy', 'mmol P/m3', 'phytoplankton')
      call self%register_state_dependency(self%id_dop, 'dop', 'mmol P/m3', 'dissolved organic phosphorus')
      call self%register_state_dependency(self%id_det, 'det', 'mmol P/m3', 'detritus')
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol O2/m3', 'oxygen')
      call self%register_state_dependency(self%id_din, 'din', 'mmol N/m3', 'dissolved inorganic nitrogen')
      call self%register_state_dependency(self%id_po4, 'pho', 'mmol P/m3', 'phosphate')
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol C/m3', 'dissolved inorganic carbon')
      call self%register_state_dependency(self%id_alk, 'alk', 'mmol/m3', 'alkalinity')

      ! LS PLUME stuff (adding some reducing factors for parameters in presence of sediment)
      call self%get_parameter(self%ACmuzoofac, 'ACmuzoofac', '1', 'reduction factor for max. grazing rate', default=1.0_rk)
      call self%get_parameter(self%ACefffac, 'ACefffac', '1', 'reduction factor for assimilation efficiency', default=1.0_rk)
      call self%get_parameter(self%AComnizfac, 'AComnizfac', '1', 'reduction factor for density dependent loss rate', default=1.0_rk)
      call self%get_parameter(self%zlambdafac, 'zlambdafac', '1', 'reduction factor mortality', default=1.0_rk)
      call self%register_state_dependency(self%id_sed, 'sed', 'g/l', 'sediment')

      ! Register environmental dependencies
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c)
      ! VS also consider total carbon and total nitrogen
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c, scale_factor=rcp)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=rnp)
      ! VS zooplankton (like phytoplankton) detritus production contributes to total detritus production by plankton
      call self%add_to_aggregate_variable(detritus_production_by_plankton, self%id_det_prod_zoo)
      ! VS an aggregate variable for all biogeochemical DIC
      call self%add_to_aggregate_variable(total_dic, self%id_dic)

      self%dt = 86400.0_rk
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_mops_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: PHY, ZOO
      real(rk) :: graz0, graz, zooexu, zooloss
      real(rk) :: SED ! PLUME
      real(rk) :: muzoofac, efffac, omnizfac, zlambdaf    ! PLUME

      _LOOP_BEGIN_

      _GET_(self%id_phy, PHY)
      _GET_(self%id_c, ZOO)

      _GET_(self%id_sed, SED) ! PLUME

       if(ZOO.gt.0.0_rk) then

         if(SED.gt.1.0e-4_rk) then ! SED > 10e-4

            muzoofac    = self%ACmuzoofac
            efffac      = self%ACefffac
            omnizfac    = self%AComnizfac
            zlambdaf    = self%zlambdafac

         else ! SED <= 10e-4

            muzoofac       = 1
            efffac         = 1
            omnizfac       = 1
            zlambdaf       = 1

         endif

         if(PHY.gt.0.0_rk) then

! Grazing of zooplankton, Holling III
           graz0=muzoofac*self%ACmuzoo*PHY*PHY/(self%ACkphy*self%ACkphy+PHY*PHY)*ZOO ! LS PLUME
! Make sure not to graze more phytoplankton than available.
           graz = MIN(PHY,graz0*bgc_dt)/bgc_dt

         else !PHY < 0

           graz=0.0_rk

         endif !ZOO

! Zooplankton exudation
          zooexu = self%AClambdaz * ZOO

! Zooplankton mortality 
          zooloss = omnizfac * self%AComniz * ZOO * ZOO ! LS PLUME

       else !ZOO < 0

           graz   =0.0_rk
           zooexu = 0.0_rk
           zooloss = 0.0_rk

       endif !ZOO

       _SET_DIAGNOSTIC_(self%id_f2, graz)
       _SET_DIAGNOSTIC_(self%id_zooexu, zooexu)
! VS detritus production by zooplankton
       _SET_DIAGNOSTIC_(self%id_det_prod_zoo, (1.0_rk-self%graztodop)*(1.0_rk-self%ACeff*efffac)*graz + (1.0_rk-self%graztodop)*zooloss)

! Collect all euphotic zone fluxes in these arrays.
        _ADD_SOURCE_(self%id_c, self%ACeff*efffac*graz-zooexu-zooloss) ! LS PLUME
        _ADD_SOURCE_(self%id_po4, zooexu)
        _ADD_SOURCE_(self%id_dop, self%graztodop*(1.0_rk-self%ACeff*efffac)*graz + self%graztodop*zooloss) ! LS PLUME
        _ADD_SOURCE_(self%id_oxy, -zooexu*ro2ut)
        _ADD_SOURCE_(self%id_phy, -graz)
        _ADD_SOURCE_(self%id_det, (1.0_rk-self%graztodop)*(1.0_rk-self%ACeff*efffac)*graz + (1.0_rk-self%graztodop)*zooloss) ! LS PLUME

        _ADD_SOURCE_(self%id_din, zooexu*rnp)
        _ADD_SOURCE_(self%id_dic, zooexu*rcp)
        _ADD_SOURCE_(self%id_alk, -zooexu*(rnp+1))

         ZOO = MAX(ZOO - alimit*alimit, 0.0_rk)
         _ADD_SOURCE_(self%id_c, -self%zlambda*zlambdaf*ZOO)   ! LS PLUME 
         _ADD_SOURCE_(self%id_dop, self%zlambda*zlambdaf*ZOO)  ! LS PLUME 

      _LOOP_END_
   end subroutine do

end module mops_zooplankton
