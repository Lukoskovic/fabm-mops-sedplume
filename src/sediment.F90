#include "fabm_driver.h"

module mops_sediment

   use fabm_types
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_sediment
      type (type_dependency_id) :: id_bgc_z, id_bgc_dz, id_sedi
      type (type_state_variable_id) :: id_sed
      type (type_diagnostic_variable_id) :: id_finput, id_ACksed

      type (type_bottom_dependency_id) :: id_bgc_z_bot
      type (type_bottom_state_variable_id) :: id_sed_burial
      type (type_bottom_diagnostic_variable_id) :: id_fsed_burial

      real(rk) :: sedlambda, sedwb, sedmartin
      real(rk) :: sedburdige_fac, sedburdige_exp
      real(rk) :: input, sedarea
      real(rk) :: ACksed
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: get_vertical_movement
      procedure :: do_bottom
      procedure :: do
   end type type_mops_sediment

   type (type_universal_standard_variable), parameter :: sediment_total = type_universal_standard_variable(name='total_sediment', conserved=.true., aggregate_variable=.true.)

contains

   subroutine initialize(self, configunit)
      class (type_mops_sediment), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      ! sediment parameters
      call self%get_parameter(self%sedwb, 'sedwb', 'm/d','sediment sinking velocity', default=140.0_rk)
      ! following sediment parameters purely taken from detritus module, which could be equivalently treated to follow a Martin curev
      ! like detritus if assumed the existence of organic matter that can be remineralized
      ! here as default: no remineralization
      call self%get_parameter(self%sedlambda, 'sedlambda', '1/d','sediment remineralization rate', default=0.0_rk) ! LS inspired by the detritus module, a remineralisation is not assumed
      call self%get_parameter(self%sedmartin, 'sedmartin', '-','exponent for Martin curve', default=0.8580_rk)      ! LS purely
      call self%get_parameter(self%sedburdige_fac, 'sedburdige_fac', '-','factor for sediment burial (see Kriest and Oschlies, 2013)', default=1.6828_rk)
      call self%get_parameter(self%sedburdige_exp, 'sedburdige_exp', '-','exponent for sediment burial (see Kriest and Oschlies, 2013)', default=0.799_rk)
      
      call self%get_parameter(self%input, 'input', 'kg/s','sediment source mass rate', default=1.0_rk) ! sediment input parameter in [kg/s]
      call self%get_parameter(self%sedarea, 'sedarea', 'm2','sediment model area', default=1.0_rk) ! assumed sediment discharge area, relevant for the 1D modelling to achieve sediment concentration per volume [kg/m3]
      call self%get_parameter(self%ACksed, 'ACksed', '1/(m*mg/l)', 'attenuation due sediment', default=0.0452_rk)

      ! register state variable for pelagic + benthic
      call self%register_state_variable(self%id_sed, 'c', 'g/l', 'concentration', minimum=0.0_rk)
      call self%register_bottom_state_variable(self%id_sed_burial, 'sedburial', 'kg/m2', 'sediment burial', minimum=0.0_rk)

      ! register diagnostic variables
      call self%register_diagnostic_variable(self%id_finput, 'fsed', 'kg/m2/d', 'sediment input rate', prefill_value=0.0_rk)
      call self%register_diagnostic_variable(self%id_fsed_burial, 'sedburial_flux', 'kg/m2/d', 'sediment burial flux')
      call self%register_diagnostic_variable(self%id_ACksed, 'att_sed', '1/m', 'attenuation due to sediment')


      ! register environmental dependencies
      call self%register_dependency(self%id_bgc_z, standard_variables%depth)
      call self%register_dependency(self%id_bgc_z_bot, standard_variables%bottom_depth)
      call self%register_dependency(self%id_bgc_dz, standard_variables%cell_thickness)
      ! register sediment depth input maske
      call self%register_dependency(self%id_sedi, 'sedinput', '-', 'mask of depth for sediment input source')

      ! aggregatable variables
      call self%add_to_aggregate_variable(sediment_total, self%id_sed)
      call self%add_to_aggregate_variable(sediment_total, self%id_sed_burial)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_sed, scale_factor=self%ACksed*1000.0_rk)

      self%dt = 86400.0_rk
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
    class (type_mops_sediment), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_

    real(rk) ::  idepth, bgc_dz, amount_rate, amount
    real(rk) :: SED, att_sed
    
    _LOOP_BEGIN_
         _GET_(self%id_sed, SED)

        att_sed = self%ACksed*1000.0_rk*SED       ! unit conversion (factor 1000): [kg/m3] = [g/l] -> [mg/l], often related for attenuation coefficient
        _SET_DIAGNOSTIC_(self%id_ACksed, att_sed) 

        _GET_(self%id_bgc_dz, bgc_dz) ! [m]
        _GET_(self%id_sedi, idepth)

        ! transformation of the input from [kg/s] to [g/l/d](assuming 1m wide boxes)
        amount_rate = self%input*idepth*self%dt          ! conversion dicharge mass rate to daily rates -> [kg/s] * [1] * [s/d]  = [kg/d] 

        amount = amount_rate/( self%sedarea*bgc_dz )   ! sediment concentration within assumed discharged area -> [kg/d] * [m-3] = [kg/m3/d]

        _SET_DIAGNOSTIC_(self%id_finput, amount_rate)                ! discharge mass rate as diagnostic [kg/d]
        _ADD_SOURCE_(self%id_sed, amount)                            ! sediemnt concentration [kg/m3/d]
      
    
    _LOOP_END_

   end subroutine do

   subroutine get_vertical_movement(self, _ARGUMENTS_DO_)
      class (type_mops_sediment), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: sedwa, bgc_z, wsed

      sedwa = self%sedlambda/self%sedmartin     ! theoretical martin curve exponent: not of relevance here (assumtion of no remineralization)
      _LOOP_BEGIN_
         _GET_(self%id_bgc_z, bgc_z)
         wsed = self%sedwb + bgc_z*sedwa        ! sediment sinking velocity, here only linear term without any further depth dependencies (no remineralisation)
         _ADD_VERTICAL_VELOCITY_(self%id_sed, -wsed)
      _LOOP_END_
   end subroutine get_vertical_movement

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_mops_sediment), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: sedwa, bgc_z, SED, wsed, fSED, sedflux_l

      ! LS 'benthic' module as implemented in same fashion like detirus
      sedwa = self%sedlambda/self%sedmartin
      _BOTTOM_LOOP_BEGIN_
         _GET_BOTTOM_(self%id_bgc_z_bot, bgc_z)                            ! [m]
         _GET_(self%id_sed, SED)                                           ! [g/l] or [kg/m3]
         wsed = self%sedwb + bgc_z*sedwa                                   ! sediment sinking velocity [m/d]
         fSED = wsed*SED                                                   ! [m/d] * [kg/m3] = [kg/m2/d]
         sedflux_l = MIN(1.0_rk,self%sedburdige_fac*fSED**self%sedburdige_exp)*fSED    ! [kg/m2/d]
         _ADD_BOTTOM_FLUX_(self%id_sed, -sedflux_l)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_fsed_burial, sedflux_l)
         _ADD_BOTTOM_SOURCE_(self%id_sed_burial, sedflux_l) ! bottom source needs to be fed as a flux, i.e. [kg/m2/d]
      _BOTTOM_LOOP_END_

   end subroutine

end module mops_sediment
