#ifndef _GODUNOV_F_H_
#define _GODUNOV_F_H_ $Id: GODUNOV_F.H,v 1.1 1997-07-08 23:08:12 vince Exp $

#ifdef BL_LANG_FORT
#    define FORT_TEST_U_RHO         test_u_rho
#    define FORT_TEST_UMAC_RHO      test_umac_rho
#    define FORT_ESTDT              estdt
#    define FORT_TRANSVEL           transvel
#    define FORT_ESTATE             estate
#    define FORT_ADV_FORCING        adv_forcing
#    define FORT_SYNC_ADV_FORCING   sync_adv_forcing
#    define FORT_SCALINVSID         scalinvsid
#    define FORT_SCALMINMAX         scalminmax
#    define FORT_UPDATE_TF          update_tf
#    define FORT_CORRECT_TF         correct_tf
#    define FORT_UPDATE_AOFS_TF     update_aofs_tf
#    define FORT_UPDATE_AOFS_TF_GP  update_aofs_tf_gp
#    define FORT_SUM_TF_GP_VISC     sum_tf_gp_visc
#    define FORT_SUM_TF_DIVU_VISC   sum_tf_divu_visc
#    define FORT_SUM_TF_GP          sum_tf_gp
#    define FORT_SUM_TF_DIVU        sum_tf_divu
#    define FORT_SET_PARAMS         set_params
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_TEST_U_RHO         TEST_U_RHO 
#    define FORT_TEST_UMAC_RHO      TEST_UMAC_RHO
#    define FORT_ESTDT              ESTDT 
#    define FORT_TRANSVEL           TRANSVEL  
#    define FORT_ESTATE             ESTATE    
#    define FORT_ADV_FORCING        ADV_FORCING    
#    define FORT_SYNC_ADV_FORCING   SYNC_ADV_FORCING   
#    define FORT_SCALINVSID         SCALINVSID
#    define FORT_SCALMINMAX         SCALMINMAX
#    define FORT_UPDATE_TF          UPDATE_TF 
#    define FORT_CORRECT_TF         CORRECT_TF 
#    define FORT_UPDATE_AOFS_TF     UPDATE_AOFS_TF 
#    define FORT_UPDATE_AOFS_TF_GP  UPDATE_AOFS_TF_GP 
#    define FORT_SUM_TF_GP_VISC     SUM_TF_GP_VISC 
#    define FORT_SUM_TF_DIVU_VISC   SUM_TF_DIVU_VISC 
#    define FORT_SUM_TF_GP          SUM_TF_GP 
#    define FORT_SUM_TF_DIVU        SUM_TF_DIVU 
#    define FORT_SET_PARAMS         SET_PARAMS
#  else
#    define FORT_TEST_U_RHO         test_u_rho_
#    define FORT_TEST_UMAC_RHO      test_umac_rho_
#    define FORT_ESTDT              estdt_
#    define FORT_TRANSVEL           transvel_
#    define FORT_ESTATE             estate_
#    define FORT_ADV_FORCING        adv_forcing_
#    define FORT_SYNC_ADV_FORCING   sync_adv_forcing_
#    define FORT_SCALINVSID         scalinvsid_
#    define FORT_SCALMINMAX         scalminmax_
#    define FORT_UPDATE_TF          update_tf_
#    define FORT_CORRECT_TF         correct_tf_
#    define FORT_UPDATE_AOFS_TF     update_aofs_tf_
#    define FORT_UPDATE_AOFS_TF_GP  update_aofs_tf_gp_
#    define FORT_SUM_TF_GP_VISC     sum_tf_gp_visc_
#    define FORT_SUM_TF_DIVU_VISC   sum_tf_divu_visc_
#    define FORT_SUM_TF_GP          sum_tf_gp_
#    define FORT_SUM_TF_DIVU        sum_tf_divu_
#    define FORT_SET_PARAMS         set_params_
#  endif

#include <ArrayLim.H>

  extern "C" {

      // fortran argument checking disabled for now
      void FORT_TEST_U_RHO( ... );   
      void FORT_TEST_UMAC_RHO( ... );   
      void FORT_ESTDT( ... );   
      void FORT_TRANSVEL( ... );   
      void FORT_ESTATE( ... );       
      void FORT_ADV_FORCING( ... );          
      void FORT_SYNC_ADV_FORCING( ... );      
      void FORT_SCALINVSID( ... );   
      void FORT_SCALMINMAX( ... );   
      void FORT_UPDATE_TF( ... );   
      void FORT_CORRECT_TF( ... );   
      void FORT_UPDATE_AOFS_TF( ... );   
      void FORT_UPDATE_AOFS_TF_GP( ... );   
      void FORT_SUM_TF_GP_VISC( ... );   
      void FORT_SUM_TF_DIVU_VISC( ... );   
      void FORT_SUM_TF_GP( ... );   
      void FORT_SUM_TF_DIVU( ... );  

      void FORT_SET_PARAMS(const int& slope_order); 
      
  };
#endif

#endif





