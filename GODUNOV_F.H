#ifndef _GODUNOV_F_H_
#define _GODUNOV_F_H_ $Id: GODUNOV_F.H,v 1.3 1997-09-29 17:06:24 lijewski Exp $

#ifdef BL_LANG_FORT
#    define FORT_TEST_U_RHO         testurho
#    define FORT_TEST_UMAC_RHO      testumacrho
#    define FORT_ESTDT              estdt
#    define FORT_TRANSVEL           transvel
#    define FORT_ESTATE             estate
#    define FORT_ADV_FORCING        advforcing
#    define FORT_SYNC_ADV_FORCING   syncadvforcing
#    define FORT_SCALINVSID         scalinvsid
#    define FORT_SCALMINMAX         scalminmax
#    define FORT_UPDATE_TF          updatetf
#    define FORT_CORRECT_TF         correcttf
#    define FORT_UPDATE_AOFS_TF     updateaofstf
#    define FORT_UPDATE_AOFS_TF_GP  updateaofstfgp
#    define FORT_SUM_TF_GP_VISC     sumtfgpvisc
#    define FORT_SUM_TF_DIVU_VISC   sumtfdivuvisc
#    define FORT_SUM_TF_GP          sumtfgp
#    define FORT_SUM_TF_DIVU        sumtfdivu
#    define FORT_SET_PARAMS         setparams
#endif

#ifdef BL_LANG_CC
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_TEST_U_RHO         TESTURHO 
#    define FORT_TEST_UMAC_RHO      TESTUMACRHO
#    define FORT_ESTDT              ESTDT 
#    define FORT_TRANSVEL           TRANSVEL  
#    define FORT_ESTATE             ESTATE    
#    define FORT_ADV_FORCING        ADVFORCING    
#    define FORT_SYNC_ADV_FORCING   SYNCADVFORCING   
#    define FORT_SCALINVSID         SCALINVSID
#    define FORT_SCALMINMAX         SCALMINMAX
#    define FORT_UPDATE_TF          UPDATETF 
#    define FORT_CORRECT_TF         CORRECTTF 
#    define FORT_UPDATE_AOFS_TF     UPDATEAOFSTF 
#    define FORT_UPDATE_AOFS_TF_GP  UPDATEAOFSTFGP 
#    define FORT_SUM_TF_GP_VISC     SUMTFGPVISC 
#    define FORT_SUM_TF_DIVU_VISC   SUMTFDIVUVISC 
#    define FORT_SUM_TF_GP          SUMTFGP 
#    define FORT_SUM_TF_DIVU        SUMTFDIVU 
#    define FORT_SET_PARAMS         SETPARAMS
#  else
#    define FORT_TEST_U_RHO         testurho_
#    define FORT_TEST_UMAC_RHO      testumacrho_
#    define FORT_ESTDT              estdt_
#    define FORT_TRANSVEL           transvel_
#    define FORT_ESTATE             estate_
#    define FORT_ADV_FORCING        advforcing_
#    define FORT_SYNC_ADV_FORCING   syncadvforcing_
#    define FORT_SCALINVSID         scalinvsid_
#    define FORT_SCALMINMAX         scalminmax_
#    define FORT_UPDATE_TF          updatetf_
#    define FORT_CORRECT_TF         correcttf_
#    define FORT_UPDATE_AOFS_TF     updateaofstf_
#    define FORT_UPDATE_AOFS_TF_GP  updateaofstfgp_
#    define FORT_SUM_TF_GP_VISC     sumtfgpvisc_
#    define FORT_SUM_TF_DIVU_VISC   sumtfdivuvisc_
#    define FORT_SUM_TF_GP          sumtfgp_
#    define FORT_SUM_TF_DIVU        sumtfdivu_
#    define FORT_SET_PARAMS         setparams_
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
      
  }
#endif

#endif





