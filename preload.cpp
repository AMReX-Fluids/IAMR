/* */

//#include "bsp.h"
#include "ParallelDescriptor.H"
#include <stdlib.h>
#include <string.h>

#ifdef BL_USE_BSP
extern int BSP_DO_STAT;
extern int BSP_DO_CGPROF;
extern int BSP_DO_PROF;
extern int BSP_NBUFFERS;
extern int BSP_BUFFER_SIZE;
extern int BSP_BUFFER_STALLS;
extern int BSP_THROTTLE_PROCS;
extern int BSP_COMM_FIFO_SIZE;
extern int BSP_OPT_CONTENTION_LEVEL;
extern int BSP_OPT_FCOMBINE_PUTS;
extern int BSP_OPT_FCOMBINE_PUTS_MAX;
extern int BSP_OPT_FCOMBINE_PUTS_MIN;
extern int BSP_OPT_BSMP_BUFFER_SIZE;
extern char *BSP_COMPILE_FLAGS;
extern char *BSP_ARCH;
extern char *BSP_INCLUDE_DIR;
extern int BSP_CHECK_SYNCS;
extern char *BSP_EXEC_FILE;
extern char BSP_LIBRARY_TYPE;
extern int  BSP_OPT_FLIBRARY_LEVEL;

extern "C" {
  void _bsp_preload_init();
  };

void _bsp_preload_init() {
   BSP_DO_CGPROF        = 0;
   BSP_DO_PROF          = 0;
   BSP_DO_STAT          = 0;
   BSP_NBUFFERS         = 2;
   BSP_BUFFER_SIZE      = 10240;
   BSP_THROTTLE_PROCS   = 0;
   BSP_COMM_FIFO_SIZE   = 100;
   BSP_BUFFER_STALLS    = 2;
   BSP_OPT_CONTENTION_LEVEL = 1;
   BSP_OPT_FCOMBINE_PUTS= 0;
   BSP_OPT_FCOMBINE_PUTS_MAX=102400;
   BSP_OPT_FCOMBINE_PUTS_MIN=5120;
   BSP_OPT_BSMP_BUFFER_SIZE =-1;
   BSP_CHECK_SYNCS  =1;
   BSP_LIBRARY_TYPE ='O';
   BSP_OPT_FLIBRARY_LEVEL=2;
 
   BSP_COMPILE_FLAGS  = (char*) malloc(1+strlen(" -flibrary-level 0 -fcontention-resolve 1"));
   BSP_ARCH=(char*) malloc(1+strlen("OSF1"));
   BSP_INCLUDE_DIR=(char*) malloc(1+strlen("/usr/people/vince/Parallel/BSP/BSP1.1/include/"));
   BSP_EXEC_FILE= (char*)malloc(1+strlen("main"));
   if (BSP_COMPILE_FLAGS==NULL || BSP_ARCH==NULL || 
       BSP_INCLUDE_DIR==NULL || BSP_EXEC_FILE==NULL)
     //bsp_abort("{bsp_start}: unable to malloc for compile flags");
     ParallelDescriptor::Abort("{bsp_start}: unable to malloc for compile flags");

   BSP_COMPILE_FLAGS=strcpy(BSP_COMPILE_FLAGS, " -flibrary-level 0 -fcontention-resolve 1");
   BSP_ARCH         =strcpy(BSP_ARCH,"OSF1");
   BSP_INCLUDE_DIR  =strcpy(BSP_INCLUDE_DIR,"/usr/people/vince/Parallel/BSP/BSP1.1/include/");
   BSP_EXEC_FILE    =strcpy(BSP_EXEC_FILE,"main");
}
#endif
