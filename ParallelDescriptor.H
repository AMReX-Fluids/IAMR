//BL_COPYRIGHT_NOTICE
// -----------------------------------------------------------------
// ParallelDescriptor.H
// -----------------------------------------------------------------
#ifndef _PARALLELDESCRIPTOR_H_
#define _PARALLELDESCRIPTOR_H_

#include <Assert.H>
#include <Boolean.H>
#include <REAL.H>
#include <Box.H>
#include <BaseFab.H>
#include <FArrayBox.H>
#include "bsp.h"
#include "bsp_level1.h"

#define StartParallel(nprocs)    bsp_begin(nprocs);
#define StartParallelAllProcs()  bsp_begin(bsp_nprocs());
#define EndParallel()            bsp_end();

// some of these can be collapsed
struct FabComTag {
  int fromProc;
  int toProc;
  int fabIndex;
  int fineIndex;
  int crseIndex;
  int srcComp;
  int destComp;
  int nComp;
  int face;
  int dataType;
  int fabArrayId;
  int fillBoxId;
  int procThatNeedsData;
  int procThatHasData;
  Box box;
  Box ovlpBox;
  Box subBox;
};


class ParallelDescriptor {

  public:
    static int MyProc() { return bsp_pid();    }
    static int NProcs() { return bsp_nprocs(); }
    static void Synchronize() { bsp_sync(); }
    static void Synchronize(const char *msg) {
	  cout << "----- " << bsp_pid() << " :  about to sync:  " << msg << endl;
					       bsp_sync();
					     }
    static bool IOProcessor() { return (bsp_pid() == ioProcessor); }
    static int  IOProcessorNumber() { return (ioProcessor); }
    static void Abort(char *msg, ...)  { bsp_abort(msg); }
    static double second()  { return bsp_time(); }

    // reduction operations
    template<class T> static void ReduceMin(T &rvar);

    // bool
    static void ReduceBoolAnd(bool &rvar);
    // Real
    static void ReduceRealPlus(Real &rvar);
    static void ReduceRealMax(Real &rvar);
    static void ReduceRealMin(Real &rvar);
    // int
    static void ReduceIntPlus(int &rvar);
    static void ReduceIntMax(int &rvar);
    static void ReduceIntMin(int &rvar);
    // long
    static void ReduceLongPlus(long &rvar);
    static void ReduceLongMax(long &rvar);
    static void ReduceLongMin(long &rvar);


    // data transfer functions
    static void ShareVar(const void *var, int bytes) {
      bsp_pushregister(var, bytes);
    }

    static void UnshareVar(const void *var) {
      bsp_popregister(var);
    }

    static void WriteData(int procnum, const void *src, void *dest,
			  int offset, int bytes)
    {
      bsp_put(procnum, src, dest, offset, bytes);
    }

    static void ReadData(int procnum, const void *src, int offset,
			 void *dest, int bytes)
    {
      bsp_get(procnum, src, offset, dest, bytes);
    }

    static void SetMessageHeaderSize(int &messageHeaderSize) {
      bsp_set_tag_size(&messageHeaderSize);
    } 

    static void GetMessageHeader(int &dataSize, void *messageHeader) {
      bsp_get_tag(&dataSize, messageHeader);
    } 

    static void SendData(int toproc, const FabComTag &messageHeader,
			 const void *data, int datasizeinbytes)
    {
      bsp_send(toproc, &messageHeader, data, datasizeinbytes);
    }

    static void ReceiveData(void *data, int datasizeinbytes) {
      bsp_move(data, datasizeinbytes);
    }

    template<class T> static void Broadcast(int toproc, T *src, T *dest) {
      bsp_bcast(toproc, (void *) src, (void *) dest, sizeof(T));
    }

    //static void SendFab(int procnum, const FArrayBox &fab);

  private:
    static const int ioProcessor = 0;
    static const int syncCount   = 0;

};

#endif /*_PARALLELDESCRIPTOR_H_*/
