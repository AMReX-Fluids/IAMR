# Microsoft Developer Studio Project File - Name="IAMRAll" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=IAMRAll - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "IAMRAll.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "IAMRAll.mak" CFG="IAMRAll - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "IAMRAll - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "IAMRAll - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=xicl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /include:"Release/" /compile_only /nologo /warn:nofileopt
# ADD F90 /include:"Release/" /compile_only /nologo /iface:cref /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "." /I "..\amrlib" /I "..\bndrylib" /I "..\pBoxLib_2" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_PARALLEL_IO" /D for="if(0);else for" /D "HG_CROSS_STENCIL" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=xilink.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /iface:cref /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /Gm /GX /Zi /Od /I "." /I "..\amrlib" /I "..\bndrylib" /I "..\pBoxLib_2" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D BL_SPACEDIM=3 /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_PARALLEL_IO" /D for="if(0);else for" /D "HG_CROSS_STENCIL" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=xilink.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "IAMRAll - Win32 Release"
# Name "IAMRAll - Win32 Debug"
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=.\ABecLaplacian.cpp
# End Source File
# Begin Source File

SOURCE=.\amr_multi.cpp
# End Source File
# Begin Source File

SOURCE=.\BndryData.cpp
# End Source File
# Begin Source File

SOURCE=.\boundary.cpp
# End Source File
# Begin Source File

SOURCE=.\cache.cpp
# End Source File
# Begin Source File

SOURCE=.\CGSolver.cpp
# End Source File
# Begin Source File

SOURCE=.\Diffusion.cpp
# End Source File
# Begin Source File

SOURCE=.\fill_patch.cpp
# End Source File
# Begin Source File

SOURCE=.\Godunov.cpp
# End Source File
# Begin Source File

SOURCE=.\hg_multi1.cpp
# End Source File
# Begin Source File

SOURCE=.\hg_multi2.cpp
# End Source File
# Begin Source File

SOURCE=.\hg_multi3.cpp
# End Source File
# Begin Source File

SOURCE=.\hg_projector.cpp
# End Source File
# Begin Source File

SOURCE=.\interface.cpp
# End Source File
# Begin Source File

SOURCE=.\InterpBndryData.cpp
# End Source File
# Begin Source File

SOURCE=.\interpolator.cpp
# End Source File
# Begin Source File

SOURCE=.\Laplacian.cpp
# End Source File
# Begin Source File

SOURCE=.\LinOp.cpp
# End Source File
# Begin Source File

SOURCE=.\MacBndry.cpp
# End Source File
# Begin Source File

SOURCE=.\MacOperator.cpp
# End Source File
# Begin Source File

SOURCE=.\MacProj.cpp
# End Source File
# Begin Source File

SOURCE=.\main.cpp
# End Source File
# Begin Source File

SOURCE=.\Mask.cpp
# End Source File
# Begin Source File

SOURCE=.\MultiGrid.cpp
# End Source File
# Begin Source File

SOURCE=.\NavierStokes.cpp
# End Source File
# Begin Source File

SOURCE=.\NS_setup.cpp
# End Source File
# Begin Source File

SOURCE=.\NSBld.cpp
# End Source File
# Begin Source File

SOURCE=.\Projection.cpp
# End Source File
# Begin Source File

SOURCE=.\restrictor.cpp
# End Source File
# Begin Source File

SOURCE=.\SyncRegister.cpp
# End Source File
# Begin Source File

SOURCE=.\ViscBndry.cpp
# End Source File
# Begin Source File

SOURCE=.\ViscBndry2D.cpp
# End Source File
# Begin Source File

SOURCE=.\WriteMultiFab.cpp
# End Source File
# End Group
# Begin Group "C++ Headers"

# PROP Default_Filter "h"
# Begin Source File

SOURCE=.\ABec_F.H
# End Source File
# Begin Source File

SOURCE=.\ABecLaplacian.H
# End Source File
# Begin Source File

SOURCE=.\amr_defs.H
# End Source File
# Begin Source File

SOURCE=.\amr_multi.H
# End Source File
# Begin Source File

SOURCE=.\ArrayView.H
# End Source File
# Begin Source File

SOURCE=.\BndryData.H
# End Source File
# Begin Source File

SOURCE=.\boundary.H
# End Source File
# Begin Source File

SOURCE=.\BoundCond.H
# End Source File
# Begin Source File

SOURCE=.\cache.H
# End Source File
# Begin Source File

SOURCE=.\CG_F.H
# End Source File
# Begin Source File

SOURCE=.\CGSolver.H
# End Source File
# Begin Source File

SOURCE=.\DatasetClient.H
# End Source File
# Begin Source File

SOURCE=.\DERIVE_F.H
# End Source File
# Begin Source File

SOURCE=.\Diffusion.H
# End Source File
# Begin Source File

SOURCE=.\DIFFUSION_F.H
# End Source File
# Begin Source File

SOURCE=.\fill_patch.H
# End Source File
# Begin Source File

SOURCE=.\Godunov.H
# End Source File
# Begin Source File

SOURCE=.\GODUNOV_F.H
# End Source File
# Begin Source File

SOURCE=.\hg_multi.H
# End Source File
# Begin Source File

SOURCE=.\hg_projector.H
# End Source File
# Begin Source File

SOURCE=.\interface.H
# End Source File
# Begin Source File

SOURCE=.\InterpBndryData.H
# End Source File
# Begin Source File

SOURCE=.\INTERPBNDRYDATA_F.H
# End Source File
# Begin Source File

SOURCE=.\interpolator.H
# End Source File
# Begin Source File

SOURCE=.\Laplacian.H
# End Source File
# Begin Source File

SOURCE=.\LinOp.H
# End Source File
# Begin Source File

SOURCE=.\LO_BCTYPES.H
# End Source File
# Begin Source File

SOURCE=.\LO_F.H
# End Source File
# Begin Source File

SOURCE=.\LP_F.H
# End Source File
# Begin Source File

SOURCE=.\MacBndry.H
# End Source File
# Begin Source File

SOURCE=.\MacOperator.H
# End Source File
# Begin Source File

SOURCE=.\MACOPERATOR_F.H
# End Source File
# Begin Source File

SOURCE=.\MacOpMacDrivers.H
# End Source File
# Begin Source File

SOURCE=.\MacOpProjDrivers.H
# End Source File
# Begin Source File

SOURCE=.\MacProj.H
# End Source File
# Begin Source File

SOURCE=.\MACPROJ_F.H
# End Source File
# Begin Source File

SOURCE=.\Mask.H
# End Source File
# Begin Source File

SOURCE=.\MG_F.H
# End Source File
# Begin Source File

SOURCE=.\MultiGrid.H
# End Source File
# Begin Source File

SOURCE=.\NavierStokes.H
# End Source File
# Begin Source File

SOURCE=.\NAVIERSTOKES_F.H
# End Source File
# Begin Source File

SOURCE=.\PROB_F.H
# End Source File
# Begin Source File

SOURCE=.\probdata.H
# End Source File
# Begin Source File

SOURCE=.\Projection.H
# End Source File
# Begin Source File

SOURCE=.\PROJECTION_F.H
# End Source File
# Begin Source File

SOURCE=.\RegType.H
# End Source File
# Begin Source File

SOURCE=.\restrictor.H
# End Source File
# Begin Source File

SOURCE=.\SYNCREG_F.H
# End Source File
# Begin Source File

SOURCE=.\SyncRegister.H
# End Source File
# Begin Source File

SOURCE=.\ViscBndry.H
# End Source File
# Begin Source File

SOURCE=.\ViscBndry2D.H
# End Source File
# Begin Source File

SOURCE=.\VISCOPERATOR_F.H
# End Source File
# Begin Source File

SOURCE=.\WriteMultiFab.H
# End Source File
# End Group
# Begin Group "Fortran"

# PROP Default_Filter ""
# Begin Group "3D"

# PROP Default_Filter "F"
# Begin Group "TEMPS"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\ABec_3D.For
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\amr_real3d.For
# PROP Exclude_From_Build 1
# End Source File
# End Group
# Begin Source File

SOURCE=.\ABec_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\amr_real3d.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\CG_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\DERIVE_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\DIFFUSION_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\GODUNOV_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\INTERPBNDRYDATA_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\LO_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\LP_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\MACOPERATOR_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\MACPROJ_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\MG_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\PROJECTION_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\SYNCREG_3D.F
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\VISCOPERATOR_3D.F
# PROP Exclude_From_Build 1
# End Source File
# End Group
# End Group
# End Target
# End Project
