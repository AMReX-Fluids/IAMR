# Microsoft Developer Studio Project File - Name="IAMRAll" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=IAMRAll - Win32 Debug3D
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "IAMRAll.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "IAMRAll.mak" CFG="IAMRAll - Win32 Debug3D"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "IAMRAll - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "IAMRAll - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE "IAMRAll - Win32 Debug3D" (based on "Win32 (x86) Console Application")
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
# ADD F90 /include:"Release/" /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN32" /D "BL_USE_NEW_HFILES" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=xilink.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib dfor.lib /nologo /subsystem:console /machine:I386

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
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /Gm /GX /Zi /Od /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN32" /D "BL_USE_NEW_HFILES" /FR /YX /FD /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=xilink.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 dfor.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug3D"
# PROP BASE Intermediate_Dir "Debug3D"
# PROP BASE Ignore_Export_Lib 0
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug3D"
# PROP Intermediate_Dir "Debug3D"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /include:"Debug3D/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD F90 /browser /include:"Debug3D/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /Gm /GX /Zi /Od /D "_CONSOLE" /D "_MBCS" /D "_DEBUG" /D BL_SPACEDIM=3 /D "WIN32" /D "BL_USE_NEW_HFILES" /D "HG_CROSS_STENCIL" /FR /YX /FD /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=xilink.exe
# ADD BASE LINK32 dfor.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 dfor.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "IAMRAll - Win32 Release"
# Name "IAMRAll - Win32 Debug"
# Name "IAMRAll - Win32 Debug3D"
# Begin Group "C Sources"

# PROP Default_Filter "*.cpp"
# Begin Source File

SOURCE=.\ABecLaplacian.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Amr.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\amr_multi.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\AmrLevel.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\BCRec.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\BndryData.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\BndryRegister.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\boundary.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\cache.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\CGSolver.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Cluster.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\CoordSys.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Derive.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Diffusion.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\ErrorList.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\FabSet.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\fill_patch.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\FluxRegister.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Geometry.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Godunov.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_multi1.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_multi2.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_multi3.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_projector.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\interface.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\InterpBndryData.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Interpolater.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\interpolator.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Laplacian.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\LinOp.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MacBndry.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MacOperator.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MacProj.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\main.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Mask.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MultiGrid.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\NavierStokes.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\NS_setup.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\NSBld.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Projection.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\RealBox.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\restrictor.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\RunStats.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\StateData.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\StateDescriptor.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\SyncRegister.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\TagBox.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\ViscBndry.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\ViscBndry2D.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\WriteMultiFab.cpp

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD CPP /I "." /I "..\pboxlib_2" /I ".\include\2d.v9" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE CPP /I "..\pboxlib_2" /I "." /I ".\include\2d.v9" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D BL_SPACEDIM=2 /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"
# ADD CPP /I "..\pboxlib_2" /I "." /I ".\include\3d.v7" /D "BL_LANG_CC" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_FORT_USE_UPPERCASE" /D for="if(0);else for"

!ENDIF 

# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "*.H"
# Begin Source File

SOURCE=.\ABec_F.H
# End Source File
# Begin Source File

SOURCE=.\ABecLaplacian.H
# End Source File
# Begin Source File

SOURCE=.\Amr.H
# End Source File
# Begin Source File

SOURCE=.\Amr_auxil.H
# End Source File
# Begin Source File

SOURCE=.\amr_defs.H
# End Source File
# Begin Source File

SOURCE=.\amr_multi.H
# End Source File
# Begin Source File

SOURCE=.\AmrLevel.H
# End Source File
# Begin Source File

SOURCE=.\ArrayView.H
# End Source File
# Begin Source File

SOURCE=.\BC_TYPES.H
# End Source File
# Begin Source File

SOURCE=.\BCRec.H
# End Source File
# Begin Source File

SOURCE=.\BndryData.H
# End Source File
# Begin Source File

SOURCE=.\BndryRegister.H
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

SOURCE=.\Cluster.H
# End Source File
# Begin Source File

SOURCE=.\CoordSys.H
# End Source File
# Begin Source File

SOURCE=.\COORDSYS_F.H
# End Source File
# Begin Source File

SOURCE=.\DatasetClient.H
# End Source File
# Begin Source File

SOURCE=.\Derive.H
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

SOURCE=.\DIMS.H
# End Source File
# Begin Source File

SOURCE=.\ErrorList.H
# End Source File
# Begin Source File

SOURCE=.\FabSet.H
# End Source File
# Begin Source File

SOURCE=.\fill_patch.H
# End Source File
# Begin Source File

SOURCE=.\FLUXREG_F.H
# End Source File
# Begin Source File

SOURCE=.\FluxRegister.H
# End Source File
# Begin Source File

SOURCE=.\Geometry.H
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

SOURCE=.\INTERP_F.H
# End Source File
# Begin Source File

SOURCE=.\InterpBndryData.H
# End Source File
# Begin Source File

SOURCE=.\INTERPBNDRYDATA_F.H
# End Source File
# Begin Source File

SOURCE=.\Interpolater.H
# End Source File
# Begin Source File

SOURCE=.\interpolator.H
# End Source File
# Begin Source File

SOURCE=.\Laplacian.H
# End Source File
# Begin Source File

SOURCE=.\LevelBld.H
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

SOURCE=.\PROB_AMR_F.H
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

SOURCE=.\RealBox.H
# End Source File
# Begin Source File

SOURCE=.\RegType.H
# End Source File
# Begin Source File

SOURCE=.\restrictor.H
# End Source File
# Begin Source File

SOURCE=.\RunStats.H
# End Source File
# Begin Source File

SOURCE=.\StateData.H
# End Source File
# Begin Source File

SOURCE=.\StateDescriptor.H
# End Source File
# Begin Source File

SOURCE=.\SYNCREG_F.H
# End Source File
# Begin Source File

SOURCE=.\SyncRegister.H
# End Source File
# Begin Source File

SOURCE=.\TagBox.H
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
# Begin Group "Fortran Files"

# PROP Default_Filter ""
# Begin Group "Fortran2D"

# PROP Default_Filter "*.for"
# Begin Source File

SOURCE=.\ABec_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\amr_real2d.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\CG_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\COORDSYS_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\DERIVE_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\DIFFUSION_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\FILCC_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\FLUXREG_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\GODUNOV_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_avg2d.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_multi2d.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_proj2d.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\INTERP_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\INTERPBNDRYDATA_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\LO_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\LP_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MACOPERATOR_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MACPROJ_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MG_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\NAVIERSTOKES_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\PROB_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\PROJECTION_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\SYNCREG_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\VISCOPERATOR_2D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP Exclude_From_Build 1
# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# End Group
# Begin Group "Fortran3D"

# PROP Default_Filter "*.for"
# Begin Source File

SOURCE=.\ABec_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\amr_real3d.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\CG_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\COORDSYS_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\DERIVE_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\DIFFUSION_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\FILCC_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\FLUXREG_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\GODUNOV_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_avg3d.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_multi3d.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\hg_proj3d.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\INTERP_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\INTERPBNDRYDATA_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\LO_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\LP_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MACOPERATOR_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MACPROJ_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\MG_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\NAVIERSTOKES_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\PROB_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\PROJECTION_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\SYNCREG_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\VISCOPERATOR_3D.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# PROP Exclude_From_Build 1
# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# PROP Exclude_From_Build 1
# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# PROP BASE Exclude_From_Build 1
# ADD BASE F90 /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# End Group
# Begin Source File

SOURCE=.\GODUNOV_F.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\LO_UTIL.for

!IF  "$(CFG)" == "IAMRAll - Win32 Release"

# ADD F90 /iface:cref

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug"

# ADD F90 /browser /iface:cref /dbglibs

!ELSEIF  "$(CFG)" == "IAMRAll - Win32 Debug3D"

# ADD BASE F90 /iface:cref /dbglibs
# ADD F90 /iface:cref /dbglibs

!ENDIF 

# End Source File
# End Group
# End Target
# End Project
