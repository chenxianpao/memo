# Microsoft Developer Studio Project File - Name="MetEnkephalinOptimizationLib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=MetEnkephalinOptimizationLib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "MetEnkephalinOptimizationLib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "MetEnkephalinOptimizationLib.mak" CFG="MetEnkephalinOptimizationLib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "MetEnkephalinOptimizationLib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "MetEnkephalinOptimizationLib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "MetEnkephalinOptimizationLib - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "MetEnkephalinOptimizationLib - Win32 Debug"

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
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ  /c
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /TP /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "MetEnkephalinOptimizationLib - Win32 Release"
# Name "MetEnkephalinOptimizationLib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\ga\GA1DArrayGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\GA1DBinStrGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\GA2DArrayGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\GA2DBinStrGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\GA3DArrayGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\GA3DBinStrGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAAllele.C
# End Source File
# Begin Source File

SOURCE=..\ga\GABaseGA.C
# End Source File
# Begin Source File

SOURCE=..\ga\GABin2DecGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\gabincvt.C
# End Source File
# Begin Source File

SOURCE=..\ga\GABinStr.C
# End Source File
# Begin Source File

SOURCE=..\ga\GADCrowdingGA.C
# End Source File
# Begin Source File

SOURCE=..\ga\GADemeGA.C
# End Source File
# Begin Source File

SOURCE=..\ga\gaerror.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAIncGA.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAList.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAListBASE.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAListGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAParameter.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAPopulation.C
# End Source File
# Begin Source File

SOURCE=..\ga\garandom.C
# End Source File
# Begin Source File

SOURCE=..\ga\GARealGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAScaling.C
# End Source File
# Begin Source File

SOURCE=..\ga\GASelector.C
# End Source File
# Begin Source File

SOURCE=..\ga\GASimpleGA.C
# End Source File
# Begin Source File

SOURCE=..\ga\GASStateGA.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAStatistics.C
# End Source File
# Begin Source File

SOURCE=..\ga\GAStringGenome.C
# End Source File
# Begin Source File

SOURCE=..\ga\GATree.C
# End Source File
# Begin Source File

SOURCE=..\ga\GATreeBASE.C
# End Source File
# Begin Source File

SOURCE=..\ga\GATreeGenome.C
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\ga\ga.h
# End Source File
# Begin Source File

SOURCE=..\ga\GA1DArrayGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\GA1DBinStrGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\GA2DArrayGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\GA2DBinStrGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\GA3DArrayGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\GA3DBinStrGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAAllele.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAArray.h
# End Source File
# Begin Source File

SOURCE=..\ga\GABaseGA.h
# End Source File
# Begin Source File

SOURCE=..\ga\GABin2DecGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\gabincvt.h
# End Source File
# Begin Source File

SOURCE=..\ga\GABinStr.h
# End Source File
# Begin Source File

SOURCE=..\ga\gaconfig.h
# End Source File
# Begin Source File

SOURCE=..\ga\GADCrowdingGA.h
# End Source File
# Begin Source File

SOURCE=..\ga\GADemeGA.h
# End Source File
# Begin Source File

SOURCE=..\ga\gaerror.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAEvalData.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\gaid.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAIncGA.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAList.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAListBASE.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAListGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAMask.h
# End Source File
# Begin Source File

SOURCE=..\ga\GANode.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAParameter.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAPopulation.h
# End Source File
# Begin Source File

SOURCE=..\ga\garandom.h
# End Source File
# Begin Source File

SOURCE=..\ga\GARealGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAScaling.h
# End Source File
# Begin Source File

SOURCE=..\ga\GASelector.h
# End Source File
# Begin Source File

SOURCE=..\ga\GASimpleGA.h
# End Source File
# Begin Source File

SOURCE=..\ga\GASStateGA.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAStatistics.h
# End Source File
# Begin Source File

SOURCE=..\ga\GAStringGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\GATree.h
# End Source File
# Begin Source File

SOURCE=..\ga\GATreeBASE.h
# End Source File
# Begin Source File

SOURCE=..\ga\GATreeGenome.h
# End Source File
# Begin Source File

SOURCE=..\ga\gatypes.h
# End Source File
# Begin Source File

SOURCE=..\ga\gaversion.h
# End Source File
# End Group
# End Target
# End Project
