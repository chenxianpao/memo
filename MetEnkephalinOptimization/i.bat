
set ECEPPAK_RES=D:\eceppak\Data\Residue
set ECEPPAK_HYD=D:\eceppak\Data\Residue
set ECEPPAK_BIN=D:\eceppak\bin
set ECEPPAK_IBEM=D:\eceppak\Data\IBEM

set RUNTYP=%1
set INPF=%2
set OUTF=%3
set ANAF=%4
set BNDF=%5

set SOLVDBS=%ECEPPAK_HYD%
set RSDBS=%ECEPPAK_RES%
set SOLVFILE=%ECEPPAK_HYD%\srfopt.set 
set FXVOLFILE=%ECEPPAK_HYD%\fix_vol_tst.parm

set RSDATA=%ECEPPAK_RES%\rsdata
set HLABELFILE=%ECEPPAK_RES%\nmr_protons.data 
set DIHDEFILE=%ECEPPAK_RES%\dihedral_def

set IBEMCHG=%ECEPPAK_IBEM%\resibem

set TGTARS=TARGET.tgtars
set RESTFILE=edmc.rstart
set RANDFILE=edmc.rndom


%ECEPPAK_BIN%\eceppak-NT.exe 


