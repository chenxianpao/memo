#if !defined(ZGJ_INFO_ZJUT_ENERGY_CALCULATOR)
#define ZGJ_INFO_ZJUT_ENERGY_CALCULATOR

#include "fstream.h"
#include "iomanip.h"
#include "stdlib.h"

class EnergyCalculator {
public:
	EnergyCalculator();
	virtual ~EnergyCalculator(){}

/******************************************************************************************/
public:
	double	singleEnergyEvaluate(double *pGEOM, int length);
	double	singleEnergyMinimize(double *pGEOM, int length, bool *SPEC);
	double	singleEnergyMinimize(double *pGEOM, int length);

	int		multipleEnergyEvaluate(double** arrayGeom, int arraySize, int angleLength, double *resultValue);
	int		multipleEnergyMinimize(double** arrayGeom, int arraySize, int angleLength, bool *spec, double *resultValue);
	int		multipleEnergyMinimize(double** arrayGeom, int arraySize, int angleLength, double *resultValue);

/******************************************************************************************/
private:
	int		singleEnergyEvaluateWriteSourceFile(double *pGEOM, int length);
	double	singleEnergyEvaluateReadResultFile( );
	void	singleEnergyEvaluateClearFiles( );

	int		singleEnergyMinimizeWriteSourceFile(double *pGEOM, int length, bool *SPEC);
	double	singleEnergyMinimizeReadResultFile(double*pGEOM);
	void	singleEnergyMinimizeClearFiles( );

	int		multipleEnergyEvaluateWriteSourceFile(double **pGEOM, int popSize, int angleLength);
	int		multipleEnergyEvaluateReadOUTOFile(double** arrayGeom, int arraySize, double* resultValue);
	int		multipleEnergyEvaluateReadMainResultFile( );
	void	multipleEnergyEvaluateClearFiles( );

	int		multipleEnergyMinimizeWriteSourceFile(double **arrayGeom, int arraySize, int angleLength, bool* spec);
	int		multipleEnergyMinimizeReadOUTOFile(double** arrayGeom, int arraySize, double* resultValue);
	int		multipleEnergyMinimizeReadMainResultFile( );
	void	multipleEnergyMinimizeClearFiles( );
/******************************************************************************************/
public:
	int		eceppFormatWrite(double *pGeom, int length); 
	int		eceppFormatWrite(ofstream& of, double *pGeom, int length);

public:
	double	getEvalCount() {return energyEvalCount;}
	double	getEvalSec() {return energyEvalSec;}
	bool	isOutputPDBFormat( ) { return outputPDBFormat;}
	bool	setOutputPDBFormat(bool flag ) {return outputPDBFormat=flag;}
	void	setPDBFile(char* str) { strcpy(pdbFile, str);}

private:
	double	energyEvalCount;
	double	energyEvalSec;
	bool	outputPDBFormat;

private:
	double	bigNumber;
	char	singleEnergyEvaluateSourceFile[50];
	char	singleEnergyEvaluateResultFile[50];
	char	singleEnergyMinimizeSourceFile[50];
	char	singleEnergyMinimizeResultFile[50];
	char	singleEnergyMinimizeOUTOFile[50];

	char	multipleEnergyEvaluateSourceFile[50];
	char	multipleEnergyEvaluateInputFile[50];
	char	multipleEnergyEvaluateResultFile[50];
	char	multipleEnergyEvaluateOUTOFile[50];

	char	multipleEnergyMinimizeSourceFile[50];
	char	multipleEnergyMinimizeInputFile[50];
	char	multipleEnergyMinimizeMainResultFile[50];
	char	multipleEnergyMinimizeOUTOResultFile[50];

	char	errFile[10];
	char	logFile[10];
	char	pdbFile[20];
/******************************************************************************************/

};

#endif 
