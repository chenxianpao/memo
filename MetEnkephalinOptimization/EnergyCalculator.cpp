
#include "stdafx.h"
#include "EnergyCalculator.h"


EnergyCalculator::EnergyCalculator(){
	energyEvalCount = 0.0;
	energyEvalSec = 0.0;
	outputPDBFormat = false;

	bigNumber = 1E25;
	strcpy(singleEnergyEvaluateSourceFile,		"SingleEnergyEvaluateSourceFile.inp");
	strcpy(singleEnergyEvaluateResultFile,		"main_out.MET_SINGLE_EVALUATE");

	strcpy(singleEnergyMinimizeSourceFile,		"SingleEnergyMinimizeSourceFile.inp");
	strcpy(singleEnergyMinimizeResultFile,		"main_out.MET_SINGLE_MINIMIZE");
	strcpy(singleEnergyMinimizeOUTOFile,		"outo.MET_SINGLE_MINIMIZE");

	strcpy(multipleEnergyEvaluateSourceFile,	"MultipleEnergyEvaluateSourceFile.inp");
	strcpy(multipleEnergyEvaluateInputFile,		"outo.MultipleEnergyEvaluateInput");
	strcpy(multipleEnergyEvaluateResultFile,	"main_out.MET_MULTIPLE_EVALUATE");
	strcpy(multipleEnergyEvaluateOUTOFile,		"outo.MET_MULTIPLE_EVALUATE");

	strcpy(multipleEnergyMinimizeSourceFile,	"MultipleEnergyMinimizeSourceFile.inp");
	strcpy(multipleEnergyMinimizeInputFile,		"outo.MultipleEnergyMinimizeInput");
	strcpy(multipleEnergyMinimizeMainResultFile,"main_out.MET_MULTIPLE_MINIMIZE");
	strcpy(multipleEnergyMinimizeOUTOResultFile,"outo.MET_MULTIPLE_MINIMIZE");

	strcpy(errFile,  "err");
	strcpy(logFile,  "log");
	strcpy(pdbFile,  "Met");
}


double EnergyCalculator::singleEnergyEvaluate(double *pGEOM, int length) {
	double ETOT = bigNumber;
	singleEnergyEvaluateWriteSourceFile(pGEOM, length);
	system("i.bat ENERGY SingleEnergyEvaluateSourceFile MET_SINGLE_EVALUATE x x  1>> log 2> err" );
	ETOT = singleEnergyEvaluateReadResultFile();
	singleEnergyEvaluateClearFiles();

	return ETOT;
}

double EnergyCalculator::singleEnergyMinimize(double *pGEOM, int length, bool *SPEC){
	double ETOT = bigNumber;
	singleEnergyMinimizeWriteSourceFile(pGEOM, length, SPEC);
	system("i.bat MINIMIZE SingleEnergyMinimizeSourceFile MET_SINGLE_MINIMIZE x x  1>> log 2> err");
	ETOT = singleEnergyMinimizeReadResultFile(pGEOM);
	singleEnergyMinimizeClearFiles();

	return ETOT;
}

double EnergyCalculator::singleEnergyMinimize(double *pGEOM, int length) {
	double ETOT = bigNumber;
	
	bool* spec = new bool[length];
	for(int i=0; i<24; i++) spec[i] = true;

	ETOT = singleEnergyMinimize(pGEOM, length, spec);

	return ETOT;
}


int EnergyCalculator::multipleEnergyEvaluate(double** pGEOM, int popSize, int angleLength, double *resultValue) {
	multipleEnergyEvaluateWriteSourceFile(pGEOM, popSize, angleLength);
	system("i.bat ENERGY MultipleEnergyEvaluateSourceFile MET_MULTIPLE_EVALUATE MultipleEnergyEvaluateInput rr  1>> log 2> err");
	multipleEnergyEvaluateReadOUTOFile(pGEOM, popSize, resultValue);
	multipleEnergyEvaluateReadMainResultFile();
	multipleEnergyEvaluateClearFiles();

	return 1;
}

int	EnergyCalculator::multipleEnergyMinimize(double** arrayGeom, int arraySize, int angleLength, bool *spec, double *resultValue) {
	multipleEnergyMinimizeWriteSourceFile(arrayGeom, arraySize, angleLength, spec);
	system("i.bat MINIMIZE MultipleEnergyMinimizeSourceFile MET_MULTIPLE_MINIMIZE MultipleEnergyMinimizeInput rr  1>> log 2> err");
	multipleEnergyMinimizeReadOUTOFile(arrayGeom, arraySize, resultValue);
	multipleEnergyMinimizeReadMainResultFile();
	multipleEnergyMinimizeClearFiles();

	return 1;
}

int	EnergyCalculator::multipleEnergyMinimize(double** arrayGeom, int arraySize, int angleLength, double *resultValue) {
	bool *spec = new bool[angleLength];
	for(int i=0; i<angleLength; i++) spec[i] = true;
	multipleEnergyMinimize(arrayGeom, arraySize, angleLength, spec, resultValue);

	return 1;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////
int	EnergyCalculator::eceppFormatWrite(double *pGeom, int length) {
	if(length != 24) return 0;

	for(int k=0; k<length; k++) {
		cout<< setw(8)								// 设定field宽度
			<< setprecision(3)						// 设置小数位置
			<< setiosflags(ios::showpoint)			// keep trailing 0s
			<< setiosflags(ios::fixed)				// 使用这些设置
			<< pGeom[k];
		if (k==5 || k==8 || k==11 || k==16 || k==23) cout << endl;
	}
	return 1;
}

int EnergyCalculator::eceppFormatWrite(ofstream& of, double *pGeom, int length){
	if(length != 24 || !of) return 0;

	for(int k=0; k<length; k++) {
		of	<< setw(8)								// 设定field宽度
			<< setprecision(3)						// 设置小数位置
			<< setiosflags(ios::showpoint)			// keep trailing 0s
			<< setiosflags(ios::fixed)				// 使用这些设置
			<< pGeom[k];
		if (k==5 || k==8 || k==11 || k==16 || k==23) of << endl;
	}
	return 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
int EnergyCalculator::singleEnergyEvaluateWriteSourceFile(double *geomData, int length) {
	ofstream outFile(singleEnergyEvaluateSourceFile, ios::out);
	
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = energy" << endl;
	if(isOutputPDBFormat()) {
		outFile << "PRINT_CART" << endl;
		outFile << "OUTFORMAT = PDB" << endl;
		outFile << "FILE  = " <<pdbFile << endl;
	}
	outFile << "res_code= one_letter" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	outFile << "YGGFM" << endl;
	outFile << "O" << endl;
	outFile << "$END" << endl;

	outFile << "$GEOM" << endl << endl;
	for(int i=0; i<6;i++) {
		outFile	<< setw(8)						// 设定field宽度
				<< setprecision(3)				// 设置小数位置
				<< setiosflags(ios::showpoint)	// keep trailing 0s
				<< setiosflags(ios::fixed)		// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;
	
	for(i=6; i<9;i++) {
		outFile	<< setw(8)						// 设定field宽度
				<< setprecision(3)				// 设置小数位置
				<< setiosflags(ios::showpoint)	// keep trailing 0s
				<< setiosflags(ios::fixed)		// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=9; i<12;i++) {
		outFile	<< setw(8)						// 设定field宽度
				<< setprecision(3)				// 设置小数位置
				<< setiosflags(ios::showpoint)	// keep trailing 0s
				<< setiosflags(ios::fixed)		// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=12; i<17;i++) {
		outFile	<< setw(8)						// 设定field宽度
				<< setprecision(3)				// 设置小数位置
				<< setiosflags(ios::showpoint)	// keep trailing 0s
				<< setiosflags(ios::fixed)		// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=17; i<24;i++) {
		outFile	<< setw(8)						// 设定field宽度
				<< setprecision(3)				// 设置小数位置
				<< setiosflags(ios::showpoint)	// keep trailing 0s
				<< setiosflags(ios::fixed)		// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;
	outFile << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();

	return 1;
}

double EnergyCalculator::singleEnergyEvaluateReadResultFile( ) {
	double ETOT = bigNumber;

	energyEvalCount = 0.0;
	energyEvalSec = 0.0;

	const int LINE_LENGTH  = 1000; 
	char wordBuf[LINE_LENGTH]={0};
	ifstream fin(singleEnergyEvaluateResultFile);
	do{
		fin	>>  wordBuf;
		if(strcmp(wordBuf, "ETOT") == 0) {
			fin >> wordBuf;
			ETOT = atof(wordBuf);
			if(strlen(wordBuf) == 0 ) {
				cout << "ETOT read error!" << endl;
				cout.flush();
				getchar();
			}
			break;
		}
	} while(!fin.eof());

	// read Evals
	if(fin.eof()) {cout << "read format error!" << endl; cout.flush(); getchar();}
	do{
		fin	>> wordBuf;
		if(strcmp(wordBuf, "=")==0) {
			fin >> wordBuf;
			energyEvalCount = atof(wordBuf);
			break;
		}
	} while(!fin.eof());
	
	// read read Evals/sec
	if(fin.eof()) {cout << "read format error!" << endl; cout.flush(); getchar();}
	do{
		fin	>> wordBuf;
		if(strcmp(wordBuf, "=")==0) {
			fin >> wordBuf;
			energyEvalSec = atof(wordBuf);
			break;
		}
	} while(!fin.eof());
	if(fin.eof()) {cout << "read format error!" << endl; cout.flush(); getchar();}

	fin.clear();
	fin.close();

	energyEvalSec = energyEvalCount/energyEvalSec;

	return ETOT;
}

void EnergyCalculator::singleEnergyEvaluateClearFiles() {
	ofstream outFile1(singleEnergyEvaluateSourceFile, ios::out);
	outFile1.flush();
	outFile1.close();

	ofstream outFile2(singleEnergyEvaluateResultFile, ios::out);
	outFile2.flush();
	outFile2.close();

	ofstream outFile3(logFile, ios::out);
	outFile3.flush();
	outFile3.close();
}


int	EnergyCalculator::singleEnergyMinimizeWriteSourceFile(double *pGEOM, int length, bool *spec) {
	ofstream outFile(singleEnergyMinimizeSourceFile, ios::out);

	outFile << "$CNTRL" << endl;
	outFile << "runtyp = minimize" << endl;
	if(isOutputPDBFormat()) {
		outFile << "PRINT_CART" << endl;
		outFile << "OUTFORMAT = PDB" << endl;
		outFile << "FILE  = " <<pdbFile << endl;
	}
	outFile << "res_code= one_letter" << endl;
	outFile << "var_angles = spec" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	outFile << "YGGFM" << endl;
	outFile << "O" << endl;
	outFile << "$END" << endl << endl;


	outFile <<"$SPEC" << endl;
	
	int varNumber = 0;
	int i = 0;

	outFile << "2" << "  ";
	for(i=0; i<6; i++) {if(spec[i]) varNumber++;}
	outFile << varNumber << "  ";
	for(i=0; i<6; i++) {if(spec[i]) outFile << i+1 << "  ";}
	outFile << endl;

	varNumber = 0;
	outFile << "3" << "  ";
	for(i=0; i<3; i++) {if(spec[i+6]) varNumber++;}
	outFile << varNumber << "  ";
	for(i=0; i<3; i++) {if(spec[i+6]) outFile << i+1 << "  ";}
	outFile << endl;

	varNumber = 0;
	outFile << "4" << "  ";
	for(i=0; i<3; i++) {if(spec[i+9]) varNumber++;}
	outFile << varNumber << "  ";
	for(i=0; i<3; i++) {if(spec[i+9]) outFile << i+1 << "  ";}
	outFile << endl;

	varNumber = 0;
	outFile << "5" << "  ";
	for(i=0; i<5; i++) {if(spec[i+12]) varNumber++;}
	outFile << varNumber << "  ";
	for(i=0; i<5; i++) {if(spec[i+12]) outFile << i+1 << "  ";}
	outFile << endl;

	varNumber = 0;
	outFile << "6" << "  ";
	for(i=0; i<7; i++) {if(spec[i+17]) varNumber++;}
	outFile << varNumber << "  ";
	for(i=0; i<7; i++) {if(spec[i+17]) outFile << i+1 << "  ";}
	outFile << endl;

	outFile <<"$END" << endl << endl;

	outFile << "$GEOM" << endl << endl;
	for(i=0; i<6;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< pGEOM[i];
	}
	outFile << endl;

	for(i=6; i<9;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< pGEOM[i];
	}
	outFile << endl;

	for(i=9; i<12;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< pGEOM[i];
	}
	outFile << endl;

	for(i=12; i<17;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< pGEOM[i];
	}
	outFile << endl;

	for(i=17; i<24;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< pGEOM[i];
	}
	outFile << endl;
	outFile << endl;
	outFile << "$END" << endl;


	outFile.flush();
	outFile.close();

	return 1;
}

double EnergyCalculator::singleEnergyMinimizeReadResultFile(double *pGEOM) {
	double ETOT = bigNumber;

	energyEvalCount = 0.0;
	energyEvalSec = 0.0;
	
	const int LINE_LENGTH = 1200; 
	char lineBuf[LINE_LENGTH]={0};
	char wordBuf[200]={0};
	char substr[9]={0};
	int  index = 0;

	// read "outo.MET_SINGLE_MINIMIZE"
	ifstream fin1(singleEnergyMinimizeOUTOFile);

	fin1 >> wordBuf;
	if(strlen(wordBuf)==0) {			// there have no minimizer!
		fin1.close();						
		return ETOT;
	}

	fin1 >> wordBuf;					// read ETOT value;
	ETOT = atof(wordBuf);
	fin1.getline(lineBuf, LINE_LENGTH); // skip first line;
	fin1.getline(lineBuf, LINE_LENGTH); // skip second line (sequence);
	fin1.getline(lineBuf, LINE_LENGTH); // skip NH2 end group line;

	fin1.getline(lineBuf, LINE_LENGTH);	// Tyr
	for(int m=0; m<6;m++) {
		for(int n=0; n<8;n++) {
			substr[n] = lineBuf[m*8+n];
		}
		pGEOM[index++] =(float)atof(substr);
	}

	fin1.getline(lineBuf, LINE_LENGTH);	// Gly1
	for(m=0; m<3;m++) {
		for(int n=0; n<8;n++) {
			substr[n] = lineBuf[m*8+n];
		}
		pGEOM[index++] =(float)atof(substr);
	}

	fin1.getline(lineBuf, LINE_LENGTH);	// Gly2
	for(m=0; m<3;m++) {
		for(int n=0; n<8;n++) {
			substr[n] = lineBuf[m*8+n];
		}
		pGEOM[index++] =(float)atof(substr);
	}

	fin1.getline(lineBuf, LINE_LENGTH);	// Phe
	for(m=0; m<5;m++) {
		for(int n=0; n<8;n++) {
			substr[n] = lineBuf[m*8+n];
		}
		pGEOM[index++] =(float)atof(substr);
	}

	fin1.getline(lineBuf, LINE_LENGTH);	// Met
	for(m=0; m<7;m++) {
		for(int n=0; n<8;n++) {
			substr[n] = lineBuf[m*8+n];
		}
		pGEOM[index++] =(float)atof(substr);
	}
	fin1.close( );

	// read "main_out.MET_SINGLE_MINIMIZE" 
	ifstream fin2(singleEnergyMinimizeResultFile);
	while(!fin2.eof()) {
		fin2.getline(lineBuf, LINE_LENGTH);
		if(!strncmp(lineBuf, " Total Number of Energy Evals. =", 30)) {
			char cEvals[20] = {0}; 
			for(int i=0; i<20; i++) {cEvals[i] = lineBuf[32+i];}
			energyEvalCount = atof(cEvals);

		}
		else if(!strncmp(lineBuf, "Energy Evals. per sec =", 20)) {
			char cSecs[20] = {0}; 
			for(int i=0; i<20; i++) {cSecs[i] = lineBuf[23+i];}
			energyEvalSec = atof(cSecs);
			energyEvalSec = energyEvalCount/ energyEvalSec;
		}
	}
	
	fin2.close( );

	return ETOT;
}

void EnergyCalculator::singleEnergyMinimizeClearFiles( ) {
	ofstream outFile1(singleEnergyMinimizeSourceFile, ios::out);
	outFile1.flush();
	outFile1.close();

	ofstream outFile2(singleEnergyMinimizeResultFile, ios::out);
	outFile2.flush();
	outFile2.close();

	ofstream outFile3(singleEnergyMinimizeOUTOFile, ios::out);
	outFile3.flush();
	outFile3.close();

	ofstream outFile4(logFile, ios::out);
	outFile4.flush();
	outFile4.close();
}






int	EnergyCalculator::multipleEnergyEvaluateWriteSourceFile(double **pGEOM, int popSize, int angleLength) {

	ofstream mainFile(multipleEnergyEvaluateSourceFile, ios::out);
	mainFile << "$CNTRL" << endl;
	mainFile << "runtyp = energy" << endl;
	
	if(isOutputPDBFormat()) 
	{
		mainFile << "PRINT_CART" << endl;
		mainFile << "OUTFORMAT = PDB" << endl;
		mainFile << "FILE  = " <<pdbFile << endl;
	}

	mainFile << "res_code= one_letter" << endl;
	mainFile << "EMINIMA = 0.100E+35" << endl;
	mainFile << "$END" << endl << endl;

	mainFile << "$SEQ" << endl;
	mainFile << "H" << endl;
	mainFile << "YGGFM" << endl;     
	mainFile << "O" << endl;
	mainFile << "$END" << endl << endl;

	mainFile << "$ENERCALC" << endl;
	mainFile << "READ_CONF" << endl;
	mainFile << "$END" << endl << endl;

	mainFile << "$GEOM" << endl;
	mainFile << endl << endl << endl << endl <<endl <<endl << endl;
	mainFile << "$END" << endl;

	mainFile.close();


	ofstream inputFile(multipleEnergyEvaluateInputFile, ios::out);
	for(int i=0; i<popSize; i++) {
		inputFile << i+1 << endl;
		inputFile << "   1  20   6   6   5  11  11" << endl;
		inputFile << endl;
		eceppFormatWrite(inputFile, pGEOM[i], angleLength);
		inputFile << endl;
	}

	inputFile.close();

	return 1;
}




int	EnergyCalculator::multipleEnergyEvaluateReadOUTOFile(double** arrayGeom, int arraySize, double* resultValue) {
	for(int i=0; i<arraySize; i++) resultValue[i] = bigNumber;

	const int LINE_LENGTH = 1000;
	const int WORD_LENGTH = 100;
	char lineBuf[LINE_LENGTH] = {0};
	char wordBuf[WORD_LENGTH] = {0};

	ifstream fin(multipleEnergyEvaluateOUTOFile);
	while(!fin.eof()) {
		fin >> wordBuf;
		int index = atoi(wordBuf); if (index==0) break;
		fin >> wordBuf;
		double value = atof(wordBuf);
		resultValue[index-1] = value;
		for(i=0; i<9; i++) fin.getline(lineBuf, LINE_LENGTH);
	}

	fin.close();
	return 1;
}


int	EnergyCalculator::multipleEnergyEvaluateReadMainResultFile( ) {
	energyEvalCount = 0;
	energyEvalSec = 0;
	
	const int LINE_LENGTH = 1200;
	char lineBuf[LINE_LENGTH] = {0};

	ifstream fin(multipleEnergyEvaluateResultFile);
	
	while(!fin.eof()) {
		fin.getline(lineBuf, LINE_LENGTH);
		if(!strncmp(lineBuf, " Total Number of Energy Evals. =", 30)) {
			char cEvals[20] = {0}; 
			for(int i=0; i<20; i++) {cEvals[i] = lineBuf[32+i];}
			energyEvalCount = atof(cEvals);

		}
		else if(!strncmp(lineBuf, "Energy Evals. per sec =", 20)) {
			char cSecs[20] = {0}; 
			for(int i=0; i<20; i++) {cSecs[i] = lineBuf[23+i];}
			energyEvalSec = atof(cSecs);
			energyEvalSec = energyEvalCount/ energyEvalSec;
		}
	}

	fin.close();
	return 1;
}

void EnergyCalculator::multipleEnergyEvaluateClearFiles( ) {
	ofstream outFile1(multipleEnergyEvaluateSourceFile, ios::out);
	outFile1.flush();
	outFile1.close();

	ofstream outFile2(multipleEnergyEvaluateInputFile, ios::out);
	outFile2.flush();
	outFile2.close();

	ofstream outFile3(multipleEnergyEvaluateResultFile, ios::out);
	outFile3.flush();
	outFile3.close();

	ofstream outFile4(multipleEnergyEvaluateOUTOFile, ios::out);
	outFile4.flush();
	outFile4.close();

	ofstream outFile5(logFile, ios::out);
	outFile5.flush();
	outFile5.close();

}

int	EnergyCalculator::multipleEnergyMinimizeWriteSourceFile(double **arrayGeom, int arraySize, int angleLength, bool* spec) {
	ofstream mainFile(multipleEnergyMinimizeSourceFile, ios::out);

	mainFile << "$CNTRL" << endl;
	mainFile << "runtyp = minimize" << endl;
	if(isOutputPDBFormat()) 
	{
		mainFile << "PRINT_CART" << endl;
		mainFile << "OUTFORMAT = PDB" << endl;
		mainFile << "FILE  = " <<pdbFile << endl;
	}
	mainFile << "res_code= one_letter" << endl;
	mainFile << "EMINIMA = 0.100E+35" << endl;
	mainFile << "var_angles = spec" << endl;
	mainFile << "$END" << endl << endl;

	mainFile << "$SEQ" << endl;
	mainFile << "H" << endl;
	mainFile << "YGGFM" << endl;
	mainFile << "O" << endl;
	mainFile << "$END" << endl << endl;

	mainFile << "$ENERCALC" << endl;
	mainFile << "READ_CONF" << endl;
	mainFile << "$END" << endl << endl;
	
	
	mainFile <<"$SPEC" << endl;
	
	int varNumber = 0;
	int i = 0;

	mainFile << "2" << "  ";
	for(i=0; i<6; i++) {if(spec[i]) varNumber++;}
	mainFile << varNumber << "  ";
	for(i=0; i<6; i++) {if(spec[i]) mainFile << i+1 << "  ";}
	mainFile << endl;

	varNumber = 0;
	mainFile << "3" << "  ";
	for(i=0; i<3; i++) {if(spec[i+6]) varNumber++;}
	mainFile << varNumber << "  ";
	for(i=0; i<3; i++) {if(spec[i+6]) mainFile << i+1 << "  ";}
	mainFile << endl;

	varNumber = 0;
	mainFile << "4" << "  ";
	for(i=0; i<3; i++) {if(spec[i+9]) varNumber++;}
	mainFile << varNumber << "  ";
	for(i=0; i<3; i++) {if(spec[i+9]) mainFile << i+1 << "  ";}
	mainFile << endl;

	varNumber = 0;
	mainFile << "5" << "  ";
	for(i=0; i<5; i++) {if(spec[i+12]) varNumber++;}
	mainFile << varNumber << "  ";
	for(i=0; i<5; i++) {if(spec[i+12]) mainFile << i+1 << "  ";}
	mainFile << endl;

	varNumber = 0;
	mainFile << "6" << "  ";
	for(i=0; i<7; i++) {if(spec[i+17]) varNumber++;}
	mainFile << varNumber << "  ";
	for(i=0; i<7; i++) {if(spec[i+17]) mainFile << i+1 << "  ";}
	mainFile << endl;

	mainFile <<"$END" << endl << endl;

	mainFile << "$GEOM" << endl;
	mainFile << endl << endl << endl << endl <<endl <<endl << endl;
	mainFile << "$END" << endl;
	
	mainFile.close();


	ofstream inputFile(multipleEnergyMinimizeInputFile, ios::out);
	for(i=0; i<arraySize; i++) {
		inputFile << i+1 << endl;
		inputFile << "   1  20   6   6   5  11  11" << endl;
		inputFile << endl;
		eceppFormatWrite(inputFile, arrayGeom[i], angleLength);
		inputFile << endl;
	}

	inputFile.close();

	return 1;
}

int	EnergyCalculator::multipleEnergyMinimizeReadOUTOFile(double** arrayGeom, int arraySize, double* resultValue){
	for(int i=0; i<arraySize; i++) resultValue[i] = bigNumber;

	const int LINE_LENGTH = 1000;
	const int WORD_LENGTH = 100;
	char lineBuf[LINE_LENGTH] = {0};
	char wordBuf[WORD_LENGTH] = {0};
	char substr[9]={0};

	ifstream fin(multipleEnergyMinimizeOUTOResultFile);
	while(!fin.eof()) {
		fin >> wordBuf;
		int index = atoi(wordBuf); if (index==0) break;
		fin >> wordBuf;
		double value = atof(wordBuf);
		resultValue[index-1] = value;
		
		fin.getline(lineBuf, LINE_LENGTH);
		fin.getline(lineBuf, LINE_LENGTH);
		fin.getline(lineBuf, LINE_LENGTH);
		
		
		double* pGEOM = arrayGeom[index-1];
		int  geneIndex = 0;

		fin.getline(lineBuf, LINE_LENGTH);	// Tyr
		for(int m=0; m<6;m++) {
			for(int n=0; n<8;n++) {
				substr[n] = lineBuf[m*8+n];
			}
			pGEOM[geneIndex++] = (float)atof(substr);
		}

		fin.getline(lineBuf, LINE_LENGTH);	// Gly1
		for(m=0; m<3;m++) {
			for(int n=0; n<8;n++) {
				substr[n] = lineBuf[m*8+n];
			}
			pGEOM[geneIndex++] = (float)atof(substr);
		}

		fin.getline(lineBuf, LINE_LENGTH);	// Gly2
		for(m=0; m<3;m++) {
			for(int n=0; n<8;n++) {
				substr[n] = lineBuf[m*8+n];
			}
			pGEOM[geneIndex++] = (float)atof(substr);
		}

		fin.getline(lineBuf, LINE_LENGTH);	// Phe
		for(m=0; m<5;m++) {
			for(int n=0; n<8;n++) {
				substr[n] = lineBuf[m*8+n];
			}
			pGEOM[geneIndex++] = (float)atof(substr);
		}

		fin.getline(lineBuf, LINE_LENGTH);	// Met
		for(m=0; m<7;m++) {
			for(int n=0; n<8;n++) {
				substr[n] = lineBuf[m*8+n];
			}
			pGEOM[geneIndex++] = (float)atof(substr);
		}

		fin.getline(lineBuf, LINE_LENGTH);
	}

	fin.close();

	return 1;
}


int	EnergyCalculator::multipleEnergyMinimizeReadMainResultFile( ) {
	energyEvalCount = 0;
	energyEvalSec = 0;
	
	const int LINE_LENGTH = 1200;
	char lineBuf[LINE_LENGTH] = {0};

	ifstream fin(multipleEnergyMinimizeMainResultFile);
	
	while(!fin.eof()) {
		fin.getline(lineBuf, LINE_LENGTH);
		if(!strncmp(lineBuf, " Total Number of Energy Evals. =", 30)) {
			char cEvals[20] = {0}; 
			for(int i=0; i<20; i++) {cEvals[i] = lineBuf[32+i];}
			energyEvalCount = atof(cEvals);

		}
		else if(!strncmp(lineBuf, "Energy Evals. per sec =", 20)) {
			char cSecs[20] = {0}; 
			for(int i=0; i<20; i++) {cSecs[i] = lineBuf[23+i];}
			energyEvalSec = atof(cSecs);
			energyEvalSec = energyEvalCount/ energyEvalSec;
		}
	}

	fin.close();
	return 1;

}

void EnergyCalculator::multipleEnergyMinimizeClearFiles( ){
	ofstream outFile1(multipleEnergyMinimizeSourceFile, ios::out);
	outFile1.flush();
	outFile1.close();

	ofstream outFile2(multipleEnergyMinimizeInputFile, ios::out);
	outFile2.flush();
	outFile2.close();

	ofstream outFile3(multipleEnergyMinimizeMainResultFile, ios::out);
	outFile3.flush();
	outFile3.close();

	ofstream outFile4(multipleEnergyMinimizeOUTOResultFile, ios::out);
	outFile4.flush();
	outFile4.close();

	ofstream outFile5(logFile, ios::out);
	outFile5.flush();
	outFile5.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////


