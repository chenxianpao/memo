#include "stdafx.h"
#include "DifferentialEvolution.h"
#include "EnergyCalculator.h"
#include "iomanip.h"
#include "math.h"
#include "rmsd.h"
#include "EnkephalinPDB.h"








CDifferentialEvolution::CDifferentialEvolution(const GAGenome& genome) : GAGeneticAlgorithm(genome) {
	mutatorFactor	= 0.5;									 // a good initial choice used by R.Storn
	crossFactor		= 0.1;									 // a good initial choice used by R.Storn

	energyEvals			= 0.0;
	energyEvalsCPUTime	= 0.0;
	index = 0;

	// Buildup initialize
	preSeed = 0;
	averageDistanceOfFirstBank = 0.0;
	numberMinimized = 0.0;
	numberN = 0.0;
}

CDifferentialEvolution::~CDifferentialEvolution() {

}

GABoolean CDifferentialEvolution::TerminateUponMaxMin(GAGeneticAlgorithm& algorithm) {
	float difference = algorithm.population().max() - algorithm.population().min();
	if (difference <= EPSILON) {
		return gaTrue; 
	}
	else { 
		return gaFalse;
	}
}

bool CDifferentialEvolution::better(const GAGenome& firstGenome, const GAGenome& secondGenome) {
	bool flag;
	if (minimaxi() == MAXIMIZE) {
		if (firstGenome.score() > secondGenome.score()) 
			flag = true; 
		else 
			flag = false;
	}
	else{
		if (firstGenome.score() < secondGenome.score()) 
			flag = true; 
		else 
			flag = false;
	}
	return flag;
}

void CDifferentialEvolution::generatePDB(double *pGEOM, int length){
	writeSingleEnergyEvaluatePDBSourceFile(SINGLE_ENERGY_EVALUATE_SOURCE_FILE, pGEOM, length);
	system("i.bat ENERGY singleEnergyEvaluateSourceFile MET_SINGLE_EVALUATE x x  1>> log 2> err");
	clearSingleEnergyEvluateResultFiles();
}


void CDifferentialEvolution::writeSingleEnergyEvaluatePDBSourceFile(const char* sourceFile, double *geomData, int length){
	ofstream outFile(sourceFile, ios::out);
	
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = energy" << endl;
	outFile << "res_code= one_letter" << endl;
	outFile << "PRINT_CART" << endl;
	outFile << "OUTFORMAT   =PDB" << endl;
	outFile << "FILE  = Met" << endl;
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
}



double CDifferentialEvolution::singleEnergyEvaluate(double *pGEOM, int length) {
	double ETOT = 1E25;
	
	writeSingleEnergyEvaluateSourceFile(SINGLE_ENERGY_EVALUATE_SOURCE_FILE, pGEOM, length);
	system("i.bat ENERGY singleEnergyEvaluateSourceFile MET_SINGLE_EVALUATE x x  1>> log 2> err");
	ETOT = readSingleEnergyEvaluateMainOutResultFile(SINGLE_ENERGY_EVALUATE_MAINOUT_FILE);
	clearSingleEnergyEvluateResultFiles();
	return ETOT;
}


void CDifferentialEvolution::singleEnergyEvaluate(GAGenome& genome) {
	GARealGenome& realGenome = dynamic_cast<GARealGenome&>(genome);
	int length = realGenome.length();
	double *pGEOM = new double[length];
	realGenome.score((float)singleEnergyEvaluate(pGEOM, length));
	delete [] pGEOM;
}




void CDifferentialEvolution::writeSingleEnergyEvaluateSourceFile(const char* sourceFile, double *geomData, int length){
	ofstream outFile(sourceFile, ios::out);
	
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = energy" << endl;
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
}

double CDifferentialEvolution::readSingleEnergyEvaluateMainOutResultFile(const char* resultFile){
	double evals		= 0.0;
	double evalsTime	= 0.0;
	double ETOT			= 1E25;

	const int LINE_LENGTH  = 1000; 
	char wordBuf[LINE_LENGTH]={0};
	ifstream fin(resultFile);
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

			fin >> wordBuf;

			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			if(strcmp(wordBuf, "=")==0) {
				fin >> wordBuf;
				evals = atof(wordBuf);
			}else {
				cout << "eval read error!" << endl;
				cout.flush();
				getchar();
			}

			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			if(strcmp(wordBuf,"=")==0) {
				fin >> wordBuf;
				evalsTime = atof(wordBuf);
			} else{
				cout << "eval sec read error!" << endl;
				cout.flush();
				getchar();
			}
			break;
		}
	} while(!fin.eof());

	fin.clear();
	fin.close();

	energyEvals = energyEvals+evals;
	if(evalsTime<0.000001) {
		cout << "Read Error!" << endl;
		cout.flush();
		getchar();
	}else {
		energyEvalsCPUTime = energyEvalsCPUTime + evals/evalsTime;
	}

	return ETOT;
}

void CDifferentialEvolution::clearSingleEnergyEvluateResultFiles( ) {
	const char *mainoutFile	= SINGLE_ENERGY_EVALUATE_MAINOUT_FILE;
	const char *sourceFile	= SINGLE_ENERGY_EVALUATE_SOURCE_FILE;
	const char *logFile		= LOG_FILE;

	ofstream outFile1(mainoutFile, ios::out);
	outFile1.flush();
	outFile1.close();

	ofstream outFile2(sourceFile, ios::out);
	outFile2.flush();
	outFile2.close();

	ofstream outFile3(logFile, ios::out);
	outFile3.flush();
	outFile3.close();
}



void CDifferentialEvolution::singleEnergyMinimize(GAGenome& genome){
	writeSingleEnergyMinimizeSourceFile(SINGLE_ENERGY_MINIMIZE_SOURCE_FILE, genome);
	system("i.bat MINIMIZE singleEnergyMinimizeSourceFile MET_SINGLE_MIN x x  1>> log 2> err");
	readSingleEnergyMinimizeResultOUTOFile(SINGLE_ENERGY_MINIMIZE_OUTO_FILE, genome);
	readSingleEnergyMinimizeResultEvalsFile(SINGLE_ENERGY_MINIMIZE_MAINOUT_FILE);
	clearSingleEnergyMinimizeResultFiles();
}
void CDifferentialEvolution::singleEnergyMinimizeGeneUnchanged(GAGenome& genome){
	writeSingleEnergyMinimizeSourceFile(SINGLE_ENERGY_MINIMIZE_SOURCE_FILE, genome);
	system("i.bat MINIMIZE singleEnergyMinimizeSourceFile MET_SINGLE_MIN x x  1>> log 2> err");
	readSingleEnergyMinimizeResultOUTOFileGeneUnchanged(SINGLE_ENERGY_MINIMIZE_OUTO_FILE, genome);
	readSingleEnergyMinimizeResultEvalsFile(SINGLE_ENERGY_MINIMIZE_MAINOUT_FILE);
	clearSingleEnergyMinimizeResultFiles();
}

double CDifferentialEvolution::singleEnergyMinimize(double* pGEOM, int length) {
	double ETOT = 1E25;
	writeSingleEnergyMinimizeSourceFile(SINGLE_ENERGY_MINIMIZE_SOURCE_FILE, pGEOM, length);
	system("i.bat MINIMIZE singleEnergyMinimizeSourceFile MET_SINGLE_MIN x x  1>> log 2> err");
	readSingleEnergyMinimizeResultEvalsFile(SINGLE_ENERGY_MINIMIZE_MAINOUT_FILE);
	ETOT = readSingleEnergyMinimizeResultETOT(SINGLE_ENERGY_MINIMIZE_OUTO_FILE, pGEOM, length);
	clearSingleEnergyMinimizeResultFiles();
	return ETOT;
}

double CDifferentialEvolution::singleEnergyMinimizeGeneUnchanged(double* pGEOM, int length) {
	double ETOT = 1E25;
	writeSingleEnergyMinimizeSourceFile(SINGLE_ENERGY_MINIMIZE_SOURCE_FILE, pGEOM, length);
	system("i.bat MINIMIZE singleEnergyMinimizeSourceFile MET_SINGLE_MIN x x  1>> log 2> err");
	readSingleEnergyMinimizeResultEvalsFile(SINGLE_ENERGY_MINIMIZE_MAINOUT_FILE);
	ETOT = readSingleEnergyMinimizeResultETOTGeneUnchanged(SINGLE_ENERGY_MINIMIZE_OUTO_FILE);
	clearSingleEnergyMinimizeResultFiles();
	return ETOT;
}

double CDifferentialEvolution::singleEnergyMinimize(double *pGEOM, bool *spec, int length) {
	double ETOT = 1E25;
	writeSingleEnergyMinimizeSpecVarAngles(SINGLE_ENERGY_MINIMIZE_SOURCE_FILE, spec, pGEOM, length);
	system("i.bat MINIMIZE singleEnergyMinimizeSourceFile MET_SINGLE_MIN x x  1>> log 2> err");
	readSingleEnergyMinimizeResultEvalsFile(SINGLE_ENERGY_MINIMIZE_MAINOUT_FILE);
	ETOT = readSingleEnergyMinimizeResultETOT(SINGLE_ENERGY_MINIMIZE_OUTO_FILE, pGEOM, length);
	clearSingleEnergyMinimizeResultFiles();
	return ETOT;

}




void CDifferentialEvolution::writeSingleEnergyMinimizeSourceFile(const char* sourceFile, GAGenome& genome) {
	GARealGenome& realGenome = dynamic_cast<GARealGenome&>(genome);
	int length = realGenome.length();

	ofstream outFile(sourceFile, ios::out);

	outFile << "$CNTRL" << endl;
	outFile << "runtyp = minimize" << endl;
	outFile << "res_code= one_letter" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	outFile << "YGGFM" << endl;
	outFile << "O" << endl;
	outFile << "$END" << endl << endl;

	double * geomData = new double[length];
	for(int k=0; k<length; k++) {
		geomData[k] = realGenome.gene(k);
	}
	outFile << "$GEOM" << endl << endl;
	for(int i=0; i<6;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=6; i<9;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=9; i<12;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=12; i<17;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=17; i<24;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;
	outFile << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();
	delete [] geomData;
}

void CDifferentialEvolution::writeSingleEnergyMinimizeSourceFile(const char* sourceFile, double* geomData, int length){

	ofstream outFile(sourceFile, ios::out);

	outFile << "$CNTRL" << endl;
	outFile << "runtyp = minimize" << endl;
	outFile << "res_code= one_letter" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	outFile << "YGGFM" << endl;
	outFile << "O" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$GEOM" << endl << endl;
	for(int i=0; i<6;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=6; i<9;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=9; i<12;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=12; i<17;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=17; i<24;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;
	outFile << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();
}

void CDifferentialEvolution::writeSingleEnergyMinimizeSpecVarAngles(const char* sourceFile, bool * spec, double* geomData, int length) {
	ofstream outFile(sourceFile, ios::out);

	outFile << "$CNTRL" << endl;
	outFile << "runtyp = minimize" << endl;
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
				<< geomData[i];
	}
	outFile << endl;

	for(i=6; i<9;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=9; i<12;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=12; i<17;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=17; i<24;i++) {
		outFile	<< setw(8)							// 设定field宽度
				<< setprecision(3)					// 设置小数位置
				<< setiosflags(ios::showpoint)		// keep trailing 0s
				<< setiosflags(ios::fixed)			// 使用这些设置
				<< geomData[i];
	}
	outFile << endl;
	outFile << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();
}



void CDifferentialEvolution::readSingleEnergyMinimizeResultOUTOFile(const char* resultFile, GAGenome& genome){
	const int LINE_LENGTH = 1200; 
	double ETOT = 1E25;
	GARealGenome& realGenome = dynamic_cast<GARealGenome&>(genome);

	ifstream fin(resultFile);

	char wordBuf[200]={0};
	char lineBuf[LINE_LENGTH]={0};
	char substr[9]={0};
	int  geneIndex = 0;

	fin >> wordBuf;
	if(strlen(wordBuf)==0) {
		fin.close();
		return;
	}

	fin >> wordBuf;
	ETOT = atof(wordBuf);

	fin.getline(lineBuf, LINE_LENGTH);
	fin.getline(lineBuf, LINE_LENGTH);
	fin.getline(lineBuf, LINE_LENGTH);

	fin.getline(lineBuf, LINE_LENGTH);	// Tyr
	for(int m=0; m<6;m++) {
		for(int n=0; n<8;n++) {
			substr[n] = lineBuf[m*8+n];
		}
		realGenome.gene(geneIndex++,(float)atof(substr));
	}

	fin.getline(lineBuf, LINE_LENGTH);	// Gly1
	for(m=0; m<3;m++) {
		for(int n=0; n<8;n++) {
			substr[n] = lineBuf[m*8+n];
		}
		realGenome.gene(geneIndex++,(float)atof(substr));
	}

	fin.getline(lineBuf, LINE_LENGTH);	// Gly2
	for(m=0; m<3;m++) {
		for(int n=0; n<8;n++) {
			substr[n] = lineBuf[m*8+n];
		}
		realGenome.gene(geneIndex++,(float)atof(substr));
	}

	fin.getline(lineBuf, LINE_LENGTH);	// Phe
	for(m=0; m<5;m++) {
		for(int n=0; n<8;n++) {
			substr[n] = lineBuf[m*8+n];
		}
		realGenome.gene(geneIndex++,(float)atof(substr));
	}

	fin.getline(lineBuf, LINE_LENGTH);	// Met
	for(m=0; m<7;m++) {
		for(int n=0; n<8;n++) {
			substr[n] = lineBuf[m*8+n];
		}
		realGenome.gene(geneIndex++,(float)atof(substr));
	}
	fin.getline(lineBuf, LINE_LENGTH);

	realGenome.score((float)ETOT);

	fin.close();
}

void CDifferentialEvolution::readSingleEnergyMinimizeResultOUTOFileGeneUnchanged(const char* resultFile, GAGenome& genome) {
	const int LINE_LENGTH = 1200; 
	double ETOT = 1E25;
	GARealGenome& realGenome = dynamic_cast<GARealGenome&>(genome);

	ifstream fin(resultFile);

	char wordBuf[200]={0};
	char lineBuf[LINE_LENGTH]={0};
	char substr[9]={0};
	int  geneIndex = 0;

	fin >> wordBuf;
	if(strlen(wordBuf)==0) {
		fin.close();
		return;
	}

	fin >> wordBuf;
	ETOT = atof(wordBuf);
	realGenome.score((float)ETOT);

	fin.close();
}

double CDifferentialEvolution::readSingleEnergyMinimizeResultETOT(const char* resultFile, double* pGEOM, int length) {
	const int LINE_LENGTH = 1200; 
	double ETOT = 1E25;

	ifstream fin(resultFile);

	char wordBuf[200]={0};
	char lineBuf[LINE_LENGTH]={0};
	char substr[9]={0};
	int  geneIndex = 0;

	fin >> wordBuf;
	if(strlen(wordBuf)==0) {
		fin.close();
		return ETOT;
	}

	fin >> wordBuf;
	ETOT = atof(wordBuf);

	fin.getline(lineBuf, LINE_LENGTH);
	fin.getline(lineBuf, LINE_LENGTH);
	fin.getline(lineBuf, LINE_LENGTH);

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

	fin.close();
	return ETOT;
}

double CDifferentialEvolution::readSingleEnergyMinimizeResultETOTGeneUnchanged(const char* resultFile) {
	const int LINE_LENGTH = 1200; 
	double ETOT = 1E25;

	ifstream fin(resultFile);

	char wordBuf[200]={0};
	char lineBuf[LINE_LENGTH]={0};
	char substr[9]={0};
	int  geneIndex = 0;

	fin >> wordBuf;
	if(strlen(wordBuf)==0) {
		fin.close();
		return ETOT;
	}

	fin >> wordBuf;
	ETOT = atof(wordBuf);
	fin.close();
	return ETOT;
}


void CDifferentialEvolution::readSingleEnergyMinimizeResultEvalsFile(const char* mainOutFile) {
	double evals = 0.0;
	double evalsTime = 0.0;

	const int LINE_LENGTH = 1000; 
	char wordBuf[LINE_LENGTH]={0};
	ifstream fin(mainOutFile);

	do{
		fin	>>  wordBuf;
		// Total Number of Energy Evals. =  0.13400E+03
		if(strcmp(wordBuf, "Total")==0) {
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			if(strcmp(wordBuf, "=")==0) {
				fin >> wordBuf;
				evals = atof(wordBuf);
				//cout << evals << endl;
				//cout.flush();
			}else {
				cout << "Read format error!" << endl;
				cout.flush();
				getchar();
			}
			// Energy Evals. per sec =  0.22568E+04
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;
			if(strcmp(wordBuf, "=")==0) {
				fin >> wordBuf;
				evalsTime = atof(wordBuf);
				//cout << evalsTime << endl;
				//cout.flush();
			}else {
				cout << "Read format error!" << endl;
				cout.flush();
				getchar();
			}
			break;
		}

	}while(!fin.eof());

	fin.clear();
	fin.close();

	energyEvals = energyEvals+evals;
	if(evalsTime<0.00001) {
		cout << "Unable to minimize!" << endl;
		cout.flush();
	}else {
		energyEvalsCPUTime = energyEvalsCPUTime + evals/evalsTime;
	}
}

void CDifferentialEvolution::clearSingleEnergyMinimizeResultFiles( ) {
	const char *outoFile	= SINGLE_ENERGY_MINIMIZE_OUTO_FILE;
	const char *mainoutFile	= SINGLE_ENERGY_MINIMIZE_MAINOUT_FILE;
	const char *sourceFile	= SINGLE_ENERGY_MINIMIZE_SOURCE_FILE;
	const char *logFile		= LOG_FILE;

	ofstream outFile1(outoFile, ios::out);
	outFile1.flush();
	outFile1.close();

	ofstream outFile2(mainoutFile, ios::out);
	outFile2.flush();
	outFile2.close();

	ofstream outFile3(sourceFile, ios::out);
	outFile3.flush();
	outFile3.close();

	ofstream outFile4(logFile, ios::out);
	outFile4.flush();
	outFile4.close();
}



void CDifferentialEvolution::multipleEnergyEvaluate(GAPopulation& p) {
	writeMetEnkephalinMultipleEvaluateSourceFile(MULTIPLE_ENERGY_EVALUATION_SOURCE_FILE);
	writeMetEnkephalinMultipleInputFile(MULTIPLE_ENERGY_EVALUATION_INPUT_FILE, p);
	system("i.bat ENERGY metEnkephalinEnergyMultipleEvaluationSourceFile MET_MULTI_EVAL MultipleEnergyEvaluateInput rr  1>> log 2> err");
	readMetEnkephalinMultipleResultOUTOFile("outo.MET_MULTI_EVAL", p);
	energyEvals = energyEvals+p.size();
	clearEnkephalinMultipleEvaluateResultFile();
}

void CDifferentialEvolution::multipleEnergyEvaluate(double** arrayGeom, int arraySize, int dihedralAngleSize) {
	writeMetEnkephalinMultipleEvaluateSourceFile(MULTIPLE_ENERGY_EVALUATION_SOURCE_FILE);
	writeMetEnkephalinMultipleInputFile(MULTIPLE_ENERGY_EVALUATION_INPUT_FILE, arrayGeom, arraySize, dihedralAngleSize);
	system("i.bat ENERGY metEnkephalinEnergyMultipleEvaluationSourceFile MET_MULTI_EVAL MultipleEnergyEvaluateInput rr  1>> log 2> err");
	readMetEnkephalinMultipleResultOUTOFile("outo.MET_MULTI_EVAL", arrayGeom, arraySize, dihedralAngleSize);
	clearEnkephalinMultipleEvaluateResultFile();
}


void CDifferentialEvolution::multipleEnergyEvaluateWritePDB(GAPopulation& p, char* pdbFileName){
	writeMetEnkephalinMultipleEvaluateSourceFile_WritePDB(MULTIPLE_ENERGY_EVALUATION_SOURCE_FILE, pdbFileName);
	writeMetEnkephalinMultipleInputFile(MULTIPLE_ENERGY_EVALUATION_INPUT_FILE, p);
	system("i.bat ENERGY metEnkephalinEnergyMultipleEvaluationSourceFile MET_MULTI_EVAL MultipleEnergyEvaluateInput rr  1>> log 2> err");
	readMetEnkephalinMultipleResultOUTOFile("outo.MET_MULTI_EVAL", p);
	clearEnkephalinMultipleEvaluateResultFile();
}


void CDifferentialEvolution::multipleEnergyEvaluateWithOneResidue(GAPopulation& p, char* residueName) {
	writeMetEnkephalinMultipleEvaluateSourceFile(MULTIPLE_ENERGY_EVALUATION_SOURCE_FILE, residueName);
	writeMetEnkephalinMultipleInputFile(MULTIPLE_ENERGY_EVALUATION_INPUT_FILE, p, residueName);
	system("i.bat ENERGY metEnkephalinEnergyMultipleEvaluationSourceFile MET_MULTI_EVAL MultipleEnergyEvaluateInput rr  1>> log 2> err");
	readMetEnkephalinMultipleResultOUTOFile("outo.MET_MULTI_EVAL", p, residueName);
	clearEnkephalinMultipleEvaluateResultFile();
}



void CDifferentialEvolution::writeMetEnkephalinMultipleEvaluateSourceFile(const char* sourceFile, char* residueName) {
	ofstream outFile(sourceFile, ios::out);
	
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = energy" << endl;
	outFile << "res_code= one_letter" << endl;
	outFile << "EMINIMA = 0.100E+35" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	if(strcmp(residueName, "Tyr") == 0) {
		outFile << "Y" << endl;     
	}
	if(strcmp(residueName, "Gly") == 0) {
		outFile << "G" << endl;     
	}

	if(strcmp(residueName, "Phe") == 0) {
		outFile << "F" << endl;     
	}

	if(strcmp(residueName, "Met") == 0) {
		outFile << "M" << endl;     
	}
	outFile << "O" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$ENERCALC" << endl;
	outFile << "READ_CONF" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$GEOM" << endl;
	outFile << endl << endl << endl; 
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();
}

void CDifferentialEvolution::writeMetEnkephalinMultipleInputFile(const char* inputFile, GAPopulation& p, char* residueName){
	
	ofstream outFile(inputFile, ios::out);

	GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(0));
	int genomeLength = realGenome.length();
	double *geomData = new double[genomeLength];

	for(int i=0; i<p.size(); i++) {
		GARealGenome& realGenomei = dynamic_cast<GARealGenome&>(p.individual(i));
		for(int j=0; j<genomeLength; j++) {
			geomData[j] = realGenomei.gene(j);
		}	

		outFile << i+1 << endl;
		if(strcmp(residueName, "Tyr")==0) {
		outFile << "   1  20  11" << endl;
		}
		if(strcmp(residueName, "Gly")==0) {
		outFile << "   1  6  11" << endl;
		}
		if(strcmp(residueName, "Phe")==0) {
		outFile << "   1  5  11" << endl;
		}
		if(strcmp(residueName, "Met")==0) {
		outFile << "   1  11  11" << endl;
		}
		outFile << endl;

		if(strcmp(residueName, "Tyr")==0) {
			for(int k=0; k<6; k++) {
				outFile	<< setw(8)						// 设定field宽度
						<< setprecision(3)				// 设置小数位置
						<< setiosflags(ios::showpoint)	// keep trailing 0s
						<< setiosflags(ios::fixed)		// 使用这些设置
						<< geomData[k];
			}
			outFile << endl;
		}

		if(strcmp(residueName, "Gly")==0) {
			for(int k=6; k<9; k++) {
				outFile	<< setw(8)						// 设定field宽度
						<< setprecision(3)				// 设置小数位置
						<< setiosflags(ios::showpoint)	// keep trailing 0s
						<< setiosflags(ios::fixed)		// 使用这些设置
						<< geomData[k];
			}
			outFile << endl;
		}


		if(strcmp(residueName, "Phe")==0) {
			for(int k=12; k<17; k++) {
				outFile	<< setw(8)						// 设定field宽度
						<< setprecision(3)				// 设置小数位置
						<< setiosflags(ios::showpoint)	// keep trailing 0s
						<< setiosflags(ios::fixed)		// 使用这些设置
						<< geomData[k];
			}
			outFile << endl; 
		}


		if(strcmp(residueName, "Met")==0) {
			for(int k=17; k<24; k++) {
				outFile	<< setw(8)						// 设定field宽度
						<< setprecision(3)				// 设置小数位置
						<< setiosflags(ios::showpoint)	// keep trailing 0s
						<< setiosflags(ios::fixed)		// 使用这些设置
						<< geomData[k];
			}
			outFile << endl; 
		}

		outFile << endl;
		outFile.flush();
	}

	delete [] geomData;
	outFile.close();
}


void CDifferentialEvolution::readMetEnkephalinMultipleResultOUTOFile(const char *resultFile, GAPopulation& p, char *residueName) {
	const int LINE_LENGTH = 1200; 
	int	popSize = p.size();
	long testLine = 0;
	
	ifstream fin(resultFile);
	char wordBuf[50]={0};
	char lineBuf[LINE_LENGTH]={0};
	char substr[9]={0};

	for(int i=0; i<popSize; i++) {
		int geneIndex = 0;
		float ETOT = 1E25;

		fin >> wordBuf;
		int index = atoi(wordBuf)-1; 
		if(index < 0 &&  i==0) {
			do {
				cout << "File Empty, Error!" << endl;
				cout.flush();
			}while(1);
		}
		if(index <0 && i>0) {
			GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(i));
			realGenome.score(ETOT);
		} else {
			while(i != index) {
				GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(i));
				realGenome.score(ETOT);
				i++;
			}
			GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(index));

			fin >> wordBuf;
			ETOT = (float)atof(wordBuf);
			fin.getline(lineBuf, LINE_LENGTH);
			fin.getline(lineBuf, LINE_LENGTH);
			fin.getline(lineBuf, LINE_LENGTH);

			if(strcmp(residueName, "Tyr")==0) {
				geneIndex = 0;
				fin.getline(lineBuf, LINE_LENGTH);	// Tyr
				for(int m=0; m<6;m++) {
					for(int n=0; n<8;n++) {
						substr[n] = lineBuf[m*8+n];
					}
					realGenome.gene(geneIndex++,(float)atof(substr));
				}
			}

			if(strcmp(residueName, "Gly") == 0) {
				geneIndex = 6;
				fin.getline(lineBuf, LINE_LENGTH);	// Gly1
				for(int m=0; m<3;m++) {
					for(int n=0; n<8;n++) {
						substr[n] = lineBuf[m*8+n];
					}
					realGenome.gene(geneIndex++,(float)atof(substr));
				}
			}

			if(strcmp(residueName, "Phe") == 0) {
				geneIndex = 12;
				fin.getline(lineBuf, LINE_LENGTH);	// Phe
				for(int m=0; m<5;m++) {
					for(int n=0; n<8;n++) {
						substr[n] = lineBuf[m*8+n];
					}
					realGenome.gene(geneIndex++,(float)atof(substr));
				}
			}

			if(strcmp(residueName, "Met") == 0) {
				geneIndex = 17;
				fin.getline(lineBuf, LINE_LENGTH);	// Met
				for(int m=0; m<7;m++) {
					for(int n=0; n<8;n++) {
						substr[n] = lineBuf[m*8+n];
					}
					realGenome.gene(geneIndex++,(float)atof(substr));
				}
			}

			fin.getline(lineBuf, LINE_LENGTH);
			realGenome.score((float)ETOT);
		}
	}

	fin.close();
}


void CDifferentialEvolution::multipleEnergyMinimizeWithOneResidue(GAPopulation& p, char* residueName) {
	writeMetEnkephalinMultipleMinimizeSourceFile(MULTIPLE_ENERGY_MINIMIZE_SOURCE_FILE, residueName);
	writeMetEnkephalinMultipleInputFile(MULTIPLE_ENERGY_MINIMIZE_INPUT_FILE, p, residueName);
	system("i.bat MINIMIZE metEnkephalinEnergyMultipleMinimizeSourceFile MET_MULTI_MIN MultipleEnergyMinimizeInput rr  1>> log 2> err");
	readMetEnkephalinMultipleResultOUTOFile("outo.MET_MULTI_MIN", p, residueName);
	clearEnkephalinMultipleMinimizeResultFile();
}

void CDifferentialEvolution::writeMetEnkephalinMultipleMinimizeSourceFile(const char* sourceFile, char* residueName) {
	ofstream outFile(sourceFile, ios::out);
	
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = minimize" << endl;
	outFile << "res_code= one_letter" << endl;
	outFile << "EMINIMA = 0.100E+35" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	if(strcmp(residueName, "Tyr") == 0) {
		outFile << "Y" << endl;     
	}
	if(strcmp(residueName, "Gly") == 0) {
		outFile << "G" << endl;     
	}

	if(strcmp(residueName, "Phe") == 0) {
		outFile << "F" << endl;     
	}

	if(strcmp(residueName, "Met") == 0) {
		outFile << "M" << endl;     
	}
	outFile << "O" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$ENERCALC" << endl;
	outFile << "READ_CONF" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$GEOM" << endl;
	outFile << endl << endl << endl; 
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();
}

void CDifferentialEvolution::multipleEnergyMinimize(GAPopulation& p) {
	writeMetEnkephalinMultipleMinimizeSourceFile(MULTIPLE_ENERGY_MINIMIZE_SOURCE_FILE);
	writeMetEnkephalinMultipleInputFile(MULTIPLE_ENERGY_MINIMIZE_INPUT_FILE, p);
	system("i.bat MINIMIZE metEnkephalinEnergyMultipleMinimizeSourceFile MET_MULTI_MIN MultipleEnergyMinimizeInput rr  1>> log 2> err");
	readMetEnkephalinMultipleResultOUTOFile("outo.MET_MULTI_MIN", p);
	readMetEnkephalinEnergyEvals(MULTIPLE_ENERGY_MINIMIZE_RESULT_FILE);
	clearEnkephalinMultipleMinimizeResultFile();
}


void CDifferentialEvolution::multipleEnergyMinimizeWritePDB(GAPopulation& p, char* pdbFileName){
	writeMetEnkephalinMultipleMinimizeSourceFile_WritePDB(MULTIPLE_ENERGY_MINIMIZE_SOURCE_FILE, pdbFileName);
	writeMetEnkephalinMultipleInputFile(MULTIPLE_ENERGY_MINIMIZE_INPUT_FILE, p);
	system("i.bat MINIMIZE metEnkephalinEnergyMultipleMinimizeSourceFile MET_MULTI_MIN MultipleEnergyMinimizeInput rr  1>> log 2> err");
	readMetEnkephalinMultipleResultOUTOFile("outo.MET_MULTI_MIN", p);
	readMetEnkephalinEnergyEvals(MULTIPLE_ENERGY_MINIMIZE_RESULT_FILE);
	clearEnkephalinMultipleMinimizeResultFile();
}

void CDifferentialEvolution::multipleEnergyMinimize(double** arrayGeom, int arraySize, int dihedralAngleSize) {
	writeMetEnkephalinMultipleMinimizeSourceFile(MULTIPLE_ENERGY_MINIMIZE_SOURCE_FILE);
	writeMetEnkephalinMultipleInputFile(MULTIPLE_ENERGY_MINIMIZE_INPUT_FILE, arrayGeom, arraySize, dihedralAngleSize);
	system("i.bat MINIMIZE metEnkephalinEnergyMultipleMinimizeSourceFile MET_MULTI_MIN MultipleEnergyMinimizeInput rr  1>> log 2> err");
	readMetEnkephalinMultipleResultOUTOFile("outo.MET_MULTI_MIN", arrayGeom, arraySize, dihedralAngleSize);
	clearEnkephalinMultipleMinimizeResultFile();
}


void CDifferentialEvolution::multipleEnergyMinimize(double** arrayGeom, bool *spec, int arraySize, int dihedralAngleSize) {
	writeMetEnkephalinMultipleMinimizeSourceFile(MULTIPLE_ENERGY_MINIMIZE_SOURCE_FILE, spec);
	writeMetEnkephalinMultipleInputFile(MULTIPLE_ENERGY_MINIMIZE_INPUT_FILE, arrayGeom, arraySize, dihedralAngleSize);
	system("i.bat MINIMIZE metEnkephalinEnergyMultipleMinimizeSourceFile MET_MULTI_MIN MultipleEnergyMinimizeInput rr  1>> log 2> err");
	readMetEnkephalinMultipleResultOUTOFile("outo.MET_MULTI_MIN", arrayGeom, arraySize, dihedralAngleSize);
	clearEnkephalinMultipleMinimizeResultFile();
}


void CDifferentialEvolution::multipleEnergyMinimizeFixedGene(GAPopulation& p) {
	writeMetEnkephalinMultipleMinimizeSourceFile(MULTIPLE_ENERGY_MINIMIZE_SOURCE_FILE);
	writeMetEnkephalinMultipleInputFile(MULTIPLE_ENERGY_MINIMIZE_INPUT_FILE, p);
	system("i.bat MINIMIZE metEnkephalinEnergyMultipleMinimizeSourceFile MET_MULTI_MIN MultipleEnergyMinimizeInput rr  1>> log 2> err");
	readMetEnkephalinMultipleResultOUTOFileGeneUnchanged("outo.MET_MULTI_MIN", p);
	readMetEnkephalinEnergyEvals(MULTIPLE_ENERGY_MINIMIZE_RESULT_FILE);
	clearEnkephalinMultipleMinimizeResultFile();
}





void CDifferentialEvolution::writeMetEnkephalinMultipleEvaluateSourceFile(const char* sourceFile){
	ofstream outFile(sourceFile, ios::out);
	
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = energy" << endl;
	outFile << "res_code= one_letter" << endl;
	outFile << "EMINIMA = 0.100E+35" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	outFile << "YGGFM" << endl;     
	outFile << "O" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$ENERCALC" << endl;
	outFile << "READ_CONF" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$GEOM" << endl;
	outFile << endl << endl << endl << endl <<endl <<endl << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();
}


void CDifferentialEvolution::writeMetEnkephalinMultipleEvaluateSourceFile_WritePDB(const char* sourceFile, char* pdbFileName){
	ofstream outFile(sourceFile, ios::out);
	
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = energy" << endl;
	outFile << "PRINT_CART" << endl;
	outFile << "OUTFORMAT = PDB" << endl;
	outFile << "FILE = ";  
	outFile << pdbFileName << endl;
	outFile << "res_code= one_letter" << endl;
	outFile << "EMINIMA = 0.100E+35" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	outFile << "YGGFM" << endl;
	outFile << "O" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$ENERCALC" << endl;
	outFile << "READ_CONF" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$GEOM" << endl;
	outFile << endl << endl << endl << endl <<endl <<endl << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();
}

void CDifferentialEvolution::writeMetEnkephalinMultipleMinimizeSourceFile(const char* sourceFile){
	ofstream outFile(sourceFile, ios::out);
	
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = minimize" << endl;
	outFile << "res_code= one_letter" << endl;
	outFile << "EMINIMA = 0.100E+35" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	outFile << "YGGFM" << endl; 

	outFile << "O" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$ENERCALC" << endl;
	outFile << "READ_CONF" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$GEOM" << endl;
	outFile << endl << endl << endl << endl <<endl <<endl << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();
}


void CDifferentialEvolution::writeMetEnkephalinMultipleMinimizeSourceFile(const char* sourceFile, bool *spec) {
	ofstream outFile(sourceFile, ios::out);
	
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = minimize" << endl;
	outFile << "res_code= one_letter" << endl;
	outFile << "var_angles = spec" << endl;
	outFile << "EMINIMA = 0.100E+35" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	outFile << "YGGFM" << endl;
	outFile << "O" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$ENERCALC" << endl;
	outFile << "READ_CONF" << endl;
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

	outFile << "$GEOM" << endl;
	outFile << endl << endl << endl << endl <<endl <<endl << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();


}

void CDifferentialEvolution::writeMetEnkephalinMultipleMinimizeSourceFile_WritePDB(const char* sourceFile, char* pdbFileName){
	ofstream outFile(sourceFile, ios::out);
	
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = minimize" << endl;
	outFile << "PRINT_CART" << endl;
	outFile << "OUTFORMAT = PDB" << endl;
	outFile << "FILE = ";  
	outFile << pdbFileName << endl;
	outFile << "res_code= one_letter" << endl;
	outFile << "EMINIMA = 0.100E+35" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	outFile << "YGGFM" << endl;
	outFile << "O" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$ENERCALC" << endl;
	outFile << "READ_CONF" << endl;
	outFile << "$END" << endl << endl;

	outFile << "$GEOM" << endl;
	outFile << endl << endl << endl << endl <<endl <<endl << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();
}


void CDifferentialEvolution::writeMetEnkephalinMultipleInputFile(const char* inputFile, GAPopulation& p){
	ofstream outFile(inputFile, ios::out);

	GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(0));
	int genomeLength = realGenome.length();
	double *geomData = new double[genomeLength];

	for(int i=0; i<p.size(); i++) {
		GARealGenome& realGenomei = dynamic_cast<GARealGenome&>(p.individual(i));
		for(int j=0; j<genomeLength; j++) {
			geomData[j] = realGenomei.gene(j);
		}	

		outFile << i+1 << endl;
		outFile << "   1  20   6   6   5  11  11" << endl;
		
		outFile << endl;

		for(int k=0; k<6; k++) {
			outFile	<< setw(8)						// 设定field宽度
					<< setprecision(3)				// 设置小数位置
					<< setiosflags(ios::showpoint)	// keep trailing 0s
					<< setiosflags(ios::fixed)		// 使用这些设置
					<< geomData[k];
		}
		outFile << endl;

		for(k=6; k<9; k++) {
			outFile	<< setw(8)						// 设定field宽度
					<< setprecision(3)				// 设置小数位置
					<< setiosflags(ios::showpoint)	// keep trailing 0s
					<< setiosflags(ios::fixed)		// 使用这些设置
					<< geomData[k];
		}
		outFile << endl;

		for(k=9; k<12; k++) {
			outFile	<< setw(8)						// 设定field宽度
					<< setprecision(3)				// 设置小数位置
					<< setiosflags(ios::showpoint)	// keep trailing 0s
					<< setiosflags(ios::fixed)		// 使用这些设置
					<< geomData[k];
		}
		outFile << endl;

		for(k=12; k<17; k++) {
			outFile	<< setw(8)						// 设定field宽度
					<< setprecision(3)				// 设置小数位置
					<< setiosflags(ios::showpoint)	// keep trailing 0s
					<< setiosflags(ios::fixed)		// 使用这些设置
					<< geomData[k];
		}
		outFile << endl;

		for(k=17; k<24; k++) {
			outFile	<< setw(8)						// 设定field宽度
					<< setprecision(3)				// 设置小数位置
					<< setiosflags(ios::showpoint)	// keep trailing 0s
					<< setiosflags(ios::fixed)		// 使用这些设置
					<< geomData[k];
		}
		outFile << endl;

		outFile << endl;
		outFile.flush();
	}

	delete [] geomData;
	outFile.close();
}


void CDifferentialEvolution::writeMetEnkephalinMultipleInputFile(const char* inputFile, 
																 double ** arrayGeom, 
																 int arraySize, 
																 int variableSize){

	ofstream outFile(inputFile, ios::out);

	for(int i=0; i<arraySize; i++) {
		double* geomData = arrayGeom[i];

		outFile << i+1 << endl;
		outFile << "   1  20   6   6   5  11  11" << endl;
		
		outFile << endl;

		for(int k=0; k<6; k++) {
			outFile	<< setw(8)						// 设定field宽度
					<< setprecision(3)				// 设置小数位置
					<< setiosflags(ios::showpoint)	// keep trailing 0s
					<< setiosflags(ios::fixed)		// 使用这些设置
					<< geomData[k];
		}
		outFile << endl;

		for(k=6; k<9; k++) {
			outFile	<< setw(8)						// 设定field宽度
					<< setprecision(3)				// 设置小数位置
					<< setiosflags(ios::showpoint)	// keep trailing 0s
					<< setiosflags(ios::fixed)		// 使用这些设置
					<< geomData[k];
		}
		outFile << endl;

		for(k=9; k<12; k++) {
			outFile	<< setw(8)						// 设定field宽度
					<< setprecision(3)				// 设置小数位置
					<< setiosflags(ios::showpoint)	// keep trailing 0s
					<< setiosflags(ios::fixed)		// 使用这些设置
					<< geomData[k];
		}
		outFile << endl;

		for(k=12; k<17; k++) {
			outFile	<< setw(8)						// 设定field宽度
					<< setprecision(3)				// 设置小数位置
					<< setiosflags(ios::showpoint)	// keep trailing 0s
					<< setiosflags(ios::fixed)		// 使用这些设置
					<< geomData[k];
		}
		outFile << endl;

		for(k=17; k<24; k++) {
			outFile	<< setw(8)						// 设定field宽度
					<< setprecision(3)				// 设置小数位置
					<< setiosflags(ios::showpoint)	// keep trailing 0s
					<< setiosflags(ios::fixed)		// 使用这些设置
					<< geomData[k];
		}
		outFile << endl;

		outFile << endl;
		outFile.flush();

	}

	outFile.close();
}


void CDifferentialEvolution::readMetEnkephalinMultipleResultOUTOFile(const char *resultFile, GAPopulation& p){
	const int LINE_LENGTH = 1200; 
	int	popSize = p.size();
	long testLine = 0;
	
	ifstream fin(resultFile);
	char wordBuf[50]={0};
	char lineBuf[LINE_LENGTH]={0};
	char substr[9]={0};

	for(int i=0; i<popSize; i++) {
		int geneIndex = 0;
		float ETOT = 1E25;

		fin >> wordBuf;
		int index = atoi(wordBuf)-1; 
		if(index < 0 &&  i==0) {
			do {
				cout << "File Empty, Error!" << endl;
				cout.flush();
			}while(1);
		}
		if(index <0 && i>0) {
			GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(i));
			realGenome.score(ETOT);
		} else {
			while(i != index) {
				GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(i));
				realGenome.score(ETOT);
				i++;
			}
			GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(index));

			fin >> wordBuf;
			ETOT = (float)atof(wordBuf);
			fin.getline(lineBuf, LINE_LENGTH);
			fin.getline(lineBuf, LINE_LENGTH);
			fin.getline(lineBuf, LINE_LENGTH);

			fin.getline(lineBuf, LINE_LENGTH);	// Tyr
			for(int m=0; m<6;m++) {
				for(int n=0; n<8;n++) {
					substr[n] = lineBuf[m*8+n];
				}
				realGenome.gene(geneIndex++,(float)atof(substr));
			}

			fin.getline(lineBuf, LINE_LENGTH);	// Gly1
			for(m=0; m<3;m++) {
				for(int n=0; n<8;n++) {
					substr[n] = lineBuf[m*8+n];
				}
				realGenome.gene(geneIndex++,(float)atof(substr));
			}

			fin.getline(lineBuf, LINE_LENGTH);	// Gly2
			for(m=0; m<3;m++) {
				for(int n=0; n<8;n++) {
					substr[n] = lineBuf[m*8+n];
				}
				realGenome.gene(geneIndex++,(float)atof(substr));
			}

			fin.getline(lineBuf, LINE_LENGTH);	// Phe
			for(m=0; m<5;m++) {
				for(int n=0; n<8;n++) {
					substr[n] = lineBuf[m*8+n];
				}
				realGenome.gene(geneIndex++,(float)atof(substr));
			}

			fin.getline(lineBuf, LINE_LENGTH);	// Met
			for(m=0; m<7;m++) {
				for(int n=0; n<8;n++) {
					substr[n] = lineBuf[m*8+n];
				}
				realGenome.gene(geneIndex++,(float)atof(substr));
			}
			fin.getline(lineBuf, LINE_LENGTH);

			realGenome.score((float)ETOT);
		}
	}

	fin.close();
}



void CDifferentialEvolution::readMetEnkephalinMultipleResultOUTOFileGeneUnchanged(const char *resultFile, GAPopulation& p) {
	const int LINE_LENGTH = 1200; 
	int	popSize = p.size();

	ifstream fin(resultFile);
	char wordBuf[50]={0};
	char lineBuf[LINE_LENGTH]={0};
	char substr[9]={0};

	for(int i=0; i<popSize; i++) {
		int geneIndex = 0;
		float ETOT = 1E25;

		fin >> wordBuf;
		int index = atoi(wordBuf)-1; 
		if(index < 0 &&  i==0) {
			do {
				cout << "File Empty, Error!" << endl;
				cout.flush();
			}while(1);
		}
		if(index <0 && i>0) {
			GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(i));
			realGenome.score(ETOT);
		} else {
			while(i != index) {
				GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(i));
				realGenome.score(ETOT);
				i++;
			}
			GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p.individual(index));

			fin >> wordBuf;
			ETOT = (float)atof(wordBuf);
			fin.getline(lineBuf, LINE_LENGTH);	// ETOT line
			fin.getline(lineBuf, LINE_LENGTH);	// sequence line 
			fin.getline(lineBuf, LINE_LENGTH);	// end group-H line
			fin.getline(lineBuf, LINE_LENGTH);	// Tyr line
			fin.getline(lineBuf, LINE_LENGTH);	// Gly1 line
			fin.getline(lineBuf, LINE_LENGTH);	// Gly2 line
			fin.getline(lineBuf, LINE_LENGTH);	// Phe line
			fin.getline(lineBuf, LINE_LENGTH);	// Met line
			fin.getline(lineBuf, LINE_LENGTH);	// eng group -COOH

			realGenome.score((float)ETOT);
		}
	}

	fin.close();
}


void CDifferentialEvolution::readMetEnkephalinMultipleResultOUTOFile(const char *resultFile, 
																	 double **arrayGeom, 
																	 int arraySize, 
																	 int variableSize) {
	const int WORD_LENGTH = 50;
	const int LINE_LENGTH = 1200; 
	char wordBuf[WORD_LENGTH] = {0};
	char lineBuf[LINE_LENGTH] = {0};
	char substr[9] = {0};


	ifstream fin(resultFile);

	for(int i=0; i<arraySize; i++) {
		fin >> wordBuf;
		int index = atoi(wordBuf)-1;
		if(index == i) {
			fin >> wordBuf;
			double ETOT = (float)atof(wordBuf);
			arrayGeom[i][variableSize] = ETOT;
			arrayGeom[i][variableSize+1] = 1;	// flag the ETOT is correct!
			
			fin.getline(lineBuf, LINE_LENGTH);	// ETOT line
			fin.getline(lineBuf, LINE_LENGTH);	// sequence line 
			fin.getline(lineBuf, LINE_LENGTH);	// end group-H line
			
			int geneIndex = 0;

			fin.getline(lineBuf, LINE_LENGTH);	// Tyr
			for(int m=0; m<6;m++) {
				for(int n=0; n<8;n++) {
					substr[n] = lineBuf[m*8+n];
				}
				arrayGeom[i][geneIndex++] = (float)atof(substr);
			}

			fin.getline(lineBuf, LINE_LENGTH);	// Gly1
			for(m=0; m<3;m++) {
				for(int n=0; n<8;n++) {
					substr[n] = lineBuf[m*8+n];
				}
				arrayGeom[i][geneIndex++] = (float)atof(substr);
			}

			fin.getline(lineBuf, LINE_LENGTH);	// Gly2
			for(m=0; m<3;m++) {
				for(int n=0; n<8;n++) {
					substr[n] = lineBuf[m*8+n];
				}
				arrayGeom[i][geneIndex++] = (float)atof(substr);
			}

			fin.getline(lineBuf, LINE_LENGTH);	// Phe
			for(m=0; m<5;m++) {
				for(int n=0; n<8;n++) {
					substr[n] = lineBuf[m*8+n];
				}
				arrayGeom[i][geneIndex++] = (float)atof(substr);
			}

			fin.getline(lineBuf, LINE_LENGTH);	// Met
			for(m=0; m<7;m++) {
				for(int n=0; n<8;n++) {
					substr[n] = lineBuf[m*8+n];
				}
				arrayGeom[i][geneIndex++] = (float)atof(substr);
			}
			
			fin.getline(lineBuf, LINE_LENGTH);	// eng group -COOH
		} else {
			// can't get the correct ETOT (maybe occures in minimized procedure)
			arrayGeom[i][variableSize] = 1E25;
			arrayGeom[i][variableSize+1] = -1;	// flag the ETOT is not correct!
		}

	}

	fin.close();
}


void CDifferentialEvolution::readMetEnkephalinEnergyEvals(const char* mainOutFile){
	double evals = 0.0;
	double evalsTime = 0.0;
	
	const int LINE_LENGTH = 1000; 
	ifstream fin(mainOutFile);
	char flagLine[LINE_LENGTH] = "          All iterations done.";
	char wordBuf[50]={0};


	bool find= false;
	do{
		char str[LINE_LENGTH] = {0};
		fin.getline(str, LINE_LENGTH);
		if(strcmp(str, flagLine)==0) {
			find = true;
			break;
		}
	} while (!fin.eof());   

	bool firstEntry = true;
	if(find) {
		do {
			fin >> wordBuf;
			if(strcmp(wordBuf, "=")==0) {
				if(firstEntry) {
					fin>>wordBuf;
					evals = atof(wordBuf);
					firstEntry = false;
				}else {
					fin>>wordBuf;
					evalsTime = atof(wordBuf);
				}
			}
		}while(!fin.eof());
	} else {
		cout << "Failure, iteration abnomral terminated!" << endl;
		cout.flush();
		getchar();


	}
	if(evals==0 || evalsTime==0) {
		cout << "Read Error, Overflow" << endl;
		cout.flush();
		getchar();
	}

	fin.close();

	energyEvals = energyEvals+evals;
	energyEvalsCPUTime = energyEvalsCPUTime + evals/evalsTime;
}



void CDifferentialEvolution::clearEnkephalinMultipleEvaluateResultFile( ){
	const char *outoFile	= "outo.MET_MULTI_EVAL";
	const char *mainoutFile	= "main_out.MET_MULTI_EVAL";
	const char *sourceFile	= MULTIPLE_ENERGY_EVALUATION_SOURCE_FILE;
	const char *inputFile	= MULTIPLE_ENERGY_EVALUATION_INPUT_FILE;
	const char *logFile		= "log";

	ofstream outFile1(outoFile, ios::out);
	outFile1.flush();
	outFile1.close();

	ofstream outFile2(mainoutFile, ios::out);
	outFile2.flush();
	outFile2.close();

	ofstream outFile3(sourceFile, ios::out);
	outFile3.flush();
	outFile3.close();

	ofstream outFile4(inputFile, ios::out);
	outFile4.flush();
	outFile4.close();

	ofstream outFile5(logFile, ios::out);
	outFile5.flush();
	outFile5.close();
}

void CDifferentialEvolution::clearEnkephalinMultipleMinimizeResultFile( ){
	const char *outoFile	= "outo.MET_MULTI_MIN";
	const char *mainoutFile	= "main_out.MET_MULTI_MIN";
	const char *sourceFile	= MULTIPLE_ENERGY_MINIMIZE_SOURCE_FILE;
	const char *inputFile	= MULTIPLE_ENERGY_MINIMIZE_INPUT_FILE;
	const char *logFile		= "log";

	ofstream outFile1(outoFile, ios::out);
	outFile1.flush();
	outFile1.close();

	ofstream outFile2(mainoutFile, ios::out);
	outFile2.flush();
	outFile2.close();

	ofstream outFile3(sourceFile, ios::out);
	outFile3.flush();
	outFile3.close();

	ofstream outFile4(inputFile, ios::out);
	outFile4.flush();
	outFile4.close();

	ofstream outFile5(logFile, ios::out);
	outFile5.flush();
	outFile5.close();
}

void CDifferentialEvolution::persistConformationData(const GAPopulation& p, char* mainFileName, int generation) {
	char charGeneration[10]	= {0};
	itoa(generation, charGeneration, 10);
	char generationsFile  [100] = {0};
	char resultSummaryFile[100] = {0};

	strcpy(generationsFile, mainFileName);
	strcat(generationsFile, charGeneration);
	strcat(generationsFile, ".txt");
	strcpy(resultSummaryFile,mainFileName);
	strcat(resultSummaryFile, "Summary.txt");

	ofstream outfile(generationsFile, ios::out);
	for(int i=0; i<p.size(); i++) {
		GARealGenome& trial = dynamic_cast<GARealGenome&>(p.individual(i));
		double fscore = trial.score();
		outfile << i+1			<< "\t"
			    << fscore		<< "\t";
		for(int j = 0; j<trial.length(); j++) {
			outfile << trial.gene(j) << "\t";
		}
		outfile << "\n";
	}
	outfile.close();


	ofstream outSummaryfile(resultSummaryFile, ios::app);

	GARealGenome& bestGenome = dynamic_cast<GARealGenome&>(p.best());
	double bestScore = bestGenome.score();
	if(generation==1) outSummaryfile << "\n";

	outSummaryfile	<< generation<< "\t"
					<< bestScore << "\t"
					<< energyEvals << "\t"
					<< energyEvalsCPUTime << "\t"
					<< stats.numrep	
					<< "\t";

	for(int j = 0; j<bestGenome.length(); j++) {
		outSummaryfile << bestGenome.gene(j) << "\t";
	}
	outSummaryfile << "\n";
	outSummaryfile.close();
}

void CDifferentialEvolution::persistAllEvolutionData(const GAPopulation& p, char* allPopulationFileName) {
	ofstream outfile(allPopulationFileName, ios::app);

	for(int i=0; i<p.size(); i++) {
		GARealGenome& trial = dynamic_cast<GARealGenome&>(p.individual(i));
		double fscore = trial.score();
		outfile << i+1			<< "\t"
			    << fscore		<< "\t";
		for(int j = 0; j<trial.length(); j++) {
			outfile << trial.gene(j) << "\t";
		}
		outfile << "\n";
	}

	outfile.close();
}


void CDifferentialEvolution::mergeConformationData(const char *mainFileName, int fileSize, int fileNumber) {
	char mergeFile[100] = {0};
	strcpy(mergeFile, mainFileName);
	strcat(mergeFile, "Merge.txt");

	char tempC[100] = {0};
	long index = 0;

	ofstream outMergeFile(mergeFile, ios::out);

	for(int i=1; i<=fileNumber; i++) {
		char charFileNumber[10]	= {0};
		itoa(i,charFileNumber,10);
		char inputFile[100] = {0};
		strcpy(inputFile, mainFileName);
		strcat(inputFile, charFileNumber);
		strcat(inputFile, ".txt");
		
		ifstream fin(inputFile);
		for(int j=1; j<=fileSize; j++) {
			double ETOT = 0;
			double geom[24] = {0};

			fin >> tempC;
			fin >> tempC;	// object value;
			ETOT = (double)atof(tempC);
			for(int k=0; k<24; k++) {
				fin>>tempC;
				geom[k] = (double)atof(tempC);
			}

			outMergeFile << ++index << "\t" << ETOT << "\t";
			for(k=0; k<24; k++) 
				outMergeFile << geom[k] << "\t";
			outMergeFile << "\n";
		}
		fin.close();
		cout << i << endl;

	}

	outMergeFile.close( );
	cout << "Files have been merged!" << "\n";
}

void CDifferentialEvolution::objectiveData(const GAEvalData& v){
	for(int i=0; i<pop->size(); i++)
		pop->individual(i).evalData(v);
}


void CDifferentialEvolution::initialize(unsigned int seed) {
	energyEvals	= 0.0;
	energyEvalsCPUTime = 0.0;

	GARandomSeed( );
	pop->initialize( );

	multipleEnergyMinimize(*pop);
	//multipleEnergyEvaluate(*pop);
	pop->evaluate(gaFalse);
	stats.reset(*pop);
	if(TCONDITION==0) terminator(TerminateUponMaxMin);
}

/**********************************************************************************************/
// Basic DE without Energy Minimization
/**********************************************************************************************/
void CDifferentialEvolution::step() {
	int popSize = pop->size();
	int	genomeLength = (dynamic_cast<GARealGenome*>(&(pop->individual(0))))->length();
	
	GAPopulation *trialPopulation = pop->clone(); 
	for (int i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));
		
		// mutation
		l= GARandomInt(0, popSize-1); while (bi==i || bi==ai);
		do ci = GARandomInt(0, popSize-1); while (ci==i || ci==bi || ci==ai);
		stats.numsel += 3;
		a = dynamic_cast<GARealGenome*>(&(pop->individual(ai)));
		b = dynamic_cast<GARealGenome*>(&(pop->individual(bi)));
		c = dynamic_cast<GARealGenome*>(&(pop->individual(ci)));
		for(int k=0; k<genomeLength; k++) {
			float newGene = c->gene(k) + mutatorFactor * (a->gene(k)-b->gene(k));
			float geneLower = trial->alleleset(k).lower();
			float geneUpper = trial->alleleset(k).upper();
			if (newGene<geneLower || newGene>geneUpper) {
				newGene = GARandomFloat(geneLower, geneUpper);
			}
			trial->gene(k, newGene);
		}
		stats.nummut++;
		
		//crossover
		int j = GARandomInt(0, genomeLength-1);
		for(k=0; k<genomeLength; k++) {
			if ( (GARandomFloat() <= crossFactor) || (k==j) ) {
				stats.numcro++;
			}
			else {
				float newGene = target->gene(k);
				trial->gene(k, newGene);
			}
		}
	}

	// acceptance
	multipleEnergyEvaluate(*trialPopulation);
	
	for( i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));
	
		stats.numeval++; 
		if (better(*trial, *target)) {
			target->copy(*trial);
			stats.numrep++;
		}
	}
	delete trialPopulation;
	
	pop->sort(gaTrue);
	pop->statistics(gaTrue);
	stats.update(*pop);
	
	persistConformationData(*pop, "data\\ConformationData", generation());
}

/**********************************************************************************************/
// Basic DE with Energy Minimization
/**********************************************************************************************/
void CDifferentialEvolution::step1() {
	int popSize = pop->size();
	int	genomeLength = (dynamic_cast<GARealGenome*>(&(pop->individual(0))))->length();
	
	GAPopulation *trialPopulation = pop->clone(); 
	for (int i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));
				
		// mutation
		GARealGenome *a, *b, *c;     
		int	ai, bi, ci; 
		do ai = GARandomInt(0, popSize-1); while (ai==i);
		do bi = GARandomInt(0, popSize-1); while (bi==i || bi==ai);
		do ci = GARandomInt(0, popSize-1); while (ci==i || ci==bi || ci==ai);
		stats.numsel += 3;
		a = dynamic_cast<GARealGenome*>(&(pop->individual(ai)));
		b = dynamic_cast<GARealGenome*>(&(pop->individual(bi)));
		c = dynamic_cast<GARealGenome*>(&(pop->individual(ci)));
		for(int k=0; k<genomeLength; k++) {
			float newGene = c->gene(k) + mutatorFactor * (a->gene(k)-b->gene(k));
			float geneLower = trial->alleleset(k).lower();
			float geneUpper = trial->alleleset(k).upper();
			if (newGene<geneLower || newGene>geneUpper) {
				newGene = GARandomFloat(geneLower, geneUpper);
			}
			trial->gene(k, newGene);
		}
		stats.nummut++;

		//crossover
		int j = GARandomInt(0, genomeLength-1);
		for(k=0; k<genomeLength; k++) {
			if ( (GARandomFloat() <= crossFactor) || (k==j) ) {
				stats.numcro++;
			}
			else {
				float newGene = target->gene(k);
				trial->gene(k, newGene);
			}
		}
	}

	// acceptance
	multipleEnergyMinimize(*trialPopulation);
	//multipleEnergyEvaluate(*trialPopulation);
	for( i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));
	
		stats.numeval++; 
		if (better(*trial, *target)) {
			target->copy(*trial);
			stats.numrep++;
		}
	}
	delete trialPopulation;
	
	pop->sort(gaTrue);
	pop->statistics(gaTrue);
	stats.update(*pop);
	
	persistConformationData(*pop, "data\\ConformationData", generation());
	persistAllEvolutionData(*pop, "data\\allPopulation.txt");
}



/**********************************************************************************************/
// Basic DE algorithm with Buildup and Energy Minimization Procedure
// ref. New optimization method for conformational energy calculations on polypeptides-conformational space annealing
/**********************************************************************************************/
void CDifferentialEvolution::step2() {
	int popSize = pop->size();
	int	genomeLength = (dynamic_cast<GARealGenome*>(&(pop->individual(0))))->length();
	
	GAPopulation *trialPopulation = pop->clone(); 
	
	// buildup operation (include 0 and 25)
	GARealGenome *seed = dynamic_cast<GARealGenome*>(&(pop->individual(GARandomInt(0, (int)0.5*genomeLength))));

	for (int i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));
				
		// mutation operation
		GARealGenome *a, *b, *c;     
		int	ai, bi, ci; 
		do ai = GARandomInt(0, popSize-1); while (ai==i);
		do bi = GARandomInt(0, popSize-1); while (bi==i || bi==ai);
		do ci = GARandomInt(0, popSize-1); while (ci==i || ci==bi || ci==ai);
		stats.numsel += 3;
		a = dynamic_cast<GARealGenome*>(&(pop->individual(ai)));
		b = dynamic_cast<GARealGenome*>(&(pop->individual(bi)));
		c = dynamic_cast<GARealGenome*>(&(pop->individual(ci)));
		for(int k=0; k<genomeLength; k++) {
			float newGene = c->gene(k) + mutatorFactor * (a->gene(k)-b->gene(k));
			float geneLower = trial->alleleset(k).lower();
			float geneUpper = trial->alleleset(k).upper();
			if (newGene<geneLower || newGene>geneUpper) {
				newGene = GARandomFloat(geneLower, geneUpper);
			}
			trial->gene(k, newGene);
		}
		stats.nummut++;
		
		
		int j = GARandomInt(1, 10);
		if(j<=5) {						// DE algorithm's crossover operator; 50%
			int jindex = GARandomInt(0, genomeLength-1);
			for(k=0; k<genomeLength; k++) {
				if ( (GARandomFloat() <= 0.1) || (k==jindex)) {
					stats.numcro++;
				}
				else {
					float newGene = target->gene(k);
					trial->gene(k, newGene);
				}
			}
		}else if(j<=7) {								// group buildup procedure; 20%
			int group = GARandomInt(1, 8);
			if(group==1) {
				trial->gene( 0, seed->gene( 0));
				trial->gene( 1, seed->gene( 1));
				trial->gene( 2, seed->gene( 2));
			}else if(group==2) {
				trial->gene( 3, seed->gene( 3));
				trial->gene( 4, seed->gene( 4));
				trial->gene( 5, seed->gene( 5));
			}else if(group==3) {
				trial->gene( 6, seed->gene( 6));
				trial->gene( 7, seed->gene( 7));
				trial->gene( 8, seed->gene( 8));
			}else if(group==4) {
				trial->gene( 9, seed->gene( 9));
				trial->gene(10, seed->gene(10));
				trial->gene(11, seed->gene(11));
			}else if(group==5) {
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
			}else if(group==6) {
				trial->gene(15, seed->gene(15));
				trial->gene(16, seed->gene(16));
			}else if(group==7) {
				trial->gene(17, seed->gene(17));
				trial->gene(18, seed->gene(18));
				trial->gene(19, seed->gene(19));
			}else {
				trial->gene(20, seed->gene(20));
				trial->gene(21, seed->gene(21));
				trial->gene(22, seed->gene(22));
				trial->gene(23, seed->gene(23));
			}
		} else {
			int connectGroup = GARandomInt(1,7);		// connect group buildup procedure; 30%
			if(connectGroup == 1) {						// 1,2,3
				trial->gene(0, seed->gene(0));
				trial->gene(1, seed->gene(1));
				trial->gene(2, seed->gene(2));
				trial->gene(3, seed->gene(3));
				trial->gene(4, seed->gene(4));
				trial->gene(5, seed->gene(5));
				trial->gene(6, seed->gene(6));
				trial->gene(7, seed->gene(7));
				trial->gene(8, seed->gene(8));
			}else if(connectGroup == 2) {				// 1,3,4
				trial->gene(0, seed->gene(0));
				trial->gene(1, seed->gene(1));
				trial->gene(2, seed->gene(2));
				trial->gene(6, seed->gene(6));
				trial->gene(7, seed->gene(7));
				trial->gene(8, seed->gene(8));
				trial->gene(9, seed->gene(9));
				trial->gene(10, seed->gene(10));
				trial->gene(11, seed->gene(11));
			}else if(connectGroup == 3) {				// 3,4,5
				trial->gene(6, seed->gene(6));
				trial->gene(7, seed->gene(7));
				trial->gene(8, seed->gene(8));
				trial->gene(9, seed->gene(9));
				trial->gene(10, seed->gene(10));
				trial->gene(11, seed->gene(11));
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
			}else if(connectGroup == 4) {				// 4,5,6
				trial->gene(9, seed->gene(9));
				trial->gene(10, seed->gene(10));
				trial->gene(11, seed->gene(11));
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
				trial->gene(15, seed->gene(15));
				trial->gene(16, seed->gene(16));
			}else if(connectGroup == 5) {				// 4,5,7
				trial->gene(9, seed->gene(9));
				trial->gene(10, seed->gene(10));
				trial->gene(11, seed->gene(11));
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
				trial->gene(17, seed->gene(17));
				trial->gene(18, seed->gene(18));
				trial->gene(19, seed->gene(19));
			}else if(connectGroup == 6) {				// 5,6,7
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
				trial->gene(15, seed->gene(15));
				trial->gene(16, seed->gene(16));
				trial->gene(17, seed->gene(17));
				trial->gene(18, seed->gene(18));
				trial->gene(19, seed->gene(19));
			}else {										// 5,7,8
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
				trial->gene(17, seed->gene(17));
				trial->gene(18, seed->gene(18));
				trial->gene(19, seed->gene(19));
				trial->gene(20, seed->gene(20));
				trial->gene(21, seed->gene(21));
				trial->gene(22, seed->gene(22));
				trial->gene(23, seed->gene(23));
			}
		}
	}

	// acceptance operation
	multipleEnergyMinimize(*trialPopulation);
	//multipleEnergyEvaluate(*trialPopulation);
	for( i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));

		stats.numeval++; 
		if (better(*trial, *target)) {
			target->copy(*trial);
			stats.numrep++;
		}
	}
	
	delete trialPopulation;
	
	pop->sort(gaTrue);
	pop->statistics(gaTrue);
	stats.update(*pop);

	persistConformationData(*pop, "data\\ConformationData", generation());
}




void CDifferentialEvolution::step3() {
	int popSize = pop->size();
	int	genomeLength = (dynamic_cast<GARealGenome*>(&(pop->individual(0))))->length();

	GAPopulation *trialPopulation = pop->clone(); 
	for (int i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));

		// mutation
		GARealGenome *a, *b, *c;    
		GARealGenome *base, *base1, *base2;

		int	ai, bi, ci; 
		do ai = GARandomInt(0, popSize-1); while (ai==i);
		do bi = GARandomInt(0, popSize-1); while (bi==i || bi==ai);
		do ci = GARandomInt(0, popSize-1); while (ci==i || ci==bi || ci==ai);
		stats.numsel += 3;
		a = dynamic_cast<GARealGenome*>(&(pop->individual(ai)));
		b = dynamic_cast<GARealGenome*>(&(pop->individual(bi)));
		c = dynamic_cast<GARealGenome*>(&(pop->individual(ci)));

		if (better(*a,*b)) {
			base = a; 
			base1 = b;
		}	
		else {
			base = b; 
			base1 = a;
		}
		
		if (better(*base,*c)) {
			base2 = c;
		}
		else {
			GARealGenome *temp; temp = base; base = c; base2 = temp;
		}

		for(int k=0; k<genomeLength; k++) {
			float newGene	= base->gene(k) + mutatorFactor * (base1->gene(k)-base2->gene(k));
			float geneLower = trial->alleleset(k).lower();
			float geneUpper = trial->alleleset(k).upper();
			if (newGene < geneLower || newGene > geneUpper) {
				newGene = GARandomFloat(geneLower, geneUpper);
			}
			trial->gene(k, newGene);
		}
		stats.nummut++;

		//crossover
		int j = GARandomInt(0, genomeLength-1);
		for(k=0; k<genomeLength; k++) {
			if ( (GARandomFloat() <= crossFactor) || (k==j) ) {
				stats.numcro++;
			}
			else {
				float newGene = target->gene(k);
				trial->gene(k, newGene);
			}
		}
	}

	// acceptance
	multipleEnergyMinimize(*trialPopulation);
	for( i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));

		stats.numeval++; 
		if (better(*trial, *target)) {
			target->copy(*trial);
			stats.numrep++;
		}
	}
	delete trialPopulation;
	
	pop->sort(gaTrue);
	pop->statistics(gaTrue);
	stats.update(*pop);

	persistConformationData(*pop, "data\\ConformationData", generation());
}

void CDifferentialEvolution::step4() {
	int popSize = pop->size();
	int	genomeLength = (dynamic_cast<GARealGenome*>(&(pop->individual(0))))->length();

	GAPopulation *trialPopulation = pop->clone(); 
	for (int i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));

		// mutation
		GARealGenome *a, *b, *c;     
		int	ai, bi, ci; 
		do ai = GARandomInt(0, popSize-1); while (ai==i);
		do bi = GARandomInt(0, popSize-1); while (bi==i || bi==ai);
		do ci = GARandomInt(0, popSize-1); while (ci==i || ci==bi || ci==ai);
		stats.numsel += 3;
		a = dynamic_cast<GARealGenome*>(&(pop->individual(ai)));
		b = dynamic_cast<GARealGenome*>(&(pop->individual(bi)));
		c = dynamic_cast<GARealGenome*>(&(pop->individual(ci)));
		for(int k=0; k<genomeLength; k++) {
			float newGene   = c->gene(k) + mutatorFactor * (a->gene(k)-b->gene(k));
			float geneLower = trial->alleleset(k).lower();
			float geneUpper = trial->alleleset(k).upper();
			if (newGene < geneLower || newGene > geneUpper) {
				newGene = GARandomFloat(geneLower, geneUpper);
			}
			trial->gene(k, newGene);
		}
		stats.nummut++;

		//crossover
		int j = GARandomInt(0, genomeLength-1);
		for(k=0; k<genomeLength; k++) {
			if ( (GARandomFloat() <= crossFactor) || (k==j) ) {
				stats.numcro++;
			}
			else {
				float newGene = target->gene(k);
				trial->gene(k, newGene);
			}
		}
	}

	// acceptance
	multipleEnergyMinimize(*trialPopulation);
	for( i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));

		stats.numeval++; 
		if (better(*trial, *target)) {
			target->copy(*trial);
			
			double *pGEOM	= new double[24];
			double *outGEOM	= new double[24];
			bool   *spec	= new bool[24];
			for(int m=0; m<24; m++) {
				pGEOM[m]	= trial->gene(m);
				outGEOM[m]	= 0.0;
				spec[m]		= false;
			}
			int index;
			do {
				index = GARandomInt(0, 23);
			}while(index==2 || index == 8|| index == 11 || index == 14 || index == 19);
			spec[index] = true;

			ofstream outf("zgj.txt", ios::out);
			outf << trial->score() << endl;
			outf  << "index=" << index << endl;
			outf << "double Tyr_phi		=" 	<< trial->gene(0) << endl;
			outf << "double Tyr_psi		=" 	<< trial->gene(1) << endl;
			outf << "double Tyr_omega	=" 	<< trial->gene(2) << endl;
			outf << "double Tyr_chi1	=" 	<< trial->gene(3) << endl;
			outf << "double Tyr_chi2	=" 	<< trial->gene(4) << endl;
			outf << "double Tyr_chi3	=" 	<< trial->gene(5) << endl;

			outf << "double Gly1_phi	=" 	<< trial->gene(6) << endl;
			outf << "double Gly1_psi	=" 	<< trial->gene(7) << endl;
			outf << "double Gly1_omega	=" 	<< trial->gene(8) << endl;

			outf << "double Gly2_phi	=" 	<< trial->gene(9) << endl;
			outf << "double Gly2_psi	=" 	<< trial->gene(10) << endl;
			outf << "double Gly2_omega	=" 	<< trial->gene(11) << endl;


			outf << "double Phe_phi		=" 	<< trial->gene(12) << endl;
			outf << "double Phe_psi		=" 	<< trial->gene(13) << endl;
			outf << "double Phe_omega	=" 	<< trial->gene(14) << endl;
			outf << "double Phe_chi1	=" 	<< trial->gene(15) << endl;
			outf << "double Phe_chi2	=" 	<< trial->gene(16) << endl;


			outf << "double Met_phi		=" 	<< trial->gene(17) << endl;
			outf << "double Met_psi		=" 	<< trial->gene(18) << endl;
			outf << "double Met_omega	=" 	<< trial->gene(19) << endl;
			outf << "double Met_chi1	=" 	<< trial->gene(20) << endl;
			outf << "double Met_chi2	=" 	<< trial->gene(21) << endl;
			outf << "double Met_chi3	=" 	<< trial->gene(22) << endl;
			outf << "double Met_chi4	=" 	<< trial->gene(23) << endl;
			outf.close();

			

			double aa = trial->score();
			double bb = target->score();
			double retValue = CAMOptimize(pGEOM, trial->score(), spec, index, 24, outGEOM);
			if (retValue < target->score()) {
				for(m=0; m<24; m++) {
					target->gene(m, outGEOM[m]);
				}
				target->score(retValue);
				stats.numrep++;

			}
			

			delete [] pGEOM;
			delete [] outGEOM;
			delete [] spec;


			//stats.numrep++;
		} else {








		}

	}
	delete trialPopulation;
	
	pop->sort(gaTrue);
	pop->statistics(gaTrue);
	stats.update(*pop);

	persistConformationData(*pop, "data\\ConformationData", generation());
}



/**********************************************************************************************/
// CrowdingDE without Energy Minimization
/**********************************************************************************************/
void CDifferentialEvolution::step5() {
	int popSize = pop->size();
	int	genomeLength = (dynamic_cast<GARealGenome*>(&(pop->individual(0))))->length();
	
	GAPopulation *trialPopulation = pop->clone(); 

// buildup operation (include 0 and 25)
	GARealGenome *seed = dynamic_cast<GARealGenome*>(&(pop->individual(GARandomInt(0, (int)0.5*genomeLength))));
	for (int i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));
		
		// mutation
	    stats.numsel +=3;
		GAGenome *a1, *b1, *c1;

		GAGenome** genomeList = 0;
		double *psum = 0;
		int newSize,removedIndex;
		//generate three target individual;
		genomeList = extractGenomeArrayList(target, i, pop);
		
		a1 = genomeList[0];
		newSize = removeGenomeFromArray(genomeList, popSize-1, removedIndex);
		

		psum = getPrepareArray(genomeList, newSize);
		removedIndex = selectIndex(psum, newSize);
		b1 = genomeList[removedIndex];
		newSize = removeGenomeFromArray(genomeList, newSize, removedIndex);
		delete [] psum;

		psum = getPrepareArray(genomeList, newSize);
		removedIndex = selectIndex(psum, newSize);
		c1 = genomeList[removedIndex];
		newSize = removeGenomeFromArray(genomeList, newSize, removedIndex);
		delete [] psum;

		GARealGenome* a = dynamic_cast<GARealGenome*>(a1);
		GARealGenome* b = dynamic_cast<GARealGenome*>(b1);
		GARealGenome* c = dynamic_cast<GARealGenome*>(c1);
		delete [] genomeList;

		// mutation
		for(int k=0; k<genomeLength; k++) {
			float newgene = c->gene(k) + mutatorFactor*(a->gene(k)-b->gene(k));
			float gLower = trial->alleleset(k).lower();
			float gUpper = trial->alleleset(k).upper();
			if (newgene<gLower || newgene>gUpper) {
				newgene = GARandomFloat(gLower,gUpper);
			}
			trial->gene(k, newgene);
		}
		stats.nummut++;
		
		int j = GARandomInt(1, 10);
		if(j<=5) {						// DE algorithm's crossover operator; 50%
			int jindex = GARandomInt(0, genomeLength-1);
			for(k=0; k<genomeLength; k++) {
				if ( (GARandomFloat() <= 0.1) || (k==jindex)) {
					stats.numcro++;
				}
				else {
					float newGene = target->gene(k);
					trial->gene(k, newGene);
				}
			}
		}else if(j<=7) {								// group buildup procedure; 20%
			int group = GARandomInt(1, 8);
			if(group==1) {
				trial->gene( 0, seed->gene( 0));
				trial->gene( 1, seed->gene( 1));
				trial->gene( 2, seed->gene( 2));
			}else if(group==2) {
				trial->gene( 3, seed->gene( 3));
				trial->gene( 4, seed->gene( 4));
				trial->gene( 5, seed->gene( 5));
			}else if(group==3) {
				trial->gene( 6, seed->gene( 6));
				trial->gene( 7, seed->gene( 7));
				trial->gene( 8, seed->gene( 8));
			}else if(group==4) {
				trial->gene( 9, seed->gene( 9));
				trial->gene(10, seed->gene(10));
				trial->gene(11, seed->gene(11));
			}else if(group==5) {
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
			}else if(group==6) {
				trial->gene(15, seed->gene(15));
				trial->gene(16, seed->gene(16));
			}else if(group==7) {
				trial->gene(17, seed->gene(17));
				trial->gene(18, seed->gene(18));
				trial->gene(19, seed->gene(19));
			}else {
				trial->gene(20, seed->gene(20));
				trial->gene(21, seed->gene(21));
				trial->gene(22, seed->gene(22));
				trial->gene(23, seed->gene(23));
			}
		} else {
			int connectGroup = GARandomInt(1,7);		// connect group buildup procedure; 30%
			if(connectGroup == 1) {						// 1,2,3
				trial->gene(0, seed->gene(0));
				trial->gene(1, seed->gene(1));
				trial->gene(2, seed->gene(2));
				trial->gene(3, seed->gene(3));
				trial->gene(4, seed->gene(4));
				trial->gene(5, seed->gene(5));
				trial->gene(6, seed->gene(6));
				trial->gene(7, seed->gene(7));
				trial->gene(8, seed->gene(8));
			}else if(connectGroup == 2) {				// 1,3,4
				trial->gene(0, seed->gene(0));
				trial->gene(1, seed->gene(1));
				trial->gene(2, seed->gene(2));
				trial->gene(6, seed->gene(6));
				trial->gene(7, seed->gene(7));
				trial->gene(8, seed->gene(8));
				trial->gene(9, seed->gene(9));
				trial->gene(10, seed->gene(10));
				trial->gene(11, seed->gene(11));
			}else if(connectGroup == 3) {				// 3,4,5
				trial->gene(6, seed->gene(6));
				trial->gene(7, seed->gene(7));
				trial->gene(8, seed->gene(8));
				trial->gene(9, seed->gene(9));
				trial->gene(10, seed->gene(10));
				trial->gene(11, seed->gene(11));
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
			}else if(connectGroup == 4) {				// 4,5,6
				trial->gene(9, seed->gene(9));
				trial->gene(10, seed->gene(10));
				trial->gene(11, seed->gene(11));
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
				trial->gene(15, seed->gene(15));
				trial->gene(16, seed->gene(16));
			}else if(connectGroup == 5) {				// 4,5,7
				trial->gene(9, seed->gene(9));
				trial->gene(10, seed->gene(10));
				trial->gene(11, seed->gene(11));
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
				trial->gene(17, seed->gene(17));
				trial->gene(18, seed->gene(18));
				trial->gene(19, seed->gene(19));
			}else if(connectGroup == 6) {				// 5,6,7
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
				trial->gene(15, seed->gene(15));
				trial->gene(16, seed->gene(16));
				trial->gene(17, seed->gene(17));
				trial->gene(18, seed->gene(18));
				trial->gene(19, seed->gene(19));
			}else {										// 5,7,8
				trial->gene(12, seed->gene(12));
				trial->gene(13, seed->gene(13));
				trial->gene(14, seed->gene(14));
				trial->gene(17, seed->gene(17));
				trial->gene(18, seed->gene(18));
				trial->gene(19, seed->gene(19));
				trial->gene(20, seed->gene(20));
				trial->gene(21, seed->gene(21));
				trial->gene(22, seed->gene(22));
				trial->gene(23, seed->gene(23));
			}
		}
	}



	// acceptance
	multipleEnergyMinimize(*trialPopulation);
  //multipleEnergyEvaluate(*trialPopulation);
	for( i=0; i<popSize; i++) {
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));
		GARealGenome *target = dynamic_cast<GARealGenome*>(findNearestIndividual(*trial, *pop));

		stats.numeval++; 
		if (better(*trial, *target)) {
			target->copy(*trial);
			stats.numrep++;
		}
	}
	delete trialPopulation;
	
	pop->sort(gaTrue);
	pop->statistics(gaTrue);
	stats.update(*pop);
	
	persistConformationData(*pop, "data\\ConformationData", generation());
}

float CDifferentialEvolution::distance(const GAGenome& first, const GAGenome& second) {
	double dist = 0.0;
	const GARealGenome &realGenome1 = dynamic_cast<const GARealGenome&>(first);
	const GARealGenome &realGenome2 = dynamic_cast<const GARealGenome&>(second);
	int dim = realGenome1.length();
	for(int i=0; i<dim; i++) {
		double gene1 = realGenome1.gene(i);
		double gene2 = realGenome2.gene(i);
		dist += (gene1-gene2)*(gene1-gene2);
	}

	dist = sqrt(dist);
	return dist;
}

GAGenome* CDifferentialEvolution::findNearestIndividual(const GAGenome& trial, const GAPopulation& population){
	GAGenome *nearestGenome = 0;
	int popSize = population.size();
	
	double dist = distance(trial, population.individual(0));
	nearestGenome = &population.individual(0);
	for(int i=1; i<popSize; i++) {
		double d = distance(trial, population.individual(i));
		if(d<dist) {
			dist = d;
			nearestGenome = &population.individual(i);
		}
	}

	return nearestGenome;
}
void CDifferentialEvolution::quickSortAscendingFitness(GAGenome **c, int l, int r) 
{
	int i,j; 
	float v; 
	GAGenome *t;
	if(r > l)
	{
		v = c[r]->fitness(); 
		i = l-1; 
		j = r;
		for(;;)
		{
			while(c[++i]->fitness() < v && i <= r);
			while(c[--j]->fitness() > v && j > 0);
			if(i >= j) break;
			t = c[i]; 
			c[i] = c[j]; 
			c[j] = t;
		}
		t = c[i]; 
		c[i] = c[r]; 
		c[r] = t;
		quickSortAscendingFitness(c,l,i-1);
		quickSortAscendingFitness(c,i+1,r);
	}
}


void CDifferentialEvolution::quickSortDescendingFitness(GAGenome **c, int l, int r) 
{
	int i,j; 
	float v; 
	GAGenome *t;
	if(r > l)
	{
		v = c[r]->fitness(); 
		i = l-1; 
		j = r;
		for(;;)
		{
			while(c[++i]->fitness() > v && i <= r);  //空循环
			while(c[--j]->fitness() < v && j > 0);
			if(i >= j) break;
			t = c[i]; c[i] = c[j]; c[j] = t;
		}
		t = c[i]; c[i] = c[r]; c[r] = t;
		quickSortDescendingFitness(c,l,i-1);
		quickSortDescendingFitness(c,i+1,r);
	}
}


GAGenome** CDifferentialEvolution::extractGenomeArrayList(GAGenome* parent, int pos, GAPopulation* pop)
{
	int popSize = pop->size();
	int lstSize = popSize-1;

	GAGenome** genomeArray = new GAGenome*[lstSize];
	int ai = 0;
	for(int i=0; i<popSize; i++) {
		GAGenome* ind = &pop->individual(i);
		if(i!=pos) {
			genomeArray[ai++] = ind;
		}
	}

	// 搜索得到最大距离
	double maxDistance = 0.0;
	for(i=0; i<lstSize; i++) {
		GAGenome* ind = genomeArray[i];
		double dist = distance(*parent, *ind);
		ind->fitness(dist);
		if (maxDistance < dist) 
			maxDistance = dist;
	}
	
	// 转换为适应度函数
	float a = 0.5;
	for(i=0; i<lstSize; i++) 
	{
		GAGenome* ind = genomeArray[i];
		double dist = ind->fitness();
		double fit = (maxDistance-dist)/maxDistance;
		fit = pow(fit, a);
		ind->fitness(fit);
	}
	// 按照适应度函数排序（从大到小）
	quickSortDescendingFitness(genomeArray, 0, lstSize-1);

	return genomeArray;
}

// genomeArrayList已按照适应度函数排序（从大到小）
//getPrepareArray:
double* CDifferentialEvolution::getPrepareArray(GAGenome** genomeArrayList, int size) {
	double *psum = new double[size];
	double amax = genomeArrayList[0]->fitness();
	double amin = genomeArrayList[size-1]->fitness( );

	psum[0] = genomeArrayList[0]->fitness();
	for(int i=1; i<size; i++) {
		psum[i] = genomeArrayList[i]->fitness()+psum[i-1];
	}

	for(i=0; i<size; i++) {
		psum[i] /= psum[size-1];
	}

	return psum;
}

GAGenome** CDifferentialEvolution::getGenomeArrayList(GAGenome* parent, GAPopulation *pop)
{
	int popSize = pop->size( );
	int arraySize = popSize-1;
	float a = 0.5;

	GAGenome** genomeArray = new GAGenome*[arraySize];
	int ai = 0;
	for(int i=0; i<popSize; i++) {
		GAGenome* ind = &pop->individual(i);
		if(parent != ind) {
			genomeArray[ai++] = ind;
		}
	}
	
	double maxDistance = 0.0;
	for(i=0; i<arraySize; i++) {
		GAGenome* ind = genomeArray[i];
		double dist = distance(*parent, *ind);
		ind->fitness(dist);
		if (maxDistance < dist) {
			maxDistance = dist;
		}
	}

	for(i=0; i<arraySize; i++) {
		GAGenome* ind = genomeArray[i];
		double dist = ind->fitness();
		double fit = (maxDistance-dist)/maxDistance;
		fit = pow(fit, a);

/*
		double sd = maxDistance/3;
		sd = 2* pow(sd,2);
		double fit = -pow(dist,2)/sd;
		fit = exp(fit);
*/

		ind->fitness(fit);
	}

	quickSortDescendingFitness(genomeArray, 0, arraySize-1);

	return genomeArray;
}

double* CDifferentialEvolution::prepareArray(GAGenome** genomeArrayList, int size) 
{
	ofstream outfile("czh.txt", ios::out);

	double *psum = new double[size];
	double amax = genomeArrayList[0]->fitness();
	double amin = genomeArrayList[size-1]->fitness( );

	for(int i=0; i<size; i++){
		outfile << genomeArrayList[i]->fitness() << endl;
	}
	outfile << endl;
	
	psum[0] = genomeArrayList[0]->fitness();
	for(i=1; i<size; i++) {
		psum[i] = genomeArrayList[i]->fitness()+psum[i-1];
	}
	for(i=0; i<size; i++) {
		psum[i] /= psum[size-1];
	}

	for(i=0; i<size; i++){
		outfile << psum[i] << endl;
	}
	outfile.close();

	return psum;

}



GAGenome* CDifferentialEvolution::select(GAGenome* parent, GAPopulation *pop) 
{
	GAGenome** genomeList = getGenomeArrayList(parent, pop);
	double* psum = prepareArray(genomeList, pop->size()-1);
	
	float cutoff;
	int i, upper, lower;
	
	cutoff = GARandomFloat( );
	lower = 0; 
	upper = pop->size()-2;
	
	while(upper >= lower){
		i = lower + (upper-lower)/2;
		if(psum[i] > cutoff)
			upper = i-1;
		else
			lower = i+1;
	}
	
	lower = GAMin(pop->size()-2, lower);
	lower = GAMax(0, lower);
	
	GAGenome *rg = genomeList[lower];
	delete [] genomeList;
	delete [] psum;
	return rg;
}


GAGenome* CDifferentialEvolution::select(GAGenome** genomeList, int listSize) 
{
	double* psum = prepareArray(genomeList, listSize);
	
	float cutoff;
	int i, upper, lower;
	
	cutoff = GARandomFloat( );
	lower = 0; 
	upper = listSize-1;
	
	while(upper >= lower)
	{
		i = lower + (upper-lower)/2;
		if(psum[i] > cutoff)
			upper = i-1;
		else
			lower = i+1;
	}
	//let the value of lower locate at [0,listSize-1];
	lower = GAMin(listSize-1, lower);
	lower = GAMax(0, lower);
	
	GAGenome *rg = genomeList[lower];
	delete [] genomeList;
	delete [] psum;
	return rg;
}

int CDifferentialEvolution::selectIndex(double* psum, int size) 
{
	float cutoff;
	int i, upper, lower;

	cutoff = GARandomFloat( );
	lower = 0; 
	upper = size-1;

	while(upper >= lower)
	{
		i = lower + (upper-lower)/2;
		if(psum[i] > cutoff)
			upper = i-1;
		else
			lower = i+1;
	}
	//let the value of lower locate at [0,listSize-1];
	lower = GAMin(size-1, lower);
	lower = GAMax(0, lower);
	return lower;
}

//removeGenomeFromeArray:remove the element in the position pos from the array list.
int CDifferentialEvolution::removeGenomeFromArray(GAGenome**list, int lstSize, int pos) 
{
	int rsize = lstSize;

	for(int i=0; i<lstSize; i++) {
		if(i==pos) {
			rsize = rsize-1;
			break;
		}
	}

	if(i<lstSize-1) {
		for(int j=i; j<lstSize-1; j++)
			list[j] = list[j+1];
	}
	return rsize;
}


int CDifferentialEvolution::removeGenomeList(GAGenome** genomeList, int size, GAGenome* removedGenome) 
{
	int rsize = 0;
	for(int i=0; i<size; i++) {
		GAGenome* trial = genomeList[i];
		if(trial==removedGenome) {
			rsize=rsize-1;
			break;
		}
	}

	if(i<size-1) {
		for(int j=i; j<size-1; j++) {
			genomeList[j] = genomeList[j+1];
		}
	}

	return rsize;
}




void CDifferentialEvolution::randomMinimize( ){
	ofstream outfile1("origin.txt",   ios::out);
	ofstream outfile2("minimize.txt", ios::out);
	long index1 = 0;
	long index2 = 0;

	int size = pop->size();
	for(int i=1; i<101; i++) {
		GARandomSeed(0);
		pop->initialize();
		char midName[20] = {0};
		itoa(i, midName, 10);
		
		/*
		multipleEnergyEvaluate(*pop);
		for(int j=0; j<size; j++) {
			GARealGenome& realGenome = dynamic_cast<GARealGenome&>(pop->individual(j));
			outfile1 << ++index1 << "\t" << realGenome.score() << "\t";
			for(int k=0; k<24; k++) {
				outfile1 << realGenome.gene(k) << "\t";
			}
			outfile1<<endl;
		}
		cout << i << endl;
		cout.flush();
		
		multipleEnergyMinimize(*pop);
		for(j=0; j<size; j++) {
			GARealGenome& realGenome = dynamic_cast<GARealGenome&>(pop->individual(j));
			outfile2 << ++index2 << "\t" << realGenome.score() << "\t";
			for(int k=0; k<24; k++) {
				outfile2 << realGenome.gene(k) << "\t";
			}
			outfile2<<endl;
		}
		cout << i << endl;
		cout.flush();
		*/
		
		multipleEnergyEvaluateWithOneResidue(*pop, "Tyr");
		for(int j=0; j<size; j++) {
			GARealGenome& realGenome = dynamic_cast<GARealGenome&>(pop->individual(j));
			outfile1 << ++index1 << "\t" << realGenome.score() << "\t";
			for(int k=0; k<6; k++) {
				outfile1 << realGenome.gene(k) << "\t";
			}
			outfile1<<endl;
		}
		cout << i << endl;
		cout.flush();

		multipleEnergyMinimizeWithOneResidue(*pop, "Tyr");
		for(j=0; j<size; j++) {
			GARealGenome& realGenome = dynamic_cast<GARealGenome&>(pop->individual(j));
			outfile2 << ++index2 << "\t" << realGenome.score() << "\t";
			for(int k=0; k<6; k++) {
				outfile2 << realGenome.gene(k) << "\t";
			}
			outfile2<<endl;
		}
		cout << i << endl;
		cout.flush();

	}


	outfile1.close();
	outfile2.close();
}








void CDifferentialEvolution::save(GAPopulation* p) {
	ofstream outfile("deminimize.txt", ios::app);

	for(int j=0; j<p->size(); j++) {
		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p->individual(j));
		outfile << ++index << "\t" << realGenome.score() << "\t";
		for(int k=0; k<24; k++) {
			outfile << realGenome.gene(k) << "\t";
		}
		outfile<<endl;
		cout<<index << endl;
		cout.flush();
	}
	outfile.close();
}

void CDifferentialEvolution::save(GAPopulation* p, char* fileName) {
	ofstream outfile(fileName, ios::out);

	for(int j=0; j<p->size(); j++) {
		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(p->individual(j));
		outfile << ++index << "\t"  << realGenome.score() << "\t";
		for(int k=0; k<24; k++) {
			outfile << realGenome.gene(k) << "\t";
		}
		outfile<<endl;
		cout.flush();
	}
	
	outfile.close();
}




void CDifferentialEvolution::dihedralAngleTest( ){
	// Reference conformation
	/*
	double Tyr_phi			= -83.466;
	double Tyr_psi			= 155.791;
	double Tyr_omega		= -177.129;
	double Tyr_chi1			= -173.181;
	double Tyr_chi2			= 79.338;
	double Tyr_chi3			= -166.326;

	double Gly1_phi			= -154.281;
	double Gly1_psi			= 85.850;
	double Gly1_omega		= 168.513;

	double Gly2_phi			= 82.965;
	double Gly2_psi			= -75.050;
	double Gly2_omega		= -169.969;

	double Phe_phi			= -136.851;
	double Phe_psi			= 19.106;
	double Phe_omega		= -174.093;
	double Phe_chi1			= 58.858;
	double Phe_chi2			= -85.471;

	double Met_phi			= -163.456;
	double Met_psi			= 160.941;
	double Met_omega		= -179.790;
	double Met_chi1			= 52.868;
	double Met_chi2			= 175.299;
	double Met_chi3			= -179.862;
	double Met_chi4			= -58.586;
	*/

	double Tyr_phi			= 163.14;
	double Tyr_psi			= 165.764;
	double Tyr_omega		= -154.655;
	double Tyr_chi1			= -113.408;
	double Tyr_chi2			= -89.181;
	double Tyr_chi3			= -6.836;

	double Gly1_phi			= -127.537;
	double Gly1_psi			= 91.741;
	double Gly1_omega		= 0.061;

	double Gly2_phi			=-24.793;
	double Gly2_psi			= -44.381;
	double Gly2_omega		= 8.21;

	double Phe_phi			=  122.411;
	double Phe_psi			= -172.159;
	double Phe_omega		= -137.778;
	double Phe_chi1			= -160.149;
	double Phe_chi2			= -119.381;

	double Met_phi			= -149.769;
	double Met_psi			= -7.885;
	double Met_omega		= -111.135;
	double Met_chi1			= 73.624;
	double Met_chi2			= 44.404;
	double Met_chi3			= -89.591;
	double Met_chi4			= 124.795;
	

	int targetVariableIndex	= 23;	
	int arraySize			= 400;
	int variableSize		= 24;
	int stepInterval		= 40;



	double * refGeom = new double[24];
	refGeom[0]	=	Tyr_phi;
	refGeom[1]	=	Tyr_psi;
	refGeom[2]	=	Tyr_omega;
	refGeom[3]	=	Tyr_chi1;
	refGeom[4]	=	Tyr_chi2;
	refGeom[5]	=	Tyr_chi3;

	refGeom[6]	=	Gly1_phi;
	refGeom[7]	=	Gly1_psi;
	refGeom[8]	=	Gly1_omega;

	refGeom[9]  =	Gly2_phi;
	refGeom[10] =	Gly2_psi;
	refGeom[11] =	Gly2_omega;

	refGeom[12] =	Phe_phi;
	refGeom[13] =	Phe_psi;
	refGeom[14] =	Phe_omega;
	refGeom[15] =	Phe_chi1;
	refGeom[16] =	Phe_chi2;

	refGeom[17] =	Met_phi;
	refGeom[18] =	Met_psi;
	refGeom[19] =	Met_omega;
	refGeom[20] =	Met_chi1;
	refGeom[21] =	Met_chi2;
	refGeom[22] =	Met_chi3;
	refGeom[23] =	Met_chi4;


	bool *spec	= new bool[variableSize];
	for(int m=0; m<24; m++) {
		spec[m]	= true;
	}
	spec[targetVariableIndex] = true;

	ofstream outFile("dihedralTest\\angleTest_1.txt", ios::out);

	for(int base=-180; base<180; base=base+stepInterval) {
		double** arrayGeom = new double* [arraySize];
		for(int ai=0; ai<arraySize; ai++) {
			arrayGeom[ai] = new double[variableSize+2];	// geom + ETOT + flag
		}

		// initializing
		for(ai=0; ai<arraySize; ai++) {
			double* pTargetGeom = arrayGeom[ai];
			for(int vs=0; vs<variableSize; vs++) {
				pTargetGeom[vs] = refGeom[vs];
			}
			pTargetGeom[variableSize]	= 0.0;		// ETOT;
			pTargetGeom[variableSize+1] = 0.0;		// is valid?;
		}

		// change variable
		for(ai=0; ai<arraySize; ai++) {
			double* pTargetGeom = arrayGeom[ai];
			pTargetGeom[targetVariableIndex] = base + ai*((double)stepInterval/(double)arraySize);
		}
		
		multipleEnergyEvaluate(arrayGeom, arraySize, variableSize);
		
		for(ai=0; ai<arraySize; ai++) {
			double x = base + ai*((double)stepInterval/(double)arraySize);
			if(arrayGeom[ai][variableSize+1]  == 1) {
				outFile << x << "\t";
				outFile << arrayGeom[ai][variableSize] << "\t";
				outFile << log(arrayGeom[ai][variableSize]+11.7074);	// loge
			} else {
				outFile << x << "\t";
				outFile << "-" <<"\t" << "-" << "\t";
			}
			outFile << "\n";
		}

		for(ai=0; ai<arraySize; ai++) {
			delete [] arrayGeom[ai];
		}
		delete [] arrayGeom;
		cout << base << " is completed!" << endl; cout.flush();
	}
	outFile.close( );


	for(m=0; m<24; m++) {
		spec[m]	= false;
	}
	spec[targetVariableIndex] = true;

	ofstream outFile2("dihedralTest\\angleTest_2.txt", ios::out);
	for(base=-180; base<180; base=base+stepInterval) {
		double** arrayGeom = new double* [arraySize];
		for(int ai=0; ai<arraySize; ai++) {
			arrayGeom[ai] = new double[variableSize+2];	// geom + ETOT + flag
		}

		// initializing
		for(ai=0; ai<arraySize; ai++) {
			double* pTargetGeom = arrayGeom[ai];
			for(int vs=0; vs<variableSize; vs++) {
				pTargetGeom[vs] = refGeom[vs];
			}
			pTargetGeom[variableSize]	= 0.0;		// ETOT;
			pTargetGeom[variableSize+1] = 0.0;		// is valid?;
		}

		// change variable
		for(ai=0; ai<arraySize; ai++) {
			double* pTargetGeom = arrayGeom[ai];
			pTargetGeom[targetVariableIndex] = base + ai*((double)stepInterval/(double)arraySize);
		}
		
		multipleEnergyMinimize(arrayGeom, spec, arraySize, variableSize);
		
		for(ai=0; ai<arraySize; ai++) {
			double x = base + ai*((double)stepInterval/(double)arraySize);
			if(arrayGeom[ai][variableSize+1]  == 1) {
				outFile2 << x << "\t";
				outFile2 << arrayGeom[ai][variableSize] << "\t";
				outFile2 << log(arrayGeom[ai][variableSize]+11.7074);	// loge
			} else {
				outFile2 << x << "\t";
				outFile2 << "-" <<"\t" << "-" << "\t";
			}
			outFile2 << "\n";
		}

		for(ai=0; ai<arraySize; ai++) {
			delete [] arrayGeom[ai];
		}
		delete [] arrayGeom;
		cout << base << " is completed!" << endl; cout.flush();
	}
	outFile2.close( );

	for(m=0; m<24; m++) {
		spec[m]	= true;
	}
	spec[targetVariableIndex] = false;

	ofstream outFile3("dihedralTest\\angleTest_3.txt", ios::out);
	for(base=-180; base<180; base=base+stepInterval) {
		double** arrayGeom = new double* [arraySize];
		for(int ai=0; ai<arraySize; ai++) {
			arrayGeom[ai] = new double[variableSize+2];	// geom + ETOT + flag
		}

		// initializing
		for(ai=0; ai<arraySize; ai++) {
			double* pTargetGeom = arrayGeom[ai];
			for(int vs=0; vs<variableSize; vs++) {
				pTargetGeom[vs] = refGeom[vs];
			}
			pTargetGeom[variableSize]	= 0.0;		// ETOT;
			pTargetGeom[variableSize+1] = 0.0;		// is valid?;
		}

		// change variable
		for(ai=0; ai<arraySize; ai++) {
			double* pTargetGeom = arrayGeom[ai];
			pTargetGeom[targetVariableIndex] = base + ai*((double)stepInterval/(double)arraySize);
		}
		
		multipleEnergyMinimize(arrayGeom, spec, arraySize, variableSize);
		
		for(ai=0; ai<arraySize; ai++) {
			double x = base + ai*((double)stepInterval/(double)arraySize);
			if(arrayGeom[ai][variableSize+1]  == 1) {
				outFile3 << x << "\t";
				outFile3 << arrayGeom[ai][variableSize] << "\t";
				outFile3 << log(arrayGeom[ai][variableSize]+11.7074);	// loge
			} else {
				outFile3 << x << "\t";
				outFile3 << "-" <<"\t" << "-" << "\t";
			}
			outFile3 << "\n";
		}

		for(ai=0; ai<arraySize; ai++) {
			delete [] arrayGeom[ai];
		}
		delete [] arrayGeom;
		cout << base << " is completed!" << endl; cout.flush();
	}
	outFile3.close( );

	for(m=0; m<24; m++) {
		spec[m]	= true;
	}
	spec[targetVariableIndex] = true;

	ofstream outFile4("dihedralTest\\angleTest_4.txt", ios::out);
	for(base=-180; base<180; base=base+stepInterval) {
		double** arrayGeom = new double* [arraySize];
		for(int ai=0; ai<arraySize; ai++) {
			arrayGeom[ai] = new double[variableSize+2];	// geom + ETOT + flag
		}

		// initializing
		for(ai=0; ai<arraySize; ai++) {
			double* pTargetGeom = arrayGeom[ai];
			for(int vs=0; vs<variableSize; vs++) {
				pTargetGeom[vs] = refGeom[vs];
			}
			pTargetGeom[variableSize]	= 0.0;		// ETOT;
			pTargetGeom[variableSize+1] = 0.0;		// is valid?;
		}

		// change variable
		for(ai=0; ai<arraySize; ai++) {
			double* pTargetGeom = arrayGeom[ai];
			pTargetGeom[targetVariableIndex] = base + ai*((double)stepInterval/(double)arraySize);
		}
		
		multipleEnergyMinimize(arrayGeom, spec, arraySize, variableSize);
		
		for(ai=0; ai<arraySize; ai++) {
			double x = base + ai*((double)stepInterval/(double)arraySize);
			if(arrayGeom[ai][variableSize+1]  == 1) {
				outFile4 << x << "\t";
				outFile4 << arrayGeom[ai][variableSize] << "\t";
				outFile4 << log(arrayGeom[ai][variableSize]+11.7074);	// loge
			} else {
				outFile4 << x << "\t";
				outFile4 << "-" <<"\t" << "-" << "\t";
			}
			outFile4 << "\n";
		}

		for(ai=0; ai<arraySize; ai++) {
			delete [] arrayGeom[ai];
		}
		delete [] arrayGeom;
		cout << base << " is completed!" << endl; cout.flush();
	}
	outFile4.close( );

	delete [] spec;
	delete [] refGeom;
}




void CDifferentialEvolution::test() {
	// Reference conformation
	double Tyr_phi		=	-83.466;
	double Tyr_psi		=	155.791;
	double Tyr_omega	=	-177.129;
	double Tyr_chi1		=	-173.181;
	double Tyr_chi2		=	79.338;
	double Tyr_chi3		=	-166.326;

	double Gly1_phi		=	-154.281;
	double Gly1_psi		=	85.850;
	double Gly1_omega	=	168.513;

	double Gly2_phi		=	82.965;
	double Gly2_psi		=	-75.050;
	double Gly2_omega	=	-169.969;

	double Phe_phi		=	-136.851;
	double Phe_psi		=	19.106;
	double Phe_omega	=	-174.093;
	double Phe_chi1		=	58.858;
	double Phe_chi2		=	-85.471;

	double Met_phi		=	-163.456;
	double Met_psi		=	160.941;
	double Met_omega	=	-179.790;
	double Met_chi1		=	52.868;
	double Met_chi2		=	175.299;
	double Met_chi3		=	-179.862;
	double Met_chi4		=	-58.586;


	double * pGeom = new double[24];
	pGeom[0] =	Tyr_phi;
	pGeom[1] =	Tyr_psi;
	pGeom[2] =	Tyr_omega;
	pGeom[3] =	Tyr_chi1;
	pGeom[4] =	Tyr_chi2;
	pGeom[5] =	Tyr_chi3;

	pGeom[6] =	Gly1_phi;
	pGeom[7] =	Gly1_psi;
	pGeom[8] =	Gly1_omega;

	pGeom[9] =	Gly2_phi;
	pGeom[10] = Gly2_psi;
	pGeom[11] = Gly2_omega;

	pGeom[12] = Phe_phi;
	pGeom[13] = Phe_psi;
	pGeom[14] = Phe_omega;
	pGeom[15] = Phe_chi1;
	pGeom[16] = Phe_chi2;

	pGeom[17] = Met_phi;
	pGeom[18] = Met_psi;
	pGeom[19] = Met_omega;
	pGeom[20] = Met_chi1;
	pGeom[21] = Met_chi2;
	pGeom[22] = Met_chi3;
	pGeom[23] = Met_chi4;

	bool* spec = new bool[24];
	for(int i=0; i<24; i++) spec[i] = false;
	spec[0] = true;

	ofstream outfile1("data\\1.dat", ios::out);
	Tyr_phi	=	-83.453;
	for(int jj=1; jj<50000; jj++) {
		pGeom[0] =  Tyr_phi;
		pGeom[1] =	GARandomFloat(-180, 180);
		pGeom[2] =	GARandomFloat(-180, 180);
		pGeom[3] =	GARandomFloat(-180, 180);
		pGeom[4] =	GARandomFloat(-180, 180);
		pGeom[5] =	GARandomFloat(-180, 180);

		pGeom[6] =	GARandomFloat(-180, 180);
		pGeom[7] =	GARandomFloat(-180, 180);
		pGeom[8] =	GARandomFloat(-180, 180);

		pGeom[9] =	GARandomFloat(-180, 180);
		pGeom[10] = GARandomFloat(-180, 180);
		pGeom[11] = GARandomFloat(-180, 180);

		pGeom[12] = GARandomFloat(-180, 180);
		pGeom[13] = GARandomFloat(-180, 180);
		pGeom[14] = GARandomFloat(-180, 180);
		pGeom[15] = GARandomFloat(-180, 180);
		pGeom[16] = GARandomFloat(-180, 180);

		pGeom[17] = GARandomFloat(-180, 180);
		pGeom[18] = GARandomFloat(-180, 180);
		pGeom[19] = GARandomFloat(-180, 180);
		pGeom[20] = GARandomFloat(-180, 180);
		pGeom[21] = GARandomFloat(-180, 180);
		pGeom[22] = GARandomFloat(-180, 180);
		pGeom[23] = GARandomFloat(-180, 180);

		double ff = singleEnergyMinimize(pGeom, 24);

		outfile1 << jj << "\t" << ff << "\t" << energyEvals << "\t";
		for(int k=0; k<24; k++) {
			outfile1 << pGeom[k] << "\t";
		}
		outfile1 << "\n";

		cout << jj << "\t" << ff << "\n";
		cout.flush();
	}
	outfile1.close();

	
	
	ofstream outfile("data\\1.dat", ios::out);
	for(double iii=-180; iii<=180; iii = iii+0.05) {
		pGeom[0] =  iii;
		pGeom[1] =	Tyr_psi;
		pGeom[2] =	Tyr_omega;
		pGeom[3] =	Tyr_chi1;
		pGeom[4] =	Tyr_chi2;
		pGeom[5] =	Tyr_chi3;

		pGeom[6] =	Gly1_phi;
		pGeom[7] =	Gly1_psi;
		pGeom[8] =	Gly1_omega;

		pGeom[9] =	Gly2_phi;
		pGeom[10] = Gly2_psi;
		pGeom[11] = Gly2_omega;

		pGeom[12] = Phe_phi;
		pGeom[13] = Phe_psi;
		pGeom[14] = Phe_omega;
		pGeom[15] = Phe_chi1;
		pGeom[16] = Phe_chi2;

		pGeom[17] = Met_phi;
		pGeom[18] = Met_psi;
		pGeom[19] = Met_omega;
		pGeom[20] = Met_chi1;
		pGeom[21] = Met_chi2;
		pGeom[22] = Met_chi3;
		pGeom[23] = Met_chi4;


		//double ff = singleEnergyEvaluate(pGeom, 24);
		//double ff = singleEnergyMinimize(pGeom, spec, 24);
		double ff = singleEnergyMinimize(pGeom, 24);


		outfile << iii << "\t" << ff << "\t" << energyEvals << "\t";
		for(int k=0; k<24; k++) {
			outfile << pGeom[k] << "\t";
		}
		outfile << "\n";

		cout << iii << "\t" << ff << "\n";
		cout.flush();
	}




	delete [] spec;
	outfile.close();
}


void CDifferentialEvolution::test2() { 
	ifstream fin("result1//DESummary.txt");
	ofstream outfile("result1//resultSummary.txt", ios::out);

	double score[600]	= {0};
	char charK[20]		= {0};
	char lineBuf[2000]	= {0};
	char wordBuf[100]	= {0};
	for(int k=0; k<50; k++) {
		char filename[50] = "result1//result";
		itoa(k+1, charK, 10);
		strcat(filename, charK);
		strcat(filename, ".txt");
		ofstream outfile1(filename, ios::out);
		for(int i=0; i<600; i++) {
			fin >> wordBuf;
			fin >> wordBuf;
			score[i] += atof(wordBuf);
			outfile1 << atof(wordBuf) << endl;
			fin.getline(lineBuf, 2000);
		}
		outfile1.close();
	}
	
	for(int j=0; j<600; j++) {
		score[j] = score[j]/50;
		outfile << j+1 << "\t"<< score[j] << endl;
		cout << j+1 << "\t"<< score[j]<<endl;
		cout.flush();
	}




	fin.close();
	outfile.close();

	return;
}




void CDifferentialEvolution::randomMinimizeDihedralAnalysis( ) {
	int anglesize		= 24;
	int popsize			= 200;
	int nsize			= 500;

	int i, j, n;

	double** popGeom = new double*[popsize];
	for(i=0; i<popsize; i++) {
		double *iGeom = new double[anglesize];
		for(j=0; j<anglesize;j++) {
			iGeom[j] = 0;
		}
		popGeom[i] = iGeom;
	}
	double *resultValue = new double[popsize];

	ofstream outFile1("dihedral1.txt", ios::out);
	ofstream outFile2("dihedral2.txt", ios::out);

	EnergyCalculator energyCal;
	for(n=0; n<nsize; n++) {
		GARandomSeed( );
		for(i=0; i<popsize; i++) {
			double *iGeom = popGeom[i];
			for(j=0; j<anglesize; j++) {
				iGeom[j] = GARandomFloat(-180, 180);
			}
		}

		energyCal.multipleEnergyEvaluate(popGeom, popsize, anglesize, resultValue);
		for(i=0; i<popsize; i++) {
			double *iGeom = popGeom[i];
			outFile1 << n*popsize+i+1 << "\t";
			outFile1 << resultValue[i] << "\t" ;
			for(j=0; j<anglesize; j++) {
				outFile1 << iGeom[j] << "\t";
			}
			outFile1 << endl;
			outFile1.flush();
		}

		energyCal.multipleEnergyMinimize(popGeom, popsize, anglesize, resultValue);
		for(i=0; i<popsize; i++) {
			double *iGeom = popGeom[i];
			outFile2 << n*popsize+i+1 << "\t";
			outFile2 << resultValue[i] << "\t" ;
			for(j=0; j<anglesize; j++) {
				outFile2 << iGeom[j] << "\t";
			}
			outFile2 << endl;
			outFile2.flush();
		}

		cout << (n+1)*popsize << "  has completed" << endl; cout.flush();
	}

	outFile1.close( );
	outFile2.close( );

	delete [] resultValue;
	for(i=0; i<popsize; i++) {
		double* iGeom = popGeom[i];
		delete [] iGeom;
	}
	delete [] popGeom;
}




double CDifferentialEvolution::calculateRMSD(double **ref_xlist, 
											 double** mov_xlist, 
											 int n_list){
	RMSD rmsdCal;
	return rmsdCal.fast_rmsd(ref_xlist, mov_xlist, n_list);
}


// rmsdType: EnkephalinPDB::AlphaC, EnkephalinPDB::Backbone, EnkephalinPDB::AllAtom;
double CDifferentialEvolution::calculateRMSD(char* refPDBFile, 
											 char* movPDBFile, 
											 EnkephalinPDB::RMSDTYPE rmsdType,
											 int number){
	EnkephalinPDB refPDB;
	EnkephalinPDB movPDB;
	if(!refPDB.readPDB(refPDBFile)) return -1.0;
	if(!movPDB.readPDB(movPDBFile)) return -1.0;

	double ** ref_xlist = refPDB.extractRMSDTypeArray(rmsdType, number);
	double ** mov_xlist = movPDB.extractRMSDTypeArray(rmsdType, number);
	double rmsd = calculateRMSD(ref_xlist, mov_xlist, number);
	refPDB.freeArray(ref_xlist, number);
	movPDB.freeArray(mov_xlist, number);

	return rmsd;
}



void CDifferentialEvolution::randomMinimizeRsmdAnalysis( ) {
	int anglesize		=	  24;			// variables;
	int popsize			=	 100;			// population size;
	int nsize			=	1000;			// iterations;
	
	// native conformation
	double Tyr_phi		=	- 83.466;
	double Tyr_psi		=	 155.791;
	double Tyr_omega	=	-177.129;
	double Tyr_chi1		=	-173.181;
	double Tyr_chi2		=	  79.338;
	double Tyr_chi3		=	-166.326;
	
	double Gly1_phi		=	-154.281;
	double Gly1_psi		=	  85.850;
	double Gly1_omega	=	 168.513;
	
	double Gly2_phi		=	  82.965;
	double Gly2_psi		=	- 75.050;
	double Gly2_omega	=	-169.969;
	
	double Phe_phi		=	-136.851;
	double Phe_psi		=	  19.106;
	double Phe_omega	=	-174.093;
	double Phe_chi1		=	  58.858;
	double Phe_chi2		=	- 85.471;
	
	double Met_phi		=	-163.456;
	double Met_psi		=	 160.941;
	double Met_omega	=	-179.790;
	double Met_chi1		=	  52.868;
	double Met_chi2		=	 175.299;
	double Met_chi3		=	-179.862;
	double Met_chi4		=	- 58.586;
	
	double * refGeom = new double[24];
	refGeom[0]			=	Tyr_phi;
	refGeom[1]			=	Tyr_psi;
	refGeom[2]			=	Tyr_omega;
	refGeom[3]			=	Tyr_chi1;
	refGeom[4]			=	Tyr_chi2;
	refGeom[5]			=	Tyr_chi3;

	refGeom[6]			=	Gly1_phi;
	refGeom[7]			=	Gly1_psi;
	refGeom[8]			=	Gly1_omega;

	refGeom[9]			=	Gly2_phi;
	refGeom[10]			=	Gly2_psi;
	refGeom[11]			=	Gly2_omega;

	refGeom[12]			=	Phe_phi;
	refGeom[13]			=	Phe_psi;
	refGeom[14]			=	Phe_omega;
	refGeom[15]			=	Phe_chi1;
	refGeom[16]			=	Phe_chi2;

	refGeom[17]			=	Met_phi;
	refGeom[18]			=	Met_psi;
	refGeom[19]			=	Met_omega;
	refGeom[20]			=	Met_chi1;
	refGeom[21]			=	Met_chi2;
	refGeom[22]			=	Met_chi3;
	refGeom[23]			=	Met_chi4;

	/***************************************************************************/
	int n, i, j, k, number;
	char mov1PDBFile[100] = {0};
	char mov2PDBFile[100] = {0};
	EnkephalinPDB refPDB;
	EnkephalinPDB::RMSDTYPE rmsdType;
	bool* spec = new bool[anglesize]; for(k=0; k<anglesize; k++) spec[k] = true;
	double* result = new double[popsize];
	/***************************************************************************/

	EnergyCalculator energyCal;
	energyCal.setOutputPDBFormat(true);
	energyCal.setPDBFile("ref");
	energyCal.singleEnergyEvaluate(refGeom, anglesize);
	if(!refPDB.readPDB("ref.pdb")) {cout << "reference pdb file read error!", cout.flush(); return;}

	ofstream outFile1("rmsd1.txt", ios::out);
	ofstream outFile2("rmsd2.txt", ios::out);
	
	for(n = 0; n < nsize; n++) {
		GARandomSeed( );
		double** popGeom = new double*[popsize];
		for(i=0; i<popsize; i++) {
			double *iGeom = new double[anglesize];
			for(j=0; j<anglesize;j++) {
				iGeom[j] = GARandomFloat(-180, 180);
			}
			popGeom[i] = iGeom;
		}
		
		/* setup the optimization variables  "false:fixed", "true: variable"; */
		// spec[ 0] = spec[ 1]	= spec[ 2] = false;
		// spec[ 3] = spec[ 4]	= spec[ 5] = false;
		// spec[ 6] = spec[ 7]	= spec[ 8] = false;
		// spec[ 9] = spec[10]	= spec[11] = false;
		// spec[12] = spec[13]	= spec[14] = false;
		// spec[15] = spec[16]	= spec[17] = false;
		// spec[18] = spec[19]	= spec[20] = false;
		// spec[21] = spec[22]	= spec[23] = false;
		
		RMSD rmsdCal;
		double rmsd1, rmsd2, rmsd3;
		
		/***********************************************************************/
		energyCal.setPDBFile("mov1");
		energyCal.multipleEnergyEvaluate(popGeom, popsize, anglesize, result);
		
		EnkephalinPDB mov1PDB;
		for(k=0; k<popsize; k++) {
			mov1PDB.combineFileName("D:\\MetEnkephalinOptimization\\", "mov1", ".pdb", k+1, mov1PDBFile);
		
			if(mov1PDB.readPDB(mov1PDBFile)) {
				double ** ref_xlist = 0;
				double ** mov_xlist = 0;
		
				rmsdType = EnkephalinPDB::AlphaC; 
				ref_xlist = refPDB.extractRMSDTypeArray(rmsdType, number);
				mov_xlist = mov1PDB.extractRMSDTypeArray(rmsdType, number);
				rmsd1 = rmsdCal.fast_rmsd(ref_xlist, mov_xlist, number);
				mov1PDB.freeArray(mov_xlist, number);
				refPDB.freeArray(ref_xlist, number);
		
				rmsdType = EnkephalinPDB::Backbone; 
				ref_xlist = refPDB.extractRMSDTypeArray(rmsdType, number);
				mov_xlist = mov1PDB.extractRMSDTypeArray(rmsdType, number);
				rmsd2 = rmsdCal.fast_rmsd(ref_xlist, mov_xlist, number);
				mov1PDB.freeArray(mov_xlist, number);
				refPDB.freeArray(ref_xlist, number);
		
				rmsdType = EnkephalinPDB::AllAtom; 
				ref_xlist = refPDB.extractRMSDTypeArray(rmsdType, number);
				mov_xlist = mov1PDB.extractRMSDTypeArray(rmsdType, number);
				rmsd3 = rmsdCal.fast_rmsd(ref_xlist, mov_xlist, number);
				mov1PDB.freeArray(mov_xlist, number);
				refPDB.freeArray(ref_xlist, number);
				
				outFile1 << popsize*n+k+1 << "\t"
						 << setw(8) << setprecision(3)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)  
						 << rmsd1 
						 << setw(8) << setprecision(3)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)  
						 << rmsd2 
						 << setw(8) << setprecision(3)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)  
						 << rmsd3 << "\t";
				outFile1.unsetf(ios::fixed);
				outFile1<< setprecision(8) << mov1PDB.getEnergy() << endl;
				outFile1.flush();
			}
		}
		/***********************************************************************/
		
		
		/***********************************************************************/
		energyCal.setPDBFile("mov2");
		energyCal.multipleEnergyMinimize(popGeom, popsize, anglesize, spec, result);

		for(k=0; k<popsize; k++) {
			EnkephalinPDB mov2PDB;
			mov2PDB.combineFileName("D:\\MetEnkephalinOptimization\\", "mov2", ".pdb", k+1, mov2PDBFile);
			if(mov2PDB.readPDB(mov2PDBFile)) {
				double ** ref_xlist = 0;
				double ** mov_xlist = 0;
				
				rmsdType = EnkephalinPDB::AlphaC; 
				ref_xlist = refPDB.extractRMSDTypeArray(rmsdType, number);
				mov_xlist = mov2PDB.extractRMSDTypeArray(rmsdType, number);
				rmsd1 = rmsdCal.fast_rmsd(ref_xlist, mov_xlist, number);
				mov2PDB.freeArray(mov_xlist, number);
				refPDB.freeArray(ref_xlist, number);

				rmsdType = EnkephalinPDB::Backbone; 
				ref_xlist = refPDB.extractRMSDTypeArray(rmsdType, number);
				mov_xlist = mov2PDB.extractRMSDTypeArray(rmsdType, number);
				rmsd2 = rmsdCal.fast_rmsd(ref_xlist, mov_xlist, number);
				mov2PDB.freeArray(mov_xlist, number);
				refPDB.freeArray(ref_xlist, number);

				rmsdType = EnkephalinPDB::AllAtom; 
				ref_xlist = refPDB.extractRMSDTypeArray(rmsdType, number);
				mov_xlist = mov2PDB.extractRMSDTypeArray(rmsdType, number);
				rmsd3 = rmsdCal.fast_rmsd(ref_xlist, mov_xlist, number);
				mov2PDB.freeArray(mov_xlist, number);
				refPDB.freeArray(ref_xlist, number);

				outFile2<< popsize*n+k+1 << "\t"
						<< setw(8) << setprecision(3)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)  
						<< rmsd1 
						<< setw(8) << setprecision(3)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)  
						<< rmsd2 
						<< setw(8) << setprecision(3)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)  
						<< rmsd3 << "\t";
				outFile2.unsetf(ios::fixed);
				outFile2<< setprecision(8) <<mov2PDB.getEnergy() << endl;
				outFile2.flush();
			}
		}
		/***********************************************************************/

		for(i=0; i<popsize; i++) {
			double *iGeom = popGeom[i];
			delete [] iGeom;
		}
		delete [] popGeom;
		
		
		cout << (n+1)*popsize << " have write to files successfully!" << endl;
	}

	outFile1.close();
	outFile2.close();
	
	delete [] refGeom;
	delete [] spec;
	delete [] result;
}


void CDifferentialEvolution::variableMinimizeSmoothAnalysis( ){
	// native conformation (trial conformation)
	double Tyr_phi		=	-83.466;
	double Tyr_psi		=	155.791;
	double Tyr_omega	=	-177.129;
	double Tyr_chi1		=	-173.181;
	double Tyr_chi2		=	79.338;
	double Tyr_chi3		=	-166.326;

	double Gly1_phi		=	-154.281;
	double Gly1_psi		=	85.850;
	double Gly1_omega	=	168.513;

	double Gly2_phi		=	82.965;
	double Gly2_psi		=	-75.050;
	double Gly2_omega	=	-169.969;

	double Phe_phi		=	-136.851;
	double Phe_psi		=	19.106;
	double Phe_omega	=	-174.093;
	double Phe_chi1		=	58.858;
	double Phe_chi2		=	-85.471;

	double Met_phi		=	-163.456;
	double Met_psi		=	160.941;
	double Met_omega	=	-179.790;
	double Met_chi1		=	52.868;
	double Met_chi2		=	175.299;
	double Met_chi3		=	-179.862;
	double Met_chi4		=	-58.586;

	double * refGeom = new double[24];
	refGeom[0]			=	Tyr_phi;
	refGeom[1]			=	Tyr_psi;
	refGeom[2]			=	Tyr_omega;
	refGeom[3]			=	Tyr_chi1;
	refGeom[4]			=	Tyr_chi2;
	refGeom[5]			=	Tyr_chi3;

	refGeom[6]			=	Gly1_phi;
	refGeom[7]			=	Gly1_psi;
	refGeom[8]			=	Gly1_omega;

	refGeom[9]			=	Gly2_phi;
	refGeom[10]			=	Gly2_psi;
	refGeom[11]			=	Gly2_omega;

	refGeom[12]			=	Phe_phi;
	refGeom[13]			=	Phe_psi;
	refGeom[14]			=	Phe_omega;
	refGeom[15]			=	Phe_chi1;
	refGeom[16]			=	Phe_chi2;

	refGeom[17]			=	Met_phi;
	refGeom[18]			=	Met_psi;
	refGeom[19]			=	Met_omega;
	refGeom[20]			=	Met_chi1;
	refGeom[21]			=	Met_chi2;
	refGeom[22]			=	Met_chi3;
	refGeom[23]			=	Met_chi4;
	
	int i;
	double angle, step;
	bool spec[24] = {true};
	step = 0.1;
	ofstream outFile1("11.txt", ios::out);
	
	EnergyCalculator energyCal;
	for(i=1; i<2; i++) {
		for(angle = -180; angle <= 180; angle=angle+step) {
			refGeom[0]			=	Tyr_phi;
			refGeom[1]			=	Tyr_psi;
			refGeom[2]			=	Tyr_omega;
			refGeom[3]			=	Tyr_chi1;
			refGeom[4]			=	Tyr_chi2;
			refGeom[5]			=	Tyr_chi3;

			refGeom[6]			=	Gly1_phi;
			refGeom[7]			=	Gly1_psi;
			refGeom[8]			=	Gly1_omega;

			refGeom[9]			=	Gly2_phi;
			refGeom[10]			=	Gly2_psi;
			refGeom[11]			=	Gly2_omega;

			refGeom[12]			=	Phe_phi;
			refGeom[13]			=	Phe_psi;
			refGeom[14]			=	Phe_omega;
			refGeom[15]			=	Phe_chi1;
			refGeom[16]			=	Phe_chi2;

			refGeom[17]			=	Met_phi;
			refGeom[18]			=	Met_psi;
			refGeom[19]			=	Met_omega;
			refGeom[20]			=	Met_chi1;
			refGeom[21]			=	Met_chi2;
			refGeom[22]			=	Met_chi3;
			refGeom[23]			=	Met_chi4;
			refGeom[i]			=	angle;
			double value1 = energyCal.singleEnergyEvaluate(refGeom, 24);

			refGeom[0]			=	Tyr_phi;
			refGeom[1]			=	Tyr_psi;
			refGeom[2]			=	Tyr_omega;
			refGeom[3]			=	Tyr_chi1;
			refGeom[4]			=	Tyr_chi2;
			refGeom[5]			=	Tyr_chi3;

			refGeom[6]			=	Gly1_phi;
			refGeom[7]			=	Gly1_psi;
			refGeom[8]			=	Gly1_omega;

			refGeom[9]			=	Gly2_phi;
			refGeom[10]			=	Gly2_psi;
			refGeom[11]			=	Gly2_omega;

			refGeom[12]			=	Phe_phi;
			refGeom[13]			=	Phe_psi;
			refGeom[14]			=	Phe_omega;
			refGeom[15]			=	Phe_chi1;
			refGeom[16]			=	Phe_chi2;

			refGeom[17]			=	Met_phi;
			refGeom[18]			=	Met_psi;
			refGeom[19]			=	Met_omega;
			refGeom[20]			=	Met_chi1;
			refGeom[21]			=	Met_chi2;
			refGeom[22]			=	Met_chi3;
			refGeom[23]			=	Met_chi4;
			refGeom[i]			=	angle;
			for(int j=0; j<24; j++) spec[i] = true;
			spec[i]				=	false;
			double value2 = energyCal.singleEnergyMinimize(refGeom, 24, spec);

			refGeom[0]			=	Tyr_phi;
			refGeom[1]			=	Tyr_psi;
			refGeom[2]			=	Tyr_omega;
			refGeom[3]			=	Tyr_chi1;
			refGeom[4]			=	Tyr_chi2;
			refGeom[5]			=	Tyr_chi3;

			refGeom[6]			=	Gly1_phi;
			refGeom[7]			=	Gly1_psi;
			refGeom[8]			=	Gly1_omega;

			refGeom[9]			=	Gly2_phi;
			refGeom[10]			=	Gly2_psi;
			refGeom[11]			=	Gly2_omega;

			refGeom[12]			=	Phe_phi;
			refGeom[13]			=	Phe_psi;
			refGeom[14]			=	Phe_omega;
			refGeom[15]			=	Phe_chi1;
			refGeom[16]			=	Phe_chi2;

			refGeom[17]			=	Met_phi;
			refGeom[18]			=	Met_psi;
			refGeom[19]			=	Met_omega;
			refGeom[20]			=	Met_chi1;
			refGeom[21]			=	Met_chi2;
			refGeom[22]			=	Met_chi3;
			refGeom[23]			=	Met_chi4;
			refGeom[i]			=	angle;
			for(j=0; j<24; j++) spec[i] = false;
			spec[i]				=	true;
			double value3 = energyCal.singleEnergyMinimize(refGeom, 24, spec);
			
			refGeom[0]			=	Tyr_phi;
			refGeom[1]			=	Tyr_psi;
			refGeom[2]			=	Tyr_omega;
			refGeom[3]			=	Tyr_chi1;
			refGeom[4]			=	Tyr_chi2;
			refGeom[5]			=	Tyr_chi3;

			refGeom[6]			=	Gly1_phi;
			refGeom[7]			=	Gly1_psi;
			refGeom[8]			=	Gly1_omega;

			refGeom[9]			=	Gly2_phi;
			refGeom[10]			=	Gly2_psi;
			refGeom[11]			=	Gly2_omega;

			refGeom[12]			=	Phe_phi;
			refGeom[13]			=	Phe_psi;
			refGeom[14]			=	Phe_omega;
			refGeom[15]			=	Phe_chi1;
			refGeom[16]			=	Phe_chi2;

			refGeom[17]			=	Met_phi;
			refGeom[18]			=	Met_psi;
			refGeom[19]			=	Met_omega;
			refGeom[20]			=	Met_chi1;
			refGeom[21]			=	Met_chi2;
			refGeom[22]			=	Met_chi3;
			refGeom[23]			=	Met_chi4;
			refGeom[i]			=	angle;
			double value4 = energyCal.singleEnergyMinimize(refGeom, 24);

			outFile1 << angle <<"\t"<< value1 << "\t" << value2 << "\t" << value3 << "\t" << value4<< endl; 
			cout << angle <<"\t"<< value1 << "\t" << value2 << "\t" << value3 << "\t" << value4<< endl; 
		}

		outFile1.close( );
	}
	delete [] refGeom;
}

void CDifferentialEvolution::variableMinimizeSmoothAnalysis1( ){
	// native conformation (trial conformation)
	double Tyr_phi		=	-83.466;
	double Tyr_psi		=	155.791;
	double Tyr_omega	=	-177.129;
	double Tyr_chi1		=	-173.181;
	double Tyr_chi2		=	79.338;
	double Tyr_chi3		=	-166.326;
	
	double Gly1_phi		=	-154.281;
	double Gly1_psi		=	85.850;
	double Gly1_omega	=	168.513;
	
	double Gly2_phi		=	82.965;
	double Gly2_psi		=	-75.050;
	double Gly2_omega	=	-169.969;
	
	double Phe_phi		=	-136.851;
	double Phe_psi		=	19.106;
	double Phe_omega	=	-174.093;
	double Phe_chi1		=	58.858;
	double Phe_chi2		=	-85.471;
	
	double Met_phi		=	-163.456;
	double Met_psi		=	160.941;
	double Met_omega	=	-179.790;
	double Met_chi1		=	52.868;
	double Met_chi2		=	175.299;
	double Met_chi3		=	-179.862;
	double Met_chi4		=	-58.586;
	
	double * refGeom = new double[24];
	refGeom[0]			=	Tyr_phi;
	refGeom[1]			=	Tyr_psi;
	refGeom[2]			=	Tyr_omega;
	refGeom[3]			=	Tyr_chi1;
	refGeom[4]			=	Tyr_chi2;
	refGeom[5]			=	Tyr_chi3;
	
	refGeom[6]			=	Gly1_phi;
	refGeom[7]			=	Gly1_psi;
	refGeom[8]			=	Gly1_omega;
	
	refGeom[9]			=	Gly2_phi;
	refGeom[10]			=	Gly2_psi;
	refGeom[11]			=	Gly2_omega;
	
	refGeom[12]			=	Phe_phi;
	refGeom[13]			=	Phe_psi;
	refGeom[14]			=	Phe_omega;
	refGeom[15]			=	Phe_chi1;
	refGeom[16]			=	Phe_chi2;
	
	refGeom[17]			=	Met_phi;
	refGeom[18]			=	Met_psi;
	refGeom[19]			=	Met_omega;
	refGeom[20]			=	Met_chi1;
	refGeom[21]			=	Met_chi2;
	refGeom[22]			=	Met_chi3;
	refGeom[23]			=	Met_chi4;
	
	int i,j, k, angle;
	int step = 30;
	int popsize = 300;
	
	double ** arrayGeom;
	arrayGeom = new double*[popsize];
	for(i=0; i<popsize;i++) {
		double *geom =new double[24];
		arrayGeom[i] = geom;
	}
	double *resultValue1 = new double[popsize];
	double *resultValue2 = new double[popsize];
	double *resultValue3 = new double[popsize];
	double *resultValue4 = new double[popsize];
	bool	spec[24] = {true};
	char	fileName[200] ={"smooth\\native"};
	char	mName[5]   = {0};
		
	EnergyCalculator energyCal;
	for(i=1; i<24; i++) {
		cout << "angle" << i+1 << "\t"; cout.flush();
		strcpy(fileName, "smooth\\native");
		strcat(fileName, itoa(i+1, mName, 10));
		strcat(fileName, ".txt");
			
		ofstream outFile(fileName, ios::out);
		
		for(angle = -180; angle < 180; angle=angle+step) {
			for(j=0; j<popsize;j++) {
				double *tempGeom = arrayGeom[j];
				for(k=0; k<24; k++) {
					tempGeom[k] = refGeom[k];
				}
				tempGeom[i] = angle+ j* (double)step/(double)popsize;
			}
			energyCal.multipleEnergyEvaluate(arrayGeom, popsize, 24, resultValue1);
			cout << "."; cout.flush();
		
		
		
			for(j=0; j<popsize;j++) {
				double *tempGeom = arrayGeom[j];
				for(k=0; k<24; k++) {
					tempGeom[k] = refGeom[k];
				}
				tempGeom[i] = angle+ j* (double)step/(double)popsize;
			}
			for(int m=0; m<24; m++) spec[m] = true; spec[i] = false;
			energyCal.multipleEnergyMinimize(arrayGeom, popsize, 24, spec, resultValue2);
			cout << ".."; cout.flush();


			for(j=0; j<popsize;j++) {
				double *tempGeom = arrayGeom[j];
				for(k=0; k<24; k++) {
					tempGeom[k] = refGeom[k];
				}
				tempGeom[i] = angle+ j* (double)step/(double)popsize;
			}
			for(m=0; m<24; m++) spec[m] = false; spec[i] = true;
			energyCal.multipleEnergyMinimize(arrayGeom, popsize, 24, spec, resultValue3);
			cout << ".."; cout.flush();


			for(j=0; j<popsize;j++) {
				double *tempGeom = arrayGeom[j];
				for(k=0; k<24; k++) {
					tempGeom[k] = refGeom[k];
				}
				tempGeom[i] = angle+ j* (double)step/(double)popsize;
			}
			energyCal.multipleEnergyMinimize(arrayGeom, popsize, 24, resultValue4);
			cout << ".."; cout.flush();


			for(j=0; j<popsize;j++) {
				double *tempGeom = arrayGeom[j];
				for(k=0; k<24; k++) {
					tempGeom[k] = refGeom[k];
				}
				tempGeom[i] = angle+ j* (double)step/(double)popsize;
			}
			for(j=0; j<popsize;j++) {
				double *temp = arrayGeom[j];
				outFile <<  temp[i] <<"\t" << resultValue1[j]<<"\t" << resultValue2[j]<<"\t"<<resultValue3[j]<<"\t"<<resultValue4[j]<< endl;
			}
		}
		cout << endl;

		outFile.close( );
	}

	delete [] refGeom;
	delete [] resultValue1;
	delete [] resultValue2;
	delete [] resultValue3;
	delete [] resultValue4;
	for(i=0; i<popsize;i++) {
		double *geom = arrayGeom[i];
		delete [] geom;
	}
	delete [] arrayGeom;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CAM optimization
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CDifferentialEvolution::CAMObjective(double x0, int index, double* pGEOM, bool *spec, int length, int type) 
{
	double retValue = 1E25;

	pGEOM[index]  = 360*x0-180;
	if(type=1) {
		retValue = singleEnergyEvaluate(pGEOM, length);
	} else {
		retValue = singleEnergyMinimize(pGEOM, spec, length);
	}
	
	return retValue;
}



/* 避免陷入松弛问题局部最优解*/
bool CDifferentialEvolution::illPosed(CL& arg) {
	bool returnValue = false;
	bool small = false;

	int dim = arg.m_dim;
	int count = 0;

	for(int i=0; i<dim; i++) {
		long m = arg.m_l[i]->getId( );
		for(int j=1; j<=dim; j++) {
			if (m==j) {
				count ++;
				break;
			}
		}
	}

	for(int j=0; j<dim; j++) {	
		if (arg.m_l[j]->getElementData(j) < 0.000001) {	//0.000001
			small = true;
			break;
		}
		if(small) break;
	} 

	if( small && count>1) {
		returnValue = true;
	}

	return returnValue;
}

void CDifferentialEvolution::updateLeaf(CLLeafArray& leafArray, Cl& l) {
	long i=0;
	CL *tempCL = 0;
	CLLeafArray tempLeafArray(l.getDim());		

	// 判断条件II,如果不满足条件II则COPY到tempLeafArray,并移除相应的V对应Leaf
	for (i=0; i<leafArray.getCount(); i++) {
		tempCL = leafArray.m_arrayCL[i];
		double *pDialogVector = new double[tempCL->m_dim];
		for(int j=0; j<tempCL->m_dim; j++) {
			pDialogVector[j] = tempCL->getDialogElement(j);
		}
		
		// 判断条件II
		bool condition2 = false;
		for(j=0; j<tempCL->m_dim; j++) {
			if (pDialogVector[j] <= l.getElementData(j)) {
				condition2 = true;
				break;
			}
		}
		delete [] pDialogVector;

		if (!condition2) {
			tempLeafArray.add(*tempCL);
			leafArray.remove(i);
			i--;
		}
	}
	

	// 生成子节点，并判断条件I
	for (i=0; i<tempLeafArray.getCount(); i++) {
		tempCL = tempLeafArray.m_arrayCL[i];
		Cl ** cl = new Cl *[tempCL->m_dim];
		for(int j=0; j<tempCL->m_dim; j++) {
			for(int m=0; m<tempCL->m_dim; m++) {
				cl[m] = tempCL->m_l[m];
			}
			cl[j] =&l;
			CL L(cl,tempCL->m_dim);
			if (L.getCondition1()) {
				if(!illPosed(L)) {
					leafArray.add(L); 
				}
			}
		}
		delete [] cl;
	}
}


double CDifferentialEvolution::CAMOptimize(double *inputGEOM,	// 输入二面角（不改变）
										   double  inputObject, // 输入能量值
										   bool   *spec,		// 优化变量列表
										   int	   index,		// 线性优化变量索引
										   int	   length,
										   double *outPutGEOM) {
	

	double	CONSTANT =100;			// 足够大正数保证Lipschitz连续性	
	int		DIM = 2;				// 问题维数
	double	f   = 0;				// 目标函数值		
	
	int		dim = DIM;
	int		k   = DIM;
	int		K   = 200;


	CL		*TargetL	= 0;
	double	TargetLMin	= 0;
	double	*tempx		= new double [dim];
	double	*templ		= new double [dim];
	int i, j;
	int m;

	for(m=0; m<length;m++) outPutGEOM[m] = inputGEOM[m];
	double bestValue = inputObject;

	ClArray		m_ClArray(dim);			// 支撑向量集合
	CLLeafArray	m_ClLeafArray(dim);		// 树叶集合


	// 初始化I:生成初始支撑向量l1,l2,...,ln
	Cl *m_pCl = 0;
	for (i=0; i<dim; i++) {
		for (int j=0; j<dim; j++) {tempx[j] = 0;}
		tempx[i] = 1;

		double *pTempGEOM = new double[length];	for(m=0; m<length; m++) pTempGEOM[m] = inputGEOM[m];
		f = CAMObjective(tempx[0], index, pTempGEOM, spec, length, 1) + CONSTANT;
		if((f-CONSTANT)<bestValue) {
			bestValue = f-CONSTANT;
			for(m=0; m<length;m++) outPutGEOM[m] = pTempGEOM[m];
		}
		delete [] pTempGEOM;

		for (j=0; j<dim; j++) {templ[j] = tempx[j]/f;}
		m_pCl = new Cl(i+1,templ,dim);
		m_ClArray.add(*m_pCl);
		delete m_pCl;
		m_pCl = 0;
	}

	// 初始化II:生成初始树根节点
	Cl **sl = new Cl* [dim];
	for(i=0; i<dim; i++) {
		sl[i] = m_ClArray.m_l[i];
	}
	CL *m_pCL = new CL(sl, dim);
	m_ClLeafArray.add(*m_pCL);
	delete m_pCL;
	delete [] sl;

	// 迭代进行...
	for (k=DIM; k<=K; k++){
		int leafNumber = m_ClLeafArray.getCount();

		// 计算获取最小树叶
		TargetL	= m_ClLeafArray.findMinimumLeaf();
		// 计算最小树叶对应最低值
		TargetLMin = TargetL->getMinimum( );
		// 计算最小树叶对应变量值
		for(j=0; j<dim; j++) {tempx[j] = TargetL->getDialogElement(j) * TargetLMin;}
		// 计算最小树叶形成支撑向量l^k+1

		double *pTempGEOM = new double[length];	for(m=0; m<length; m++) pTempGEOM[m] = inputGEOM[m];
		f = CAMObjective(tempx[0], index, pTempGEOM, spec, length, 1) + CONSTANT;
		if((f-CONSTANT)<bestValue) {
			bestValue = f-CONSTANT;
			for(m=0; m<length;m++) outPutGEOM[m] = pTempGEOM[m];
		}
		delete [] pTempGEOM;

		if (TargetLMin>=f) break;


		for(j=0; j<dim; j++) {templ[j] = tempx[j]/f;}
		
		m_pCl = new Cl(k+1, templ, dim);
		
		m_ClArray.add(*m_pCl);

		delete m_pCl;

		/*************************************************************************/
		// 控制台输出信息
		cout.setf(ios::fixed|ios::showpoint); 
		cout << setw(6) << k << setw(6) << m_ClLeafArray.getCount();
		for(int j=0; j<dim; j++) {
			cout << setw(8) << TargetL->getClLable(j);
		}
		cout << setw(10) << setprecision(3) <<360*tempx[0]-180;
		cout << setw(10) << setprecision(3) << f-CONSTANT;
		cout << setw(16) << setprecision(6) << TargetLMin-CONSTANT<<"\t";
		cout << endl; 
		/*************************************************************************/
		// 更新树叶
		updateLeaf(m_ClLeafArray, *m_ClArray.m_l[k]);
		if (leafNumber== m_ClLeafArray.getCount()) break;
	}

	delete [] tempx;
	delete [] templ;

	return bestValue;
}








/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Conformational space annealing optimization
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CDifferentialEvolution::dihedralAngleDistance(double i, double j) { 	// -180<=i, j<=180;
	double returnValue = 0.0;

	double dif1 = fabs(i-j);
	double dif2 = fabs(360-dif1);

	if (dif1<=dif2) {
		returnValue = dif1;
	} else {
		returnValue = dif2;
	}

	return returnValue;
}

double CDifferentialEvolution::geometryDistance(double* gi, double* gj, int angleLength) {
	double returnValue = 0.0;

	for(int k=0; k<angleLength; k++) {
		returnValue += dihedralAngleDistance(gi[k], gj[k]);
	}

	return returnValue;
}

double CDifferentialEvolution::geometryDistance(GAGenome& gi, GAGenome& gj){
	double returnValue = 0.0;

	GARealGenome& realGenome1 = dynamic_cast<GARealGenome&>(gi);
	GARealGenome& realGenome2 = dynamic_cast<GARealGenome&>(gj);

	int length = realGenome1.length( );
	for(int k=0; k<length; k++) {
		double i = realGenome1.gene(k);
		double j = realGenome2.gene(k);

		returnValue += dihedralAngleDistance(i, j);
	}

	return returnValue;
}


double CDifferentialEvolution::populationAverageDistance(const GAPopulation& trialPop){
	double returnValue = 0.0;

	int popSize = trialPop.size();
	for(int i=0; i<popSize; i++) {
		for(int j=i+1; j<popSize; j++) {
			GAGenome& gi = trialPop.individual(i);
			GAGenome& gj = trialPop.individual(j);
			returnValue += geometryDistance(gi, gj);
		}
	}

	returnValue = returnValue/(popSize*(popSize-1)/2);
	return returnValue;
}

double CDifferentialEvolution::computeAnnealingX(double aveD, int n, double terminate) {
	double x = 0.0;
	double a = log10((2*terminate)/aveD);
	a = a/(double)n;
	x = pow(10, a);
	return x;
}

double CDifferentialEvolution::distanceCutoff(double aveD, double x, int n){
	double cuttingD = 0.0;
	cuttingD = (aveD/2.0) * pow(x, n);
	if (n>thresholdN) cuttingD = thresholdDistanceCutoff;
	return cuttingD;
}

void CDifferentialEvolution::distanceCuttingCurve(const char* filename, double initialDistanceCutoff, double x) {
	ofstream outFile(filename, ios::out);
	for(int n=1; n<10000; n++) {
		double distCutoff = distanceCutoff(initialDistanceCutoff, x, n);
		outFile << n << "\t";
		outFile << distCutoff << endl;
	}

	outFile.close();
}


SeedStatus CDifferentialEvolution::getGenomeSeedStatus(GAGenome& genome) {
	MyEvalData* pMyEvalData = dynamic_cast<MyEvalData*>(genome.evalData());
	return pMyEvalData->getSeedStatus();
}

void CDifferentialEvolution::setGenomeSeedStatus(GAGenome& genome, SeedStatus status) {
	MyEvalData* pMyEvalData = dynamic_cast<MyEvalData*>(genome.evalData());
	pMyEvalData->setSeedStatus(status);
}

void CDifferentialEvolution::setPopSeedStatus(GAPopulation& triPop, SeedStatus status) {
	for(int i=0; i<triPop.size(); i++) {
		GAGenome& genome = triPop.individual(i);
		setGenomeSeedStatus(genome, status);
	}
}

GAGenome* CDifferentialEvolution::findNearestCorformations(GAGenome& genome, GAPopulation& comparePop) {
	GAGenome* retGenome = 0;
	double distance = -1;

	for(int i=0; i<comparePop.size(); i++) {
		GAGenome& compareGenome = comparePop.individual(i);
		double dist = geometryDistance(genome, compareGenome);
		if(distance == -1 || distance > dist) {
			retGenome = &compareGenome;
			distance = dist;
		}
	}

	return retGenome;
}

GAGenome* CDifferentialEvolution::findNearestCorformations(GAGenome& genome, GAPopulation& comparePop, int index) {
	GAGenome* retGenome = 0;
	double distance = -1;

	for(int i=0; i<index+1; i++) {
		GAGenome& compareGenome = comparePop.individual(i);
		double dist = geometryDistance(genome, compareGenome);
		if(distance == -1 || distance > dist) {
			retGenome = &compareGenome;
			distance = dist;
		}
	}

	return retGenome;
}


GAGenome* CDifferentialEvolution::getSeedConformation(GAPopulation& bank) {
	GAGenome* seed = 0;
	bank.sort(gaTrue);
	int bankSize = bank.size();
	double distAverage = populationAverageDistance(bank);


	if(preSeed == 0) {		
		GAGenome& seedGenome = bank.individual(0);	// the lowest energy conformation in the current bank;
		setGenomeSeedStatus(seedGenome, DETrue);
		preSeed = seed = &seedGenome;
	} else {	// subsequent seed conformation is selected in a different manner from firstly selected seed conformation
		for(int i=0; i<bankSize; i++) {
			GAGenome& genome = bank.individual(i);
			if(getGenomeSeedStatus(genome)==DEUnknown) {		// condition1: it has not yet been used as a seed conformation
				double dist = geometryDistance(genome, *preSeed);
				if(dist >= distAverage) {	// condition2: it lies in the conformation space sufficiently far away from the previous seed conformation (ave distance)
					setGenomeSeedStatus(genome,DETrue);
					seed = &genome;
					preSeed= seed;
					break;
				}
				else {
					//setGenomeSeedStatus(genome, DEFalse);
					continue;
				}
			}
		}
	}

	return seed;
}

bool CDifferentialEvolution::hasSeedConformation(GAPopulation& bank) {
	bool bSeed = false;

	bank.sort(gaTrue);
	int bankSize = bank.size();
	double distAverage = populationAverageDistance(bank);


	if(preSeed == 0) {
		bSeed = true;
	} else {	// subsequent seed conformation is selected in a different manner from firstly selected seed conformation
		for(int i=0; i<bankSize; i++) {
			GAGenome& genome = bank.individual(i);
			if(getGenomeSeedStatus(genome)==DEUnknown) {		// condition1: it has not yet been used as a seed conformation
				double dist = geometryDistance(genome, *preSeed);
				if(dist >= distAverage) {	// condition2: it lies in the conformation space sufficiently far away from the previous seed conformation (ave distance)
					bSeed = true;
					break;
				}
				else {
					//setGenomeSeedStatus(genome, DEFalse);
					continue;
				}
			}
		}
	}

	return bSeed;
}


GAPopulation* CDifferentialEvolution::getTrialConformations(GAGenome& seedGenome, 
															const GAPopulation& firstBank, 
															const GAPopulation& bank) {
	GARealGenome& realSeedGenome = dynamic_cast<GARealGenome&>(seedGenome);
	GAPopulation* trialPop = new GAPopulation(realSeedGenome, 10);
	trialPop->order(GAPopulation::LOW_IS_BEST);
	
	int tryNumber = 100;


	for(int i=0; i<trialPop->size(); i++) {
		GARealGenome& genome = dynamic_cast<GARealGenome&>(trialPop->individual(i));
		genome.copy(realSeedGenome);
		setGenomeSeedStatus(genome, DEUnknown);
	}

	// Step1: three conformations are generated by changing a few(one to four) randomly selected dihedral angles of the seed conformations
	for(i=0; i<3; i++) {
		int dihedralAngleVaried = GARandomInt(1,4);
		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(trialPop->individual(i));
		int length = realGenome.length( );
		for(int j=0; j<dihedralAngleVaried; j++) {
			GARealGenome& randomGenome = dynamic_cast<GARealGenome&>(firstBank.individual(GARandomInt(0, firstBank.size()-1)));
			int site = GARandomInt(0, length-1);
			realGenome.gene(site, randomGenome.gene(site));
		}
	
		if(i>0 && tryNumber>0) {	// assure that all 10 trial conformations are significantly different (at least 30)
			GAGenome* nearestGenome = findNearestCorformations(realGenome, *trialPop, i-1);
			if(geometryDistance(realGenome, *nearestGenome) < thresholdDifference ) {
				i--;
				tryNumber--;
			}
		}
	}

	// Step2: three conformations are generated as above but the dihedral angles are selected only from phi, psi, chi1, that is, the dihedral angles around the alpha carbon atom.
	int indexArray[] = {0, 1, 3, 6, 7, 9, 10, 12, 13, 15, 17, 18, 20};
	for(i=3; i<6; i++) {
		int dihedralAngleVaried = GARandomInt(1,4);
		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(trialPop->individual(i));
		int length = realGenome.length( );
		for(int j=0; j<dihedralAngleVaried; j++) {
			GARealGenome& randomGenome = dynamic_cast<GARealGenome&>(firstBank.individual(GARandomInt(0, firstBank.size()-1)));
			int site = GARandomInt(0, 12);
			realGenome.gene(indexArray[site], randomGenome.gene(indexArray[site]));
		}

		// assure that all 10 trial conformations are significantly different (at least 30)
		if(tryNumber>0) {
			GAGenome* nearestGenome = findNearestCorformations(realGenome, *trialPop, i-1);
			if(geometryDistance(realGenome, *nearestGenome) < thresholdDifference ) {
				i--;
				tryNumber--;
			}
		}
	}

	// Step3: two conformations are generated by randomly replacing one of the eight groups of dihedral angles of the seed conformation with a corresponding group of dihedral angles of a randomly selected in the bank
	for(i = 6; i<8; i++) {
		int group = GARandomInt(1,8);

		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(trialPop->individual(i));
		GARealGenome& randGenome = dynamic_cast<GARealGenome&>(bank.individual(GARandomInt(0, bank.size()-1)));
		
		if(group==1) {
			realGenome.gene(0, randGenome.gene(0));
			realGenome.gene(1, randGenome.gene(1));
			realGenome.gene(2, randGenome.gene(2));
		} else if(group == 2) {
			realGenome.gene(3, randGenome.gene(3));
			realGenome.gene(4, randGenome.gene(4));
			realGenome.gene(5, randGenome.gene(5));
		} else if(group == 3) {
			realGenome.gene(6, randGenome.gene(6));
			realGenome.gene(7, randGenome.gene(7));
			realGenome.gene(8, randGenome.gene(8));
		} else if(group == 4) {
			realGenome.gene( 9, randGenome.gene( 9));
			realGenome.gene(10, randGenome.gene(10));
			realGenome.gene(11, randGenome.gene(11));
		} else if(group == 5) {
			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));
		} else if(group == 6) {
			realGenome.gene(15, randGenome.gene(15));
			realGenome.gene(16, randGenome.gene(16));
		} else if(group == 7) {
			realGenome.gene(17, randGenome.gene(17));
			realGenome.gene(18, randGenome.gene(18));
			realGenome.gene(19, randGenome.gene(19));
		} else {
			realGenome.gene(20, randGenome.gene(20));
			realGenome.gene(21, randGenome.gene(21));
			realGenome.gene(22, randGenome.gene(22));
			realGenome.gene(23, randGenome.gene(23));
		}

		// assure that all 10 trial conformations are significantly different (at least 30)
		if(tryNumber>0) {
			GAGenome* nearestGenome = findNearestCorformations(realGenome, *trialPop, i-1);
			if(geometryDistance(realGenome, *nearestGenome) < thresholdDifference ) {
				i--;
				tryNumber--;
			}
		}
	}

	// Step4: two conformations are generated as above but using connected groups instead of groups. 
	for(i=8; i<10; i++) {
		int connectGroup = GARandomInt(1, 7);
		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(trialPop->individual(i));
		GARealGenome& randGenome = dynamic_cast<GARealGenome&>(bank.individual(GARandomInt(0, bank.size()-1)));
		if (connectGroup ==1) {			// 1, 2, 3
			realGenome.gene(0, randGenome.gene(0));
			realGenome.gene(1, randGenome.gene(1));
			realGenome.gene(2, randGenome.gene(2));

			realGenome.gene(3, randGenome.gene(3));
			realGenome.gene(4, randGenome.gene(4));
			realGenome.gene(5, randGenome.gene(5));
		
			realGenome.gene(6, randGenome.gene(6));
			realGenome.gene(7, randGenome.gene(7));
			realGenome.gene(8, randGenome.gene(8));
		} else if(connectGroup == 2) {	// 1, 3, 4
			realGenome.gene(0, randGenome.gene(0));
			realGenome.gene(1, randGenome.gene(1));
			realGenome.gene(2, randGenome.gene(2));
		
			realGenome.gene(6, randGenome.gene(6));
			realGenome.gene(7, randGenome.gene(7));
			realGenome.gene(8, randGenome.gene(8));
		
			realGenome.gene(9, randGenome.gene(9));
			realGenome.gene(10, randGenome.gene(10));
			realGenome.gene(11, randGenome.gene(11));
		} else if(connectGroup == 3) {	// 3, 4, 5
			realGenome.gene(6, randGenome.gene(6));
			realGenome.gene(7, randGenome.gene(7));
			realGenome.gene(8, randGenome.gene(8));
		
			realGenome.gene(9, randGenome.gene(9));
			realGenome.gene(10, randGenome.gene(10));
			realGenome.gene(11, randGenome.gene(11));

			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));
		} else if(connectGroup == 4) {	//4, 5, 6
			realGenome.gene(9,  randGenome.gene(9));
			realGenome.gene(10, randGenome.gene(10));
			realGenome.gene(11, randGenome.gene(11));

			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));

			realGenome.gene(15, randGenome.gene(15));
			realGenome.gene(16, randGenome.gene(16));
		} else if(connectGroup == 5) {	// 4, 5, 7
			realGenome.gene(9,  randGenome.gene(9));
			realGenome.gene(10, randGenome.gene(10));
			realGenome.gene(11, randGenome.gene(11));

			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));

			realGenome.gene(17, randGenome.gene(17));
			realGenome.gene(18, randGenome.gene(18));
			realGenome.gene(19, randGenome.gene(19));

		} else if(connectGroup == 6) {	// 5, 6, 7
			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));

			realGenome.gene(15, randGenome.gene(15));
			realGenome.gene(16, randGenome.gene(16));

			realGenome.gene(17, randGenome.gene(17));
			realGenome.gene(18, randGenome.gene(18));
			realGenome.gene(19, randGenome.gene(19));
		} else {						// 5, 7, 8
			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));

			realGenome.gene(17, randGenome.gene(17));
			realGenome.gene(18, randGenome.gene(18));
			realGenome.gene(19, randGenome.gene(19));

			realGenome.gene(20, randGenome.gene(20));
			realGenome.gene(21, randGenome.gene(21));
			realGenome.gene(22, randGenome.gene(22));
			realGenome.gene(23, randGenome.gene(23));
		}

		// assure that all 10 trial conformations are significantly different (at least 30)
		if(tryNumber>0) {
			GAGenome* nearestGenome = findNearestCorformations(realGenome, *trialPop, i-1);
			if(geometryDistance(realGenome, *nearestGenome) < thresholdDifference ) {
				i--;
				tryNumber--;
			}
		}
	}
	multipleEnergyMinimize(*trialPop);
	numberMinimized += trialPop->size();
	currentDistanceCutoff = distanceCutoff(averageDistanceOfFirstBank, xCutoff, numberMinimized);

	return trialPop;
}

GAPopulation* CDifferentialEvolution::getTrialConformationsWithoutMinimized(GAGenome& seedGenome, 
															const GAPopulation& firstBank, 
															const GAPopulation& bank) {
	GARealGenome& realSeedGenome = dynamic_cast<GARealGenome&>(seedGenome);
	GAPopulation* trialPop = new GAPopulation(realSeedGenome, 10);
	trialPop->order(GAPopulation::LOW_IS_BEST);
	
	int tryNumber = 100;


	for(int i=0; i<trialPop->size(); i++) {
		GARealGenome& genome = dynamic_cast<GARealGenome&>(trialPop->individual(i));
		genome.copy(realSeedGenome);
		setGenomeSeedStatus(genome, DEUnknown);
	}

	// Step1: three conformations are generated by changing a few(one to four) randomly selected dihedral angles of the seed conformations
	for(i=0; i<3; i++) {
		int dihedralAngleVaried = GARandomInt(1,4);
		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(trialPop->individual(i));
		int length = realGenome.length( );
		for(int j=0; j<dihedralAngleVaried; j++) {
			GARealGenome& randomGenome = dynamic_cast<GARealGenome&>(firstBank.individual(GARandomInt(0, firstBank.size()-1)));
			int site = GARandomInt(0, length-1);
			realGenome.gene(site, randomGenome.gene(site));
		}
	
		if(i>0 && tryNumber>0) {	// assure that all 10 trial conformations are significantly different (at least 30)
			GAGenome* nearestGenome = findNearestCorformations(realGenome, *trialPop, i-1);
			if(geometryDistance(realGenome, *nearestGenome) < thresholdDifference ) {
				i--;
				tryNumber--;
			}
		}
	}

	// Step2: three conformations are generated as above but the dihedral angles are selected only from phi, psi, chi1, that is, the dihedral angles around the alpha carbon atom.
	int indexArray[] = {0, 1, 3, 6, 7, 9, 10, 12, 13, 15, 17, 18, 20};
	for(i=3; i<6; i++) {
		int dihedralAngleVaried = GARandomInt(1,4);
		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(trialPop->individual(i));
		int length = realGenome.length( );
		for(int j=0; j<dihedralAngleVaried; j++) {
			GARealGenome& randomGenome = dynamic_cast<GARealGenome&>(firstBank.individual(GARandomInt(0, firstBank.size()-1)));
			int site = GARandomInt(0, 12);
			realGenome.gene(indexArray[site], randomGenome.gene(indexArray[site]));
		}

		// assure that all 10 trial conformations are significantly different (at least 30)
		if(tryNumber>0) {
			GAGenome* nearestGenome = findNearestCorformations(realGenome, *trialPop, i-1);
			if(geometryDistance(realGenome, *nearestGenome) < thresholdDifference ) {
				i--;
				tryNumber--;
			}
		}
	}

	// Step3: two conformations are generated by randomly replacing one of the eight groups of dihedral angles of the seed conformation with a corresponding group of dihedral angles of a randomly selected in the bank
	for(i = 6; i<8; i++) {
		int group = GARandomInt(1,8);

		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(trialPop->individual(i));
		GARealGenome& randGenome = dynamic_cast<GARealGenome&>(bank.individual(GARandomInt(0, bank.size()-1)));
		
		if(group==1) {
			realGenome.gene(0, randGenome.gene(0));
			realGenome.gene(1, randGenome.gene(1));
			realGenome.gene(2, randGenome.gene(2));
		} else if(group == 2) {
			realGenome.gene(3, randGenome.gene(3));
			realGenome.gene(4, randGenome.gene(4));
			realGenome.gene(5, randGenome.gene(5));
		} else if(group == 3) {
			realGenome.gene(6, randGenome.gene(6));
			realGenome.gene(7, randGenome.gene(7));
			realGenome.gene(8, randGenome.gene(8));
		} else if(group == 4) {
			realGenome.gene( 9, randGenome.gene( 9));
			realGenome.gene(10, randGenome.gene(10));
			realGenome.gene(11, randGenome.gene(11));
		} else if(group == 5) {
			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));
		} else if(group == 6) {
			realGenome.gene(15, randGenome.gene(15));
			realGenome.gene(16, randGenome.gene(16));
		} else if(group == 7) {
			realGenome.gene(17, randGenome.gene(17));
			realGenome.gene(18, randGenome.gene(18));
			realGenome.gene(19, randGenome.gene(19));
		} else {
			realGenome.gene(20, randGenome.gene(20));
			realGenome.gene(21, randGenome.gene(21));
			realGenome.gene(22, randGenome.gene(22));
			realGenome.gene(23, randGenome.gene(23));
		}

		// assure that all 10 trial conformations are significantly different (at least 30)
		if(tryNumber>0) {
			GAGenome* nearestGenome = findNearestCorformations(realGenome, *trialPop, i-1);
			if(geometryDistance(realGenome, *nearestGenome) < thresholdDifference ) {
				i--;
				tryNumber--;
			}
		}
	}

	// Step4: two conformations are generated as above but using connected groups instead of groups. 
	for(i=8; i<10; i++) {
		int connectGroup = GARandomInt(1, 7);
		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(trialPop->individual(i));
		GARealGenome& randGenome = dynamic_cast<GARealGenome&>(bank.individual(GARandomInt(0, bank.size()-1)));
		if (connectGroup ==1) {			// 1, 2, 3
			realGenome.gene(0, randGenome.gene(0));
			realGenome.gene(1, randGenome.gene(1));
			realGenome.gene(2, randGenome.gene(2));

			realGenome.gene(3, randGenome.gene(3));
			realGenome.gene(4, randGenome.gene(4));
			realGenome.gene(5, randGenome.gene(5));
		
			realGenome.gene(6, randGenome.gene(6));
			realGenome.gene(7, randGenome.gene(7));
			realGenome.gene(8, randGenome.gene(8));
		} else if(connectGroup == 2) {	// 1, 3, 4
			realGenome.gene(0, randGenome.gene(0));
			realGenome.gene(1, randGenome.gene(1));
			realGenome.gene(2, randGenome.gene(2));
		
			realGenome.gene(6, randGenome.gene(6));
			realGenome.gene(7, randGenome.gene(7));
			realGenome.gene(8, randGenome.gene(8));
		
			realGenome.gene(9, randGenome.gene(9));
			realGenome.gene(10, randGenome.gene(10));
			realGenome.gene(11, randGenome.gene(11));
		} else if(connectGroup == 3) {	// 3, 4, 5
			realGenome.gene(6, randGenome.gene(6));
			realGenome.gene(7, randGenome.gene(7));
			realGenome.gene(8, randGenome.gene(8));
		
			realGenome.gene(9, randGenome.gene(9));
			realGenome.gene(10, randGenome.gene(10));
			realGenome.gene(11, randGenome.gene(11));

			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));
		} else if(connectGroup == 4) {	//4, 5, 6
			realGenome.gene(9,  randGenome.gene(9));
			realGenome.gene(10, randGenome.gene(10));
			realGenome.gene(11, randGenome.gene(11));

			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));

			realGenome.gene(15, randGenome.gene(15));
			realGenome.gene(16, randGenome.gene(16));
		} else if(connectGroup == 5) {	// 4, 5, 7
			realGenome.gene(9,  randGenome.gene(9));
			realGenome.gene(10, randGenome.gene(10));
			realGenome.gene(11, randGenome.gene(11));

			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));

			realGenome.gene(17, randGenome.gene(17));
			realGenome.gene(18, randGenome.gene(18));
			realGenome.gene(19, randGenome.gene(19));

		} else if(connectGroup == 6) {	// 5, 6, 7
			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));

			realGenome.gene(15, randGenome.gene(15));
			realGenome.gene(16, randGenome.gene(16));

			realGenome.gene(17, randGenome.gene(17));
			realGenome.gene(18, randGenome.gene(18));
			realGenome.gene(19, randGenome.gene(19));
		} else {						// 5, 7, 8
			realGenome.gene(12, randGenome.gene(12));
			realGenome.gene(13, randGenome.gene(13));
			realGenome.gene(14, randGenome.gene(14));

			realGenome.gene(17, randGenome.gene(17));
			realGenome.gene(18, randGenome.gene(18));
			realGenome.gene(19, randGenome.gene(19));

			realGenome.gene(20, randGenome.gene(20));
			realGenome.gene(21, randGenome.gene(21));
			realGenome.gene(22, randGenome.gene(22));
			realGenome.gene(23, randGenome.gene(23));
		}

		// assure that all 10 trial conformations are significantly different (at least 30)
		if(tryNumber>0) {
			GAGenome* nearestGenome = findNearestCorformations(realGenome, *trialPop, i-1);
			if(geometryDistance(realGenome, *nearestGenome) < thresholdDifference ) {
				i--;
				tryNumber--;
			}
		}
	}

	return trialPop;
}

void CDifferentialEvolution::updateBank(GAPopulation& trialPop, GAPopulation& bank) {
	// use the global minimize number variable to compute the current distance cutoff (niching)
	for(int i=0; i<trialPop.size( ); i++) {
		GARealGenome& genome = dynamic_cast<GARealGenome&>(trialPop.individual(i));
		GAGenome* nearestGenome = findNearestCorformations(genome, bank);
		double minDistance = geometryDistance(genome, *nearestGenome);
		if(minDistance <= currentDistanceCutoff) {
			if (genome.score() < nearestGenome->score()) {	// updating group, otherwise discard it
				nearestGenome->copy(genome);
			}
		}else {
			bank.sort(gaTrue);
			GAGenome& worstGenome = bank.individual(bank.size()-1);
			if(genome.score( ) < worstGenome.score( )) {
				worstGenome.copy(genome);
			}
		}
	}
}


void CDifferentialEvolution::updateBank2(GAPopulation& trialPop, GAPopulation& bank) {
	for(int i=0; i<trialPop.size( ); i++) {
		double curDistanceCutoff = distanceCutoff(averageDistanceOfFirstBank, xCutoff, numberN++);
		GARealGenome& genome = dynamic_cast<GARealGenome&>(trialPop.individual(i));
		GAGenome* nearestGenome = findNearestCorformations(genome, bank);
		double minDistance = geometryDistance(genome, *nearestGenome);
		if(minDistance <= curDistanceCutoff) {
			if (genome.score() < nearestGenome->score()) {	// updating group, otherwise discard it
				nearestGenome->copy(genome);
			}
		}else {
			bank.sort(gaTrue);
			GAGenome& worstGenome = bank.individual(bank.size()-1);
			if(genome.score( ) < worstGenome.score( )) {
				worstGenome.copy(genome);
			}
		}
	}
}


void CDifferentialEvolution::conformationSpaceAnnealingContinuousOptimize( ) {
	for(int k=0; k<50; k++) {
		conformationSpaceAnnealingOptimize( );
	}
}

void CDifferentialEvolution::conformationSpaceAnnealingOptimize( ) {
	////////////////////////////////////////////////////////////////////////////////
	ofstream outfile("CSASummary.txt", ios::app);
	////////////////////////////////////////////////////////////////////////////////

	// parameters setup
	thresholdN = 5000;
	thresholdDistanceCutoff = 90;
	thresholdDifference = 30;

	// variables initialize
	energyEvals			= 0.0;
	energyEvalsCPUTime	= 0.0;
	numberMinimized		= 0.0;
	numberN             = 0.0;
	preSeed             = 0;

	// "FirstBank" initialize
	GARandomSeed(0);
	pop->initialize( );
	multipleEnergyMinimize(*pop);
	numberMinimized += pop->size();
	numberN += pop->size();

	averageDistanceOfFirstBank = populationAverageDistance(*pop);
	xCutoff = computeAnnealingX(averageDistanceOfFirstBank, thresholdN, thresholdDistanceCutoff);
	currentDistanceCutoff = distanceCutoff(averageDistanceOfFirstBank, xCutoff, numberMinimized);

	// create the "Bank" and copy the "FirstBank" to it
	GAPopulation *bankPopulation = pop->clone(); 
	
	bool flag = false;
	for(int k=0; k<9; k++) {
		preSeed = 0;
		setPopSeedStatus(*bankPopulation, DEUnknown);

		GAGenome* seedGenome;
		while(seedGenome = getSeedConformation(*bankPopulation)) {
			GAPopulation * pTrialConformations = getTrialConformations(*seedGenome, *pop, *bankPopulation);
			updateBank(*pTrialConformations, *bankPopulation);
		
			GARealGenome& bestGenome = dynamic_cast<GARealGenome&>(bankPopulation->best());
			cout<< setw(3) 
				<< k+1;
			cout<< setw(10) 
				<< setprecision(4) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< bestGenome.score( );
		
			cout<< setw(12) 
				<< setprecision(0) 
				<< setiosflags(ios::fixed)
				<< numberMinimized;

			cout<< setw(12) 
				<< setprecision(0) 
				<< setiosflags(ios::fixed)
				<< energyEvals;
			cout<< setw(12) 
				<< setprecision(3) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< energyEvalsCPUTime;

			cout<< setw(12) 
				<< setprecision(0) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< currentDistanceCutoff;

			cout<< setw(12) 
				<< setprecision(0) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< populationAverageDistance(*bankPopulation);

			cout<< endl;
			cout.flush();

			delete pTrialConformations;
			
			if(bestGenome.score() <-11.707) {flag = true; break;}
		}
		if(flag) break;
	}

	////////////////////////////////////////////////////////////////////////////////
	GARealGenome& bestGenome = dynamic_cast<GARealGenome&>(bankPopulation->best());
	if(flag) {
		outfile << 1 << "\t";
		outfile << bestGenome.score() << "\t";
		outfile << numberMinimized<< "\t";
		outfile << energyEvals<< "\t";
		outfile << energyEvalsCPUTime<< "\t";
		outfile << currentDistanceCutoff<< "\t";
		outfile << populationAverageDistance(*bankPopulation)<< "\t";
		for(int ii=0; ii<bestGenome.length(); ii++) {
			outfile<< bestGenome.gene(ii) << "\t";
		}
		outfile << endl;
	} else {
		outfile << 0 << "\t";
		outfile << bestGenome.score() << "\t";
		outfile << numberMinimized<< "\t";
		outfile << energyEvals<< "\t";
		outfile << energyEvalsCPUTime<< "\t";
		outfile << currentDistanceCutoff<< "\t";
		outfile << populationAverageDistance(*bankPopulation)<< "\t";
		for(int ii=0; ii<bestGenome.length(); ii++) {
			outfile<< bestGenome.gene(ii) << "\t";
		}
		outfile << endl;
	}
	outfile.close( );
	////////////////////////////////////////////////////////////////////////////////

	delete bankPopulation;
}


void CDifferentialEvolution::conformationSpaceAnnealingOptimize2( ) {
	// parameters setup
	thresholdN = 5000;
	thresholdDistanceCutoff = 90;
	thresholdDifference = 30;

	// variables initialize
	energyEvals			= 0.0;
	energyEvalsCPUTime	= 0.0;
	numberMinimized		= 0.0;
	numberN             = 0.0;
	preSeed             = 0;

	// "FirstBank" initialize
	GARandomSeed(0);
	pop->initialize( );
	multipleEnergyMinimize(*pop);
	numberMinimized += pop->size();
	numberN += pop->size();

	averageDistanceOfFirstBank = populationAverageDistance(*pop);
	xCutoff = computeAnnealingX(averageDistanceOfFirstBank, thresholdN, thresholdDistanceCutoff);
	currentDistanceCutoff = distanceCutoff(averageDistanceOfFirstBank, xCutoff, numberMinimized);

	// create the "Bank" and copy the "FirstBank" to it
	GAPopulation *bankPopulation = pop->clone(); 
	
	for(int k=0; k<9; k++) {
		preSeed = 0;
		setPopSeedStatus(*bankPopulation, DEUnknown);

		GAGenome* seedGenome;
		while(hasSeedConformation(*bankPopulation)) {
			GAPopulation* pAllTrials = 0;

			while(seedGenome = getSeedConformation(*bankPopulation)) {
				GAPopulation * pTrialConformations = getTrialConformationsWithoutMinimized(*seedGenome, *pop, *bankPopulation);
				if(pAllTrials == 0) {
					pAllTrials = new GAPopulation(*pTrialConformations);
				} else {
					for(int kk=0; kk<pTrialConformations->size(); kk++) {
						GAGenome& pTrialGenome = pTrialConformations->individual(kk);
						pAllTrials->add(pTrialGenome);
					}
				}
				delete pTrialConformations;
			}
			
			multipleEnergyMinimize(*pAllTrials);
			numberMinimized += pAllTrials->size();

			updateBank2(*pAllTrials, *bankPopulation);
			currentDistanceCutoff = distanceCutoff(averageDistanceOfFirstBank, xCutoff, numberMinimized);

			//////////////////////////////////////////////////////////////////////////////
			GARealGenome& bestGenome = dynamic_cast<GARealGenome&>(bankPopulation->best());
			cout<< setw(3) 
				<< k+1;
			cout<< setw(10) 
				<< setprecision(4) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< bestGenome.score( );
		
			cout<< setw(12) 
				<< setprecision(0) 
				<< setiosflags(ios::fixed)
				<< numberMinimized;

			cout<< setw(12) 
				<< setprecision(0) 
				<< setiosflags(ios::fixed)
				<< energyEvals;
			cout<< setw(12) 
				<< setprecision(3) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< energyEvalsCPUTime;

			cout<< setw(12) 
				<< setprecision(0) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< currentDistanceCutoff;

			cout<< setw(12) 
				<< setprecision(0) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< populationAverageDistance(*bankPopulation);

			cout<< endl;
			cout.flush();
			//////////////////////////////////////////////////////////////////////////////
			

			delete pAllTrials;
		}
	}

	delete bankPopulation;
}
