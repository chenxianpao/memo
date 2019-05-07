// MetEnkephalinOptimization.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <stdio.h>
#include <iostream.h>
#include <ga/ga.h>
#include <ga/GARealGenome.h>
#include <ga/GARealGenome.C>
#include <delb.h>
#include <differentialevolution.h>
#include <SAAlgorithm.h>
#include <stdio.h>
#include <io.h> 
#include <iomanip.h>
#include <cam.h>
#include <math.h>
#include <dataStruct.h>
#include <MyEvalData.h>

#include  "MOCrowdingDEvolution.h"







int number = 0;	 //计算目标函数评价次数；

float genomeEvaluate(GAGenome& genome) {
	number++;	 // 评价次数

	double result = 0;
	
	GARealGenome& g = dynamic_cast<GARealGenome&>(genome);
	FunctionStruct *pFS = (FunctionStruct*)(g.userData());
	ObjectiveFunction objectiveFunction = pFS->objectiveFunction;
	
	int length = g.length();
	double *x = new double[length];
	for(int i=0; i<length; i++) x[i] = g.gene(i);
	result = objectiveFunction(length, x);
	delete[] x;
	
	return result;
}

double metEnkePhalinSingleEnergyEvaluation(int dim, double *geomData) {
	const char* filename = "metEnkephalinEnergyEvaluationSourceFile.inp";
	const char* readfile = "main_out.MET_ENKEPHALIN_ENERGY_EVALUATION";
	const char* logfile  = "log";

	double ETOT = 1E20;

	// write the configuration file used for ECEPP/3 model engine
	ofstream outFile(filename, ios::out);
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
		outFile	<< setw(8) // 设定field宽度
				<< setprecision(3) // 设置小数位置
				<< setiosflags(ios::showpoint) // keep trailing 0s
				<< setiosflags(ios::fixed) // 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=6; i<9;i++) {
		outFile	<< setw(8) // 设定field宽度
				<< setprecision(3) // 设置小数位置
				<< setiosflags(ios::showpoint) // keep trailing 0s
				<< setiosflags(ios::fixed) // 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=9; i<12;i++) {
		outFile	<< setw(8) // 设定field宽度
				<< setprecision(3) // 设置小数位置
				<< setiosflags(ios::showpoint) // keep trailing 0s
				<< setiosflags(ios::fixed) // 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=12; i<17;i++) {
		outFile	<< setw(8) // 设定field宽度
				<< setprecision(3) // 设置小数位置
				<< setiosflags(ios::showpoint) // keep trailing 0s
				<< setiosflags(ios::fixed) // 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=17; i<24;i++) {
		outFile	<< setw(8) // 设定field宽度
				<< setprecision(3) // 设置小数位置
				<< setiosflags(ios::showpoint) // keep trailing 0s
				<< setiosflags(ios::fixed) // 使用这些设置
				<< geomData[i];
	}
	outFile << endl;
	outFile << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();


	// run ECEPP/3 model engine
	system("i.bat ENERGY metEnkephalinEnergyEvaluationSourceFile MET_ENKEPHALIN_ENERGY_EVALUATION x x  1>> log 2> err");
	
	// read ETOT from main_out.MET_ENKEPHALIN_ENERGY_EVALUATION file
	ifstream fin(readfile); 
	char buf[1000]={0};
	bool find = false;
	do{
		fin	>>  buf;
		if(find) {ETOT=atof(buf); break;}
		if(buf[0]=='E' && buf[1]=='T' && buf[2]=='O' && buf[3]=='T') {
			find = true;
		}

	}while(!fin.eof());
	fin.clear();
	fin.close();

	if(ETOT == 1E20) cout << "error!" << endl;

	// clear file content
	ofstream logFile(logfile, ios::out);
	logFile.close();

	return ETOT;
}


double getRMSD(double *reference, double *mobile, int length) { // Root mean square deviation

	double rms = 0.0;
	double mes = 0.0;
	for(int i=0; i<length; i++) {
		mes += (fabs(reference[i])-fabs(mobile[i]))*(fabs(reference[i])-fabs(mobile[i]));
	}
	rms = sqrt(mes/length);
	return rms;
}

void writePDB(double* geomData, char* fileName_PDB) {
	const char* inpFile = "writePDB.inp";
	const char* mainoutFile = "main_out.MetPDB";
	const char* logFile  = "log";
	const char* errorFile = "err";


	// write the configuration file used for ECEPP/3 model engine
	ofstream outFile(inpFile, ios::out);
	outFile << "$CNTRL" << endl;
	outFile << "runtyp = energy" << endl;
	outFile << "PRINT_CART" << endl;
	outFile << "OUTFORMAT = PDB" << endl;
	outFile << "FILE = ";  
	outFile << fileName_PDB << endl;
	outFile << "res_code= one_letter" << endl;

	outFile << "$END" << endl << endl;

	outFile << "$SEQ" << endl;
	outFile << "H" << endl;
	outFile << "YGGFM" << endl;
	outFile << "O" << endl;
	outFile << "$END" << endl;


	outFile << "$GEOM" << endl << endl;
	for(int i=0; i<6;i++) {
		outFile	<< setw(8) // 设定field宽度
				<< setprecision(3) // 设置小数位置
				<< setiosflags(ios::showpoint) // keep trailing 0s
				<< setiosflags(ios::fixed) // 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=6; i<9;i++) {
		outFile	<< setw(8) // 设定field宽度
				<< setprecision(3) // 设置小数位置
				<< setiosflags(ios::showpoint) // keep trailing 0s
				<< setiosflags(ios::fixed) // 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=9; i<12;i++) {
		outFile	<< setw(8) // 设定field宽度
				<< setprecision(3) // 设置小数位置
				<< setiosflags(ios::showpoint) // keep trailing 0s
				<< setiosflags(ios::fixed) // 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=12; i<17;i++) {
		outFile	<< setw(8) // 设定field宽度
				<< setprecision(3) // 设置小数位置
				<< setiosflags(ios::showpoint) // keep trailing 0s
				<< setiosflags(ios::fixed) // 使用这些设置
				<< geomData[i];
	}
	outFile << endl;

	for(i=17; i<24;i++) {
		outFile	<< setw(8) // 设定field宽度
				<< setprecision(3) // 设置小数位置
				<< setiosflags(ios::showpoint) // keep trailing 0s
				<< setiosflags(ios::fixed) // 使用这些设置
				<< geomData[i];
	}
	outFile << endl;
	outFile << endl;
	outFile << "$END" << endl;

	outFile.flush();
	outFile.close();
	system("i.bat ENERGY writePDB MetPDB x x 1>> log 2> err");

	ofstream outFile1(inpFile, ios::out);
	outFile1.flush();
	outFile1.close();

	ofstream outFile2(mainoutFile, ios::out);
	outFile2.flush();
	outFile2.close();
	
	ofstream outFile3(logFile, ios::out);
	outFile3.flush();
	outFile3.close();

	ofstream outFile4(errorFile, ios::out);
	outFile4.flush();
	outFile4.close();
}
double computeRMSD(char* ref) {
	
	char str[500] = {0};
	char numstr[100];
	strcpy(str, "python c:/python24/match.py -n pdb/ref.pdb pdb/minimize/met");
	double rmsd = 0;

	
	for(int i=1; i<3; i++) {
		itoa(i, numstr, 10);
		strcat(str, numstr);
		strcat(str, ".pdb s");

		system(str);
		cout << ".";
		/*
		ifstream fin("rmsd.txt");
		char wordBuf[200]={0};
		fin >> wordBuf;
		rmsd = atof(wordBuf);
		cout << "rmsd=" << rmsd << endl;
		fin.close();
		*/
	}

	return rmsd;
}


void readGeomFromText(char *fileName)  {
	ifstream fin(fileName);
	double* mobile = new double[24];
	char wordBuf[200]={0};

	while(!fin.eof()) {
		fin >> wordBuf;
		if(!fin.eof()) {
			char newFileName[100] = {0};
			strcpy(newFileName, "Met");
			strcat(newFileName, wordBuf);

			double ordinate = atof(wordBuf);
			fin >> wordBuf;
			double object   = atof(wordBuf);
			
			for(int i=0; i<24; i++) {
				fin>>wordBuf;
				mobile[i] =atof(wordBuf); 
			}
			
			writePDB(mobile, newFileName);
			cout << ordinate << "\t" << object << endl;
		}
	}
	delete [] mobile;
	fin.close();
}






void readMobile(char* fileName, int variableLength) {

	double *mobile = new double[variableLength];
	double *reference = new double[variableLength];

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
	reference[0] = Tyr_phi;
	reference[1] = Tyr_psi;
	reference[2] = Tyr_omega;
	reference[3] = Tyr_chi1;
	reference[4] = Tyr_chi2;
	reference[5] = Tyr_chi3;

	reference[6] = Gly1_phi;
	reference[7] = Gly1_psi;
	reference[8] = Gly1_omega;

	reference[9] = Gly2_phi;
	reference[10]= Gly2_psi;
	reference[11]= Gly2_omega;

	reference[12]= Phe_phi;
	reference[13]= Phe_psi;
	reference[14]= Phe_omega;
	reference[15]= Phe_chi1;
	reference[16]= Phe_chi2;

	reference[17]= Met_phi;
	reference[18]= Met_psi;
	reference[19]= Met_omega;
	reference[20]= Met_chi1;
	reference[21]= Met_chi2;
	reference[22]= Met_chi3;
	reference[23]= Met_chi4;

	ifstream fin(fileName);
	ofstream outfile("process.txt",   ios::out);


	char wordBuf[200]={0};
	while(!fin.eof()) {
		fin >> wordBuf;
		if(!fin.eof()) {
			double ordinate = atof(wordBuf);
			fin >> wordBuf;
			double object   = atof(wordBuf);
			
			fin >> wordBuf;
			fin >> wordBuf;
			fin >> wordBuf;

			for(int i=0; i<24; i++) {
				fin>>wordBuf;
				mobile[i] =atof(wordBuf); 
			}
			outfile << ordinate << "\t" << object << "\t\t" << getRMSD(reference, mobile, variableLength) << "\t";
			for(i=0; i<24; i++) {
				outfile << mobile[i] << "\t";
			}
			outfile << endl;
		}
	}
	cout.flush();

	delete [] reference;
	delete [] mobile;

	outfile.close();
	fin.close();
}




void generatePDBFromFiles(const char* fileName) {
	char newFileName[50] = {0};
	strcpy(newFileName, fileName);
	strcat(newFileName, "_");

	ifstream inFile  (fileName);					
	ofstream outFile (newFileName, ios::out);
	long outFileIndex = 0;

	const int LINE_LENGTH = 1200; 
	char   wordBuf[100]  = {0};
	char   lineBuf[LINE_LENGTH] = {0};
	char   substr[9]={0};
	double geomData[24] = {0};
	double ETOT = 0;
	
	int readTime = 201;
	int readSize = 501;

	for(int m = 1; m<readTime; m++) {
		/***********************************************************
		// write inpFile
		************************************************************/
		char midName[10] = {0};
		char pdbFileName[30] = {0};
		itoa(m, midName, 10);
		strcpy(pdbFileName, "Met_");
		strcat(pdbFileName, midName);
		strcat(pdbFileName, "_");

		ofstream inpFile(MULTIPLE_ENERGY_EVALUATION_SOURCE_FILE, ios::out);
		inpFile << "$CNTRL" << endl;
		inpFile << "runtyp = energy" << endl;
		inpFile << "PRINT_CART" << endl;
		inpFile << "OUTFORMAT = PDB" << endl;
		inpFile << "FILE = ";  
		inpFile << pdbFileName << endl;
		inpFile << "res_code= one_letter" << endl;
		inpFile << "EMINIMA = 0.100E+35" << endl;
		inpFile << "$END" << endl << endl;

		inpFile << "$SEQ" << endl;
		inpFile << "H" << endl;
		inpFile << "YGGFM" << endl;
		inpFile << "O" << endl;
		inpFile << "$END" << endl << endl;

		inpFile << "$ENERCALC" << endl;
		inpFile << "READ_CONF" << endl;
		inpFile << "$END" << endl << endl;

		inpFile << "$GEOM" << endl;
		inpFile << endl << endl << endl << endl <<endl <<endl << endl;
		inpFile << "$END" << endl;

		inpFile.flush( );
		inpFile.close( );

		/***********************************************************
		// write input outo file
		************************************************************/
		ofstream inputFile(MULTIPLE_ENERGY_EVALUATION_INPUT_FILE, ios::out);
		for(int i=1; i<readSize; i++) {
			inFile >> wordBuf;								// 序号
			inFile >> wordBuf;								// 评价值
			ETOT = atof(wordBuf);
			for(int k=0; k<24; k++) {
				inFile >> wordBuf;
				geomData[k] = atof(wordBuf);
			}

			inputFile << i << endl;
			inputFile << "   1  20   6   6   5  11  11" << endl;
		
			inputFile << endl;

			for(k=0; k<6; k++) {
				inputFile << setw(8)						// 设定field宽度
						  << setprecision(3)				// 设置小数位置
						  << setiosflags(ios::showpoint)	// keep trailing 0s
						  << setiosflags(ios::fixed)		// 使用这些设置
						  << geomData[k];
			}
			inputFile << endl;

			for(k=6; k<9; k++) {
				inputFile << setw(8)						// 设定field宽度
						  << setprecision(3)				// 设置小数位置
						  << setiosflags(ios::showpoint)	// keep trailing 0s
						  << setiosflags(ios::fixed)		// 使用这些设置
						  << geomData[k];
			}
			inputFile << endl;

			for(k=9; k<12; k++) {
				inputFile << setw(8)						// 设定field宽度
						  << setprecision(3)				// 设置小数位置
						  << setiosflags(ios::showpoint)	// keep trailing 0s
						  << setiosflags(ios::fixed)		// 使用这些设置
						  << geomData[k];
			}
			inputFile << endl;

			for(k=12; k<17; k++) {
				inputFile << setw(8)						// 设定field宽度
						  << setprecision(3)				// 设置小数位置
						  << setiosflags(ios::showpoint)	// keep trailing 0s
						  << setiosflags(ios::fixed)		// 使用这些设置
						  << geomData[k];
			}
			inputFile << endl;

			for(k=17; k<24; k++) {
				inputFile << setw(8)						// 设定field宽度
						  << setprecision(3)				// 设置小数位置
						  << setiosflags(ios::showpoint)	// keep trailing 0s
						  << setiosflags(ios::fixed)		// 使用这些设置
						  << geomData[k];
			}
			inputFile << endl;
			inputFile << endl;
			inputFile.flush();
		}
		inputFile.close();	
		
		/***********************************************************
		// call ecepp model script
		************************************************************/
		system("i.bat ENERGY metEnkephalinEnergyMultipleEvaluationSourceFile MET_MULTI_EVAL MultipleEnergyEvaluateInput rr  1>> log 2> err");
		
		/***********************************************************
		// read result and save as newFileName
		************************************************************/
		ifstream resultFile("outo.MET_MULTI_EVAL");

		for(i=1; i<readSize; i++) {
			double saveGeom[24] = {0};
			int index =0;

			resultFile >> wordBuf;
			resultFile >> wordBuf;
			ETOT = (float)atof(wordBuf);
			resultFile.getline(lineBuf, LINE_LENGTH);	// ETOT line
			resultFile.getline(lineBuf, LINE_LENGTH);	// sequence line 
			resultFile.getline(lineBuf, LINE_LENGTH);	// end group-H line

			resultFile.getline(lineBuf, LINE_LENGTH);	// Tyr
			for(int m1=0; m1<6;m1++) {
				for(int n1=0; n1<8;n1++) {
					substr[n1] = lineBuf[m1*8+n1];
				}
				saveGeom[index++] = atof(substr);
			}

			resultFile.getline(lineBuf, LINE_LENGTH);	// Gly1
			for(m1=0; m1<3;m1++) {
				for(int n1=0; n1<8;n1++) {
					substr[n1] = lineBuf[m1*8+n1];
				}
				saveGeom[index++] = atof(substr);
			}

			resultFile.getline(lineBuf, LINE_LENGTH);	// Gly2
			for(m1=0; m1<3;m1++) {
				for(int n1=0; n1<8;n1++) {
					substr[n1] = lineBuf[m1*8+n1];
				}
				saveGeom[index++] = atof(substr);
			}

			resultFile.getline(lineBuf, LINE_LENGTH);	// Phe
			for(m1=0; m1<5;m1++) {
				for(int n1=0; n1<8;n1++) {
					substr[n1] = lineBuf[m1*8+n1];
				}
				saveGeom[index++] = atof(substr);
			}

			resultFile.getline(lineBuf, LINE_LENGTH);	// Met
			for(m1=0; m1<7;m1++) {
				for(int n1=0; n1<8;n1++) {
					substr[n1] = lineBuf[m1*8+n1];
				}
				saveGeom[index++] = atof(substr);
			}
			
			resultFile.getline(lineBuf, LINE_LENGTH);	// eng group -COOH

			outFile << ++outFileIndex << "\t" << ETOT << "\t";
			for(int oi=0; oi<24; oi++) {
				outFile << saveGeom[oi] << "\t";
			}
			outFile<<endl;
			outFile.flush();
		}
		resultFile.close( );

		/***********************************************************
		// clear trace
		************************************************************/
		ofstream outFile1("outo.MET_MULTI_EVAL", ios::out);
		outFile1.flush();
		outFile1.close();

		ofstream outFile2("main_out.MET_MULTI_EVAL", ios::out);
		outFile2.flush();
		outFile2.close();

		ofstream outFile3(MULTIPLE_ENERGY_EVALUATION_SOURCE_FILE, ios::out);
		outFile3.flush();
		outFile3.close();

		ofstream outFile4(MULTIPLE_ENERGY_EVALUATION_INPUT_FILE, ios::out);
		outFile4.flush();
		outFile4.close();

		ofstream outFile5("log", ios::out);
		outFile5.flush();
		outFile5.close();

		/***********************************************************
		// moving .pdb to file fold
		************************************************************/
		system("md Met");
		system("move Met_*.pdb Met");
	}

	inFile.close();
	outFile.close();
}




int	main(int argc, char** argv) {
	cout << "MetEnkephalin Optimization Beging...\n"; cout.flush();
		
	GAGenome *_theGenome = 0;
	GAGeneticAlgorithm *_theGA = 0;
	
	unsigned int seed = 0;
	for(int ii=1; ii<argc; ii++) {
		if(strcmp(argv[ii++],"seed") == 0) {
		seed = atoi(argv[ii]);
		}
	}
	
	// navtive conformation
	double Tyr_phi		=  -83.466;
	double Tyr_psi		=	155.791;
	double Tyr_omega	=  -177.129;
	double Tyr_chi1		=  -173.181;
	double Tyr_chi2		=	 79.338;
	double Tyr_chi3		=  -166.326;
	
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
	
	GARealAlleleSetArray array;
	array.add(-180,	180);	// 酪氨酸:   Tyr_phi
	array.add(-180,	180);	// 酪氨酸:   Tyr_psi
	array.add(-180,	180);	// 酪氨酸:   Tyr_omega
	array.add(-180,	180);	// 酪氨酸:   Tyr_chi1
	array.add(-180,	180);	// 酪氨酸:   Tyr_chi2
	array.add(-180,	180);	// 酪氨酸:   Tyr_chi3
	
	array.add(-180,	180);	// 甘氨酸:   Gly1_phi
	array.add(-180,	180);	// 甘氨酸:   Gly1_psi
	array.add(-180,	180);	// 甘氨酸:   Gly1_omega
	
	array.add(-180,	180);	// 甘氨酸:   Gly2_phi
	array.add(-180,	180);	// 甘氨酸:   Gly2_psi
	array.add(-180,	180);	// 甘氨酸:   Gly2_omega

	array.add(-180,	180);	// 苯丙氨酸: Phe_phi
	array.add(-180,	180);	// 苯丙氨酸: Phe_psi
	array.add(-180,	180);	// 苯丙氨酸: Phe_omega
	array.add(-180,	180);	// 苯丙氨酸: Phe_chi1
	array.add(-180,	180);	// 苯丙氨酸: Phe_chi2

	array.add(-180,	180);	// 甲硫氨酸: Met_phi
	array.add(-180,	180);	// 甲硫氨酸: Met_psi
	array.add(-180,	180);	// 甲硫氨酸: Met_omega
	array.add(-180,	180);	// 甲硫氨酸: Met_chi1
	array.add(-180,	180);	// 甲硫氨酸: Met_chi2
	array.add(-180,	180);	// 甲硫氨酸: Met_chi3
	array.add(-180,	180);	// 甲硫氨酸: Met_chi4

	// each genome contains a generic pointer to user-specifiable data.
	// all genomes in a population refer to the same user data.
	FunctionStruct FS;
	FS.objectiveFunction = metEnkePhalinSingleEnergyEvaluation;
	FS.optimum = FS.lowerBound = -100;
	// used to store genome-specific evaluation data. 
	MyEvalData med;

	_theGenome = new GARealGenome(array, genomeEvaluate, &FS);
	_theGenome->initializer(GARealUniformInitializer);
	
	_theGA = new CDifferentialEvolution(*_theGenome);
	_theGA->minimize();
	_theGA->pMutation(0.2);
	_theGA->pCrossover(0.8);
	_theGA->nGenerations(400);
	_theGA->populationSize(100);
	CDifferentialEvolution& de = dynamic_cast<CDifferentialEvolution&>(*_theGA);
	de.objectiveData(med);



/*
	CSAAlgorithm *csaAlgorithm = new CSAAlgorithm(*_theGenome);
	csaAlgorithm->minimize();
	csaAlgorithm->nGenerations(10);
	csaAlgorithm->populationSize(5);
	csaAlgorithm->objectiveData(med);
	csaAlgorithm->initialize(seed); 
	csaAlgorithm->step();

	cout<<"abnormal exited!"<<endl; return 0; 
*/

	//de.randomMinimizeRsmdAnalysis();
	//return 0;

	//de.dihedralAngleTest( );
	//de.randomMinimize( );

	//double x = de.computeAnnealingX(1677, 5000, 90);
	//de.distanceCuttingCurve("testData.txt", 1677, x);
	//return 0;

	//de.conformationSpaceAnnealingOptimize();
	//return 0;

	//de.randomMinimizeRsmdAnalysis( );
	//return 0;

	//de.variableMinimizeSmoothAnalysis1( );
	//return 0;

	//de.randomMinimizeDihedralAnalysis( );
	//return 0;

	int runNumbers = 1;
	double* averageScore					= new double[_theGA->nGenerations()];
	double* averageEnergyEvals				= new double[_theGA->nGenerations()];
	double* averageEnergyEvalsCPUTime		= new double[_theGA->nGenerations()];
	double* averageNumReps					= new double[_theGA->nGenerations()];
	for(int asIndex = 0; asIndex<_theGA->nGenerations(); asIndex++) {
		averageScore[asIndex]				= 0.0;
		averageEnergyEvals[asIndex]			= 0.0;
		averageEnergyEvalsCPUTime[asIndex]	= 0.0;
		averageNumReps[asIndex]				= 0.0;
	}
	
	for(int i=0; i<runNumbers; i++) {
		cout << "runNumber: " << i+1 <<  endl;
		_theGA->initialize(seed); 
		while(!_theGA->done()){
			
			CDifferentialEvolution& de = dynamic_cast<CDifferentialEvolution&>(*_theGA);
            de.step5();
			cout<<_theGA->generation()<<"\t";			//输出“迭代次数”信息；
			cout<< setw(8)								//输出“最好分值”信息；
				<< setprecision(4) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< _theGA->population().best().score();
			cout<< setw(12)								//输出“评价次数”信息；
				<< setprecision(0) 
				<< setiosflags(ios::fixed)
				<< de.energyEvals;
			cout<< setw(12)								//输出“计算时间”信息； 
				<< setprecision(3) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< de.energyEvalsCPUTime;
			cout<< setw(12)								//输出“替换次数”信息
				<< setprecision(0) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< _theGA->statistics().numrep;
			cout<< setw(12)								//输出“平均距离”信息 
				<< setprecision(0) 
				<< setiosflags(ios::showpoint) 
				<< setiosflags(ios::fixed)
				<< de.populationAverageDistance(de.population());
			cout<< endl;
			cout.flush();
			
			averageScore[_theGA->generation()-1] += _theGA->population().best().score();
			averageEnergyEvals[_theGA->generation()-1] += de.energyEvals;
			averageEnergyEvalsCPUTime[_theGA->generation()-1] += de.energyEvalsCPUTime;
			averageNumReps[_theGA->generation()-1] += _theGA->statistics().numrep;
		}
	}

	
	ofstream averageFile("D:\\MetEnkephalinOptimization\\data\\averageIndicators.txt", ios::out);
	averageFile << "gen" << "\t" << "score" << "\t" << "evals" << "\t" << "times" << "\t" << "reps" << "\n";
	for(int ai=0; ai<_theGA->nGenerations(); ai++) {
		averageScore[ai] = averageScore[ai]/runNumbers;
		averageEnergyEvals[ai] = averageEnergyEvals[ai]/runNumbers;
		averageEnergyEvalsCPUTime[ai] = averageEnergyEvalsCPUTime[ai]/runNumbers;
		averageNumReps[ai] = averageNumReps[ai]/runNumbers;

		averageFile << ai+1								<< "\t" 
					<< averageScore[ai]					<< "\t"
					<< averageEnergyEvals[ai]			<< "\t"
					<< averageEnergyEvalsCPUTime[ai]	<< "\t"
					<< averageNumReps[ai]				<< "\n";
	}
	averageFile.close();


	delete [] averageScore;
	delete [] averageEnergyEvals;
	delete [] averageEnergyEvalsCPUTime;
	delete [] averageNumReps;

	return 0;
}
 




