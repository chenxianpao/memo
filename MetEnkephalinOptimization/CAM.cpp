// CAM.cpp: implementation of the CCAM class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "CAM.h"
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include <fstream.h>







//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


CCAM::CCAM()
{

}

CCAM::~CCAM()
{

}

double CCAM::function(double *x, int dim) {
	return 1;
}


void CCAM::updateLeafArray(CLLeafArray& leafArray, Cl& l) {
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
			CL L(cl, tempCL->m_dim);
			if (L.getCondition1()) {
				if(!illposed(L)) {
					leafArray.add(L); 
				}
				else {
					int ii  = 0;
				}
			}
		}
		delete [] cl;
	}
}

bool CCAM::illposed(CL& arg) {
	bool ret = false;
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

	for(int j=0; j<dim; j++) {		// 避免陷入松弛问题局部最优解
		if (arg.m_l[j]->getElementData(j) < 0.0000000001) {			//0.0000000001
			small = true;
			break;
		}
		if(small) break;
	} 

	//if( small && count==2) {
	//	ret = true;
	//}
	if (small) ret = true;
	return ret;
}

void CCAM::optimize( ){
	Evaluator FuncEvaluator	= function; 
	
	double	CONSTANT =20;			// 50000	
	int		DIM = 2;				// 问题维数
	
	int		dim = DIM;
	int		k   = DIM;
	int		K   = 10000;
	double	f   = 0;				// 目标函数		

	CL		*TargetL	= 0;
	double	TargetLMin	= 0;
	double	*tempx		= new double [dim];
	double	*templ		= new double [dim];
	int i, j;

	ClArray		m_ClArray(dim);			// 支撑向量集合
	CLLeafArray	m_ClLeafArray(dim);		// 树叶集合


	// 初始化I:生成初始支撑向量l1,l2,...,ln
	Cl *m_pCl = 0;
	for (i=0; i<dim; i++) {
		for (int j=0; j<dim; j++) {tempx[j] = 0;}
		tempx[i] = 1;
		f = energyFunction(tempx, dim)+CONSTANT;

		//f = FuncEvaluator(tempx, dim)+CONSTANT;
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
		// 计算获取最小树叶
		TargetL	= m_ClLeafArray.findMinimumLeaf();
		// 计算最小树叶对应最低值
		TargetLMin = TargetL->getMinimum( );
		// 计算最小树叶对应变量值
		for(j=0; j<dim; j++) {tempx[j] = TargetL->getDialogElement(j) * TargetLMin;}
		// 计算最小树叶形成支撑向量l^k+1
//		f = FuncEvaluator(tempx, dim) + CONSTANT;
		f = energyFunction(tempx, dim) + CONSTANT;
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
		updateLeafArray(m_ClLeafArray, *m_ClArray.m_l[k]);
	}
	delete [] tempx;
	delete [] templ;
}

double CCAM::energyFunction(double *x, int dim) {
	double ret = 0;

	double Tyr_phi = 360*x[0]-180;
	//double Tyr_phi	=	-73.453;
/*	
	double Tyr_psi	=	55.802;
	double Tyr_omega=	-177.155;
	double Tyr_chi1	=	-73.178;
	double Tyr_chi2	=	79.396;
	double Tyr_chi3	=	-166.361;

	double Gly1_phi	=	-54.320;
	double Gly1_psi	=	80.925;
	double Gly1_omega=	68.504;

	double Gly2_phi	=	82.967;
	double Gly2_psi	=	-5.095;
	double Gly2_omega=	-19.983;

	double Phe_phi	=	-136.8;
	double Phe_psi	=	19.100;
	double Phe_omega=	-14.090;
	double Phe_chi1	=	58.864;
	double Phe_chi2	=	-85.474;

	double Met_phi	=	-13.448;
	double Met_psi	=	11.000;
	double Met_omega=	-19.791;
	double Met_chi1	=	52.873;
	double Met_chi2	=	125.297;
	double Met_chi3	=	-149.869;
	double Met_chi4	=	-56.583;
*/
	double Tyr_psi	=	155.802;
	double Tyr_omega=	-177.155;
	double Tyr_chi1	=	-173.178;
	double Tyr_chi2	=	79.396;
	double Tyr_chi3	=	-166.361;

	double Gly1_phi	=	-154.320;
	double Gly1_psi	=	85.925;
	double Gly1_omega=	168.504;

	double Gly2_phi	=	82.967;
	double Gly2_psi	=	-75.095;
	double Gly2_omega=	-169.983;

	double Phe_phi	=	-136.839;
	double Phe_psi	=	19.100;
	double Phe_omega=	-174.090;
	double Phe_chi1	=	58.864;
	double Phe_chi2	=	-85.474;

	double Met_phi	=	-163.448;
	double Met_psi	=	161.000;
	double Met_omega=	-179.791;
	double Met_chi1	=	52.873;
	double Met_chi2	=	175.297;
	double Met_chi3	=	-179.869;
	double Met_chi4	=	-58.583;

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
	
	ret = singleEnergyMinimize(pGeom, 24);

	delete [] pGeom;
	return ret;

}


double CCAM::singleEnergyMinimize(double* pGEOM, int length) {
	double ETOT = 1E25;
	writeSingleEnergyMinimizeSourceFile(SINGLE_ENERGY_MINIMIZE_SOURCE_FILE, pGEOM, length);
	system("i.bat MINIMIZE singleEnergyMinimizeSourceFile MET_SINGLE_MIN x x  1>> log 2> err");
	ETOT = readSingleEnergyMinimizeResultETOT(SINGLE_ENERGY_MINIMIZE_OUTO_FILE, pGEOM, length);
	clearSingleEnergyMinimizeResultFiles();
	return ETOT;
}

void CCAM::writeSingleEnergyMinimizeSourceFile(const char* sourceFile, double* geomData, int length){

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

double CCAM::readSingleEnergyMinimizeResultETOT(const char* resultFile, double* pGEOM, int length) {
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


void CCAM::clearSingleEnergyMinimizeResultFiles( ) {
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
