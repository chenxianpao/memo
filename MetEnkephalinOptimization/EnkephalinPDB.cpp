#include "stdafx.h"
#include "EnkephalinPDB.h"


EnkephalinPDB::EnkephalinPDB() {
	atomNumberOfTotalPDB = 0;
	energy = 0;
	pAtomEntry=new CAtomEntry* [75];					// at most have 75 atoms, for simple just set 75;
	for(int i=0; i<24; i++) dihedralAngle[i] = 0.0;
}

EnkephalinPDB::~EnkephalinPDB() {
	if(pAtomEntry) {
		for(int i=0; i<atomNumberOfTotalPDB;i++) {
			delete pAtomEntry[i];
		}
	}
	delete [] pAtomEntry;
	pAtomEntry = 0;
}

int EnkephalinPDB::readPDB(char* PDBFileName) {
	const int LINE_LENGTH = 80;
	const int CHAR_LENGTH = 10;
	char lineBuf[LINE_LENGTH] ={0};		
	char charBuf[CHAR_LENGTH] ={0};

	ifstream fin(PDBFileName);
	setPDBFileName(PDBFileName);

	atomNumberOfTotalPDB = 0;
	while( fin.getline(lineBuf, LINE_LENGTH) ) {
		strncpy(charBuf, lineBuf, 6);
		if(!strncmp(charBuf, "REMARK", 6)) {			// strict format, can be improved further...
			char   strvalue[13] = {0};
			double energy = 0;
			for(int m=0; m<13; m++) {
				strvalue[m] = lineBuf[m+31];
			}
			setEnergy(atof(strvalue));
		} 

		if(!strncmp(charBuf, "ATOM  " , 6)) {
			CAtomEntry* pAtomEntry= new CAtomEntry( );	// managed by container; do delete it!!!
			
			pAtomEntry->setRecName(charBuf);
			
			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+6, 5);
			pAtomEntry->setSerial(charBuf);

			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+12, 4);
			pAtomEntry->setAtomName(charBuf);

			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+16,1);
			pAtomEntry->setAltLoc(charBuf);

			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+17,3);
			pAtomEntry->setResName(charBuf);

			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+21,1);
			pAtomEntry->setChainID(charBuf);

			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+22,4);
			pAtomEntry->setSeqNo(charBuf);

			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+26,1);
			pAtomEntry->setInsertCode(charBuf);

			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+30,8);
			pAtomEntry->setX(charBuf);
			
			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+38,8);
			pAtomEntry->setY(charBuf);
	
			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+46,8);
			pAtomEntry->setZ(charBuf);
			
			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+54,6);
			pAtomEntry->setOccupancy(charBuf);

			strncpy(charBuf, "\0", CHAR_LENGTH);
			strncpy(charBuf, lineBuf+60,6);
			pAtomEntry->setTempFactor(charBuf);
			add(pAtomEntry);
		}
	}

	fin.close( );
	return atomNumberOfTotalPDB;
}

int EnkephalinPDB::writePDB( ) {
	if(atomNumberOfTotalPDB < 1) return 0;

	cout << "OUTPUT FORMAT: PDB MIRROR" << endl;
	cout << "Energy: "  << getEnergy( ) << endl;

	for(int i=0; i<atomNumberOfTotalPDB; i++) {
		cout << getEntryByIndex(i)->getSerial() << "\t"
			 << getEntryByIndex(i)->getResName()<< "-"
			 << getEntryByIndex(i)->getAtomName()<<"\t"
			 << setw(8)<< setprecision(3)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)<< getEntryByIndex(i)->getX()
			 << setw(8)<< setprecision(3)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)<< getEntryByIndex(i)->getY()
			 << setw(8)<< setprecision(3)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)<< getEntryByIndex(i)->getZ()
			 << setw(6)<< setprecision(2)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)<< getEntryByIndex(i)->getOccupancy() 
			 << setw(6)<< setprecision(2)<< setiosflags(ios::showpoint)<< setiosflags(ios::fixed)<< getEntryByIndex(i)->getTempFactor()
			 << endl;
	}
	cout.flush();
	return 1;
}


double ** EnkephalinPDB::extractRMSDTypeArray(RMSDTYPE rmsdType, int& numberOfTypeAtom) {
	double** array = 0;
	int i, j, k;

	if (atomNumberOfTotalPDB<1) return 0;

	if(rmsdType == AlphaC) {
		numberOfTypeAtom = 5;
		array = new double*[numberOfTypeAtom];
		for(i=0; i<numberOfTypeAtom; i++) {
			array[i] = new double[3];
			for(j=0; j<3;j++) array[i][j] = 0.0;
		}
		k = 0;
		array[k][0]		=(getEntryByType(TYR_CA))->getX();
		array[k][1]		=(getEntryByType(TYR_CA))->getY();
		array[k++][2]	=(getEntryByType(TYR_CA))->getZ();

		array[k][0]		=(getEntryByType(GLY1_CA))->getX();
		array[k][1]		=(getEntryByType(GLY1_CA))->getY();
		array[k++][2]	=(getEntryByType(GLY1_CA))->getZ();

		array[k][0]		=(getEntryByType(GLY2_CA))->getX();
		array[k][1]		=(getEntryByType(GLY2_CA))->getY();
		array[k++][2]	=(getEntryByType(GLY2_CA))->getZ();

		array[k][0]		=(getEntryByType(PHE_CA))->getX();
		array[k][1]		=(getEntryByType(PHE_CA))->getY();
		array[k++][2]	=(getEntryByType(PHE_CA))->getZ();

		array[k][0]		=(getEntryByType(MET_CA))->getX();
		array[k][1]		=(getEntryByType(MET_CA))->getY();
		array[k++][2]	=(getEntryByType(MET_CA))->getZ();
	}else if(rmsdType == Backbone) {
		numberOfTypeAtom = 15;
		array = new double*[numberOfTypeAtom];
		for(i=0; i<numberOfTypeAtom; i++) {
			array[i] = new double[3];
			for(j=0; j<3;j++) array[i][j] = 0.0;
		}
		k = 0;
		array[k][0]		=(getEntryByType(TYR_N))->getX();
		array[k][1]		=(getEntryByType(TYR_N))->getY();
		array[k++][2]	=(getEntryByType(TYR_N))->getZ();
		array[k][0]		=(getEntryByType(TYR_CA))->getX();
		array[k][1]		=(getEntryByType(TYR_CA))->getY();
		array[k++][2]	=(getEntryByType(TYR_CA))->getZ();
		array[k][0]		=(getEntryByType(TYR_C))->getX();
		array[k][1]		=(getEntryByType(TYR_C))->getY();
		array[k++][2]	=(getEntryByType(TYR_C))->getZ();

		array[k][0]		=(getEntryByType(GLY1_N))->getX();
		array[k][1]		=(getEntryByType(GLY1_N))->getY();
		array[k++][2]	=(getEntryByType(GLY1_N))->getZ();
		array[k][0]		=(getEntryByType(GLY1_CA))->getX();
		array[k][1]		=(getEntryByType(GLY1_CA))->getY();
		array[k++][2]	=(getEntryByType(GLY1_CA))->getZ();
		array[k][0]		=(getEntryByType(GLY1_C))->getX();
		array[k][1]		=(getEntryByType(GLY1_C))->getY();
		array[k++][2]	=(getEntryByType(GLY1_C))->getZ();

		array[k][0]		=(getEntryByType(GLY2_N))->getX();
		array[k][1]		=(getEntryByType(GLY2_N))->getY();
		array[k++][2]	=(getEntryByType(GLY2_N))->getZ();
		array[k][0]		=(getEntryByType(GLY2_CA))->getX();
		array[k][1]		=(getEntryByType(GLY2_CA))->getY();
		array[k++][2]	=(getEntryByType(GLY2_CA))->getZ();
		array[k][0]		=(getEntryByType(GLY2_C))->getX();
		array[k][1]		=(getEntryByType(GLY2_C))->getY();
		array[k++][2]	=(getEntryByType(GLY2_C))->getZ();

		array[k][0]		=(getEntryByType(PHE_N))->getX();
		array[k][1]		=(getEntryByType(PHE_N))->getY();
		array[k++][2]	=(getEntryByType(PHE_N))->getZ();
		array[k][0]		=(getEntryByType(PHE_CA))->getX();
		array[k][1]		=(getEntryByType(PHE_CA))->getY();
		array[k++][2]	=(getEntryByType(PHE_CA))->getZ();
		array[k][0]		=(getEntryByType(PHE_C))->getX();
		array[k][1]		=(getEntryByType(PHE_C))->getY();
		array[k++][2]	=(getEntryByType(PHE_C))->getZ();

		array[k][0]		=(getEntryByType(MET_N))->getX();
		array[k][1]		=(getEntryByType(MET_N))->getY();
		array[k++][2]	=(getEntryByType(MET_N))->getZ();
		array[k][0]		=(getEntryByType(MET_CA))->getX();
		array[k][1]		=(getEntryByType(MET_CA))->getY();
		array[k++][2]	=(getEntryByType(MET_CA))->getZ();
		array[k][0]		=(getEntryByType(MET_C))->getX();
		array[k][1]		=(getEntryByType(MET_C))->getY();
		array[k++][2]	=(getEntryByType(MET_C))->getZ();
	}else if(rmsdType == HeavyAtom) {
			// do support by now!

	}else {
		numberOfTypeAtom = 75;
		array = new double*[numberOfTypeAtom];
		for(i=0; i<numberOfTypeAtom; i++) {
			array[i] = new double[3];
			for(j=0; j<3;j++) array[i][j] = 0.0;
		}
		for(i=0; i<numberOfTypeAtom; i++) {
			array[i][0] = (getEntryByIndex(i))->getX( );
			array[i][1] = (getEntryByIndex(i))->getY( );
			array[i][2] = (getEntryByIndex(i))->getZ( );
		}
	}

	return array;
}

void EnkephalinPDB::freeArray(double ** xlist, int number) {
	for(int i=0; i<number; i++) {
		double *a = xlist[i];
		delete [] a;
	}
	delete [] xlist;
}

void EnkephalinPDB::combineFileName(char* dictory, 
									char* firstName, 
									char*suffixName, 
									int middleIndex, 
									char* fullFileName) 
{
	char midName[5]	= {0};
	char mName[5]   = {0};
	itoa(middleIndex, mName, 10); 
	if(middleIndex<10) {
		strcpy(midName, "00");
		strcat(midName, mName);
	}
	else if(middleIndex<100) {
		strcpy(midName, "0");
		strcat(midName, mName);
	}else if(middleIndex<1000) {
		strcat(midName, mName);
	}

	strcpy(fullFileName, dictory);
	strcat(fullFileName, firstName);
	strcat(fullFileName, midName);
	strcat(fullFileName, suffixName);
}


int EnkephalinPDB::computeDihedral( ) {
	if (atomNumberOfTotalPDB<1) return 0;

	CAtomEntry *pEntry = 0;
	double a1[3] = {0};
	double a2[3] = {0};
	double a3[3] = {0};
	double a4[3] = {0};

	/*
	! TYROSINE
	20   6
    0 C         1 N         1 CA        1 C          1   phi
    1 N         1 CA        1 C         2 N          2   psi
    1 CA        1 C         2 N         2 CA         3   omega
    1 N         1 CA        1 CB        1 CG         4   chi1
    1 CA        1 CB        1 CG        1 CD1        5   chi2
    1 CE1       1 CZ        1 OH        1 HH         6   chi3
	*/

	// TYR
	initArray(a1, *getEntryByType(TYR_NH2));
	initArray(a2, *getEntryByType(TYR_N));
	initArray(a3, *getEntryByType(TYR_CA));
	initArray(a4, *getEntryByType(TYR_C));
	dihedralAngle[0] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(TYR_N));
	initArray(a2, *getEntryByType(TYR_CA));
	initArray(a3, *getEntryByType(TYR_C));
	initArray(a4, *getEntryByType(GLY1_N));
	dihedralAngle[1] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(TYR_CA));
	initArray(a2, *getEntryByType(TYR_C));
	initArray(a3, *getEntryByType(GLY1_N));
	initArray(a4, *getEntryByType(GLY1_CA));
	dihedralAngle[2] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(TYR_N));
	initArray(a2, *getEntryByType(TYR_CA));
	initArray(a3, *getEntryByType(TYR_CB));
	initArray(a4, *getEntryByType(TYR_CG));
	dihedralAngle[3] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(TYR_CA));
	initArray(a2, *getEntryByType(TYR_CB));
	initArray(a3, *getEntryByType(TYR_CG));
	initArray(a4, *getEntryByType(TYR_CD1));
	dihedralAngle[4] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(TYR_CE1));
	initArray(a2, *getEntryByType(TYR_CZ));
	initArray(a3, *getEntryByType(TYR_OH));
	initArray(a4, *getEntryByType(TYR_HH));
	dihedralAngle[5] = computerDihedral(a1,a2,a3,a4);

	
	// GLY1
	initArray(a1, *getEntryByType(TYR_C));
	initArray(a2, *getEntryByType(GLY1_N));
	initArray(a3, *getEntryByType(GLY1_CA));
	initArray(a4, *getEntryByType(GLY1_C));
	dihedralAngle[6] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(GLY1_N));
	initArray(a2, *getEntryByType(GLY1_CA));
	initArray(a3, *getEntryByType(GLY1_C));
	initArray(a4, *getEntryByType(GLY2_N));
	dihedralAngle[7] = computerDihedral(a1,a2,a3,a4);
		
	initArray(a1, *getEntryByType(GLY1_CA));
	initArray(a2, *getEntryByType(GLY1_C));
	initArray(a3, *getEntryByType(GLY2_N));
	initArray(a4, *getEntryByType(GLY2_CA));
	dihedralAngle[8] = computerDihedral(a1,a2,a3,a4);

	// GLY2
	initArray(a1, *getEntryByType(GLY1_C));
	initArray(a2, *getEntryByType(GLY2_N));
	initArray(a3, *getEntryByType(GLY2_CA));
	initArray(a4, *getEntryByType(GLY2_C));
	dihedralAngle[9] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(GLY2_N));
	initArray(a2, *getEntryByType(GLY2_CA));
	initArray(a3, *getEntryByType(GLY2_C));
	initArray(a4, *getEntryByType(PHE_N));
	dihedralAngle[10] = computerDihedral(a1,a2,a3,a4);
	
	initArray(a1, *getEntryByType(GLY2_CA));
	initArray(a2, *getEntryByType(GLY2_C));
	initArray(a3, *getEntryByType(PHE_N));
	initArray(a4, *getEntryByType(PHE_CA));
	dihedralAngle[11] = computerDihedral(a1,a2,a3,a4);
	
	// PHE
	initArray(a1, *getEntryByType(GLY2_C));
	initArray(a2, *getEntryByType(PHE_N));
	initArray(a3, *getEntryByType(PHE_CA));
	initArray(a4, *getEntryByType(PHE_C));
	dihedralAngle[12] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(PHE_N));
	initArray(a2, *getEntryByType(PHE_CA));
	initArray(a3, *getEntryByType(PHE_C));
	initArray(a4, *getEntryByType(MET_N));
	dihedralAngle[13] = computerDihedral(a1,a2,a3,a4);
		
	initArray(a1, *getEntryByType(PHE_CA));
	initArray(a2, *getEntryByType(PHE_C));
	initArray(a3, *getEntryByType(MET_N));
	initArray(a4, *getEntryByType(MET_CA));
	dihedralAngle[14] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(PHE_N));
	initArray(a2, *getEntryByType(PHE_CA));
	initArray(a3, *getEntryByType(PHE_CB));
	initArray(a4, *getEntryByType(PHE_CG));
	dihedralAngle[15] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(PHE_CA));
	initArray(a2, *getEntryByType(PHE_CB));
	initArray(a3, *getEntryByType(PHE_CG));
	initArray(a4, *getEntryByType(PHE_CD1));
	dihedralAngle[16] = computerDihedral(a1,a2,a3,a4);

	// MET
	initArray(a1, *getEntryByType(PHE_C));
	initArray(a2, *getEntryByType(MET_N));
	initArray(a3, *getEntryByType(MET_CA));
	initArray(a4, *getEntryByType(MET_C));
	dihedralAngle[17] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(MET_N));
	initArray(a2, *getEntryByType(MET_CA));
	initArray(a3, *getEntryByType(MET_C));
	initArray(a4, *getEntryByType(CXH_O));
	dihedralAngle[18] = computerDihedral(a1,a2,a3,a4);
	
	initArray(a1, *getEntryByType(MET_CA));
	initArray(a2, *getEntryByType(MET_C));
	initArray(a3, *getEntryByType(CXH_O));
	initArray(a4, *getEntryByType(CXH_H));
	dihedralAngle[19] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(MET_N));
	initArray(a2, *getEntryByType(MET_CA));
	initArray(a3, *getEntryByType(MET_CB));
	initArray(a4, *getEntryByType(MET_CG));
	dihedralAngle[20] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(MET_CA));
	initArray(a2, *getEntryByType(MET_CB));
	initArray(a3, *getEntryByType(MET_CG));
	initArray(a4, *getEntryByType(MET_SD));
	dihedralAngle[21] = computerDihedral(a1,a2,a3,a4);

	initArray(a1, *getEntryByType(MET_CB));
	initArray(a2, *getEntryByType(MET_CG));
	initArray(a3, *getEntryByType(MET_SD));
	initArray(a4, *getEntryByType(MET_CE));
	dihedralAngle[22] = computerDihedral(a1,a2,a3,a4);


	initArray(a1, *getEntryByType(MET_CG));
	initArray(a2, *getEntryByType(MET_SD));
	initArray(a3, *getEntryByType(MET_CE));
	initArray(a4, *getEntryByType(MET_HE3));
	dihedralAngle[23] = computerDihedral(a1,a2,a3,a4);

	return 1;
}

int EnkephalinPDB::writeDihedral( ) {
	if(atomNumberOfTotalPDB < 1) return 0;
	
	cout << "OUTPUT FORMAT: ECEPPAK" << endl;
	cout << "Energy: "  << setprecision(3)  << getEnergy( ) << endl;
	for(int i=0; i<24; i++) {
		cout<< setw(8) 
			<< setprecision(3) 
			<< setiosflags(ios::showpoint) 
			<< setiosflags(ios::fixed)
			<< dihedralAngle[i];
		if (i==5 || i==8 || i==11 || i== 16 || i==23) cout << endl;
	}

	cout.flush();

	return 1;
}
