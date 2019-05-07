#ifndef	ENKEPHALINPDB_ZGJ_H
#define ENKEPHALINPDB_ZGJ_H

#include "AtomEntry.h"
#include "fstream.h"
#include "iostream.h"
#include "iomanip.h"
#include "math.h"


class EnkephalinPDB {
/*---------------------------------------------------------------------------------*/
public:
	EnkephalinPDB();
	virtual ~EnkephalinPDB();
/*---------------------------------------------------------------------------------*/
public:
	typedef enum {
		AlphaC		= 0,
		Backbone	= 1,
		HeavyAtom	= 2,
		AllAtom		= 3
	} RMSDTYPE;

/*---------------------------------------------------------------------------------*/
public:
	void	setEnergy(double value) { energy = value;}
	double	getEnergy( ) {return energy;}

	void	setPDBFileName(const char* file) {strcpy(PDBFileName, file); }
	char*	getPDBFileName( ) {return PDBFileName;}

	CAtomEntry* getEntryByIndex(int index) {return pAtomEntry[index];}
	CAtomEntry* getEntryByType(ATOMTYPE atomType) {return pAtomEntry[atomType];}

	int getAtomNumber( ) {return atomNumberOfTotalPDB; }
/*---------------------------------------------------------------------------------*/
public:
	int  readPDB(char* PDBFileName);
	int  writePDB( ); 
	double** extractRMSDTypeArray(RMSDTYPE rmsdType, int& numberOfTypeAtom);
	void freeArray(double ** xlist, int number);
public:
	void combineFileName(char* dictory, char* firstName, char*suffixName, int middleIndex, char* fullFileName);

private:
	// the atomEntry is delete by system automatically.
	void add(CAtomEntry* atomEntry) {pAtomEntry[atomNumberOfTotalPDB++] = atomEntry;};
/*---------------------------------------------------------------------------------*/
public:
	int computeDihedral( );		
	int writeDihedral( );

public:
	void initArray(double *a, CAtomEntry& atomEntry){
		a[0] = atomEntry.getX();
		a[1] = atomEntry.getY();
		a[2] = atomEntry.getZ();
	}

	void vec_sub(double *a, const double *b, const double *c) {
		a[0] = b[0]-c[0];
		a[1] = b[1]-c[1];
		a[2] = b[2]-c[2];
	}

	double dot_prod(const double *v1, const double *v2) {
		return v1[0]* v2[0] + v1[1]* v2[1] + v1[2] * v2[2];
	}

	double * cross_prod(double *x1, const double *x2, const double *x3) {
		x1[0] =  x2[1]*x3[2] - x3[1]*x2[2];
		x1[1] = -x2[0]*x3[2] + x3[0]*x2[2];
		x1[2] =  x2[0]*x3[1] - x3[0]*x2[1];

		return x1;
	}
	
	//compute the angle between two vectors a & b (in degrees [0,180]).
	double angle(const double *a, const double *b) {
		double ab[3];
		cross_prod(ab, a, b);
		double psin = sqrt(dot_prod(ab,ab));
		double pcos = dot_prod(a,b);
		return 57.2958f * (double) atan2(psin, pcos);
	}
	//compute the angle formed between three consecutive atoms across two bonds
	//atom a2 joins two other atoms a1, a3. 
	double angle(const double *a1, const double *a2, const double *a3) {
		double v1[3];
		double v2[3];
		vec_sub(v1, a1, a2);
		vec_sub(v2, a3, a2);
		return angle(v1, v2);
	}


	//compute the dihedral angle for the given atoms, returning a value between -180 and 180
	//using faster, cleaner implementation based on atan2
	double computerDihedral(const double *a1, const double *a2, 
						    const double *a3, const double *a4) 
	{
		double r1[3], r2[3], r3[3], n1[3], n2[3];
		
		vec_sub(r1, a2, a1);
		vec_sub(r2, a3, a2);
		vec_sub(r3, a4, a3);

		cross_prod(n1, r1, r2);
		cross_prod(n2, r2, r3);

		double pcos = dot_prod(n1, n2);
		double psin = dot_prod(n1, r3) * sqrt(dot_prod(r2, r2));

		return 57.2958f * (double)atan2(psin, pcos);
	}
/*************************************************************************************/

private:
	CAtomEntry ** pAtomEntry;
	double energy;
	double dihedralAngle[24];
	int atomNumberOfTotalPDB;
	char PDBFileName[200];
};

#endif 