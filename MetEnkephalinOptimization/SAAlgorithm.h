#if !defined(MET_ENKEPHALIN_CSA_ZGJ_ZJUT)
#define MET_ENKEPHALIN_CSA_ZGJ_ZJUT

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <ga/GABaseGA.h>
#include <ga/gatypes.h>
#include <ga/gaRealGenome.h>
#include <EnergyCalculator.h>
#include <MyEvalData.h>
#include <math.h>
#include <EnkephalinPDB.h>
#include <rmsd.h>


class CSAAlgorithm : public GAGeneticAlgorithm {
public:
	GADefineIdentity("CSAAlgrithm", 1001);
public:
	CSAAlgorithm(const GAGenome&);
	virtual ~CSAAlgorithm();

public:
	virtual void initialize(unsigned int seed=0);
	CSAAlgorithm & operator++() { step(); return *this; }
	virtual void step ( );
public:
	GAPopulation	*firstBank;					// reference bank
	GAPopulation	*bank;						// dynamic bank
	int				nBank;						// bank size;
	double			dCut;						// distance cutoff 
	double			dInitialAve;				// initial average distance
private:
	double			annealingScheduleX;			
	double			energyMinimizedCutOff;
	double			nicheCutoff;
public:	// 性能指标
	double			energyEval;
	double			energyEvalSec;
	double			energyMinimizedValue;		

public:
	GAGenome* selectSeedConformation(GAPopulation* bank);
	GAGenome* selectSeedConformation(GAGenome* previousSeed, GAPopulation *bank);
	GAPopulation* getTrialConformations(GAGenome* seedConformation, GAPopulation* firstBank, GAPopulation* bank);
	void upDateBank(GAPopulation* trialPopulation, GAPopulation* bankPopulation);

private:
	void multipleEnergyMinimize(GAPopulation& p);
	double computeModLeastPositiveValue(double numerator, double denominator) {
		while(numerator<0) {
			numerator += denominator;
		}
		return fmod(numerator, denominator);      //numerator-quotient*denominator
	}
	double computeModLeastPositiveValueVector(double *geom1, double *geom2, int dim){
		double sum = 0.0;
		for(int i=0; i<dim; i++) {
			double r1 = computeModLeastPositiveValue(geom1[i]-geom2[i], 360); //non-symmetry
			double r2 = 360-r1;
			if(r1<r2) sum+=r1; else sum+=r2;
		}
		return sum;
	}
	double computerGenomeDistance(GAGenome* g1, GAGenome* g2);
	double computeRMSDDistance(double *geomAngle1, double*geomAngle2, int angleSize);
	double computerGenomeRMSDDistance(GAGenome* g1, GAGenome* g2);
	double distanceAverage(GAPopulation& p);

	GAGenome* findNearestConformationFromBank(GAGenome* genome, GAPopulation* bank);
	GAGenome* findWorstConformationFromBank(GAPopulation* bank);

private:
	bool distinctBeforePositionL(int*array, int l) {
		for(int i=0; i<l; i++) {
			if(array[i] == array[l]) return false;
		}
		return true;
	}
	
	double getCurrentDcut( ) {
		double retValue; 
		if(energyMinimizedValue<energyMinimizedCutOff) 
			retValue = dInitialAve/2.0 * pow(annealingScheduleX, energyMinimizedValue);
		else
			retValue = nicheCutoff;
		return retValue;
	}
	
	void setGenomeSeedStatus(GAGenome* g, SeedStatus status) {
		MyEvalData* pMyEvalData = dynamic_cast<MyEvalData*>(g->evalData());
		pMyEvalData->setSeedStatus(status);
	}
	
	SeedStatus getGenomeSeedStatus(GAGenome* g) {
		MyEvalData* pMyEvalData = dynamic_cast<MyEvalData*>(g->evalData());
		return pMyEvalData->getSeedStatus();
	}

	void setGenomeSeedStatus(GAPopulation* p, SeedStatus status) {
		for(int i=0; i<p->size(); i++) {
			GAGenome* g = &p->individual(i);
			setGenomeSeedStatus(g, status);
		}
	}

	void resetBank(GAPopulation* pop) {
		for(int i=0; i<pop->size(); i++) {
			setGenomeSeedStatus(&pop->individual(i), DEUnknown);
		}
	}

private:
	void screenPrint(GAGenome* g);
	void screenPrint(GAPopulation* p);
};

#endif 

