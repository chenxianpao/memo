//////////////////////////////////////////////////////////////////////
// DELB.h: interface for the CDELB class.
//////////////////////////////////////////////////////////////////////

#if !defined(MET_ENKEPHALIN_DELB_ZGJ_ZJUT)
#define MET_ENKEPHALIN_DELB_ZGJ_ZJUT

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <ga/GABaseGA.h>
#include <ga/gaRealGenome.h>
#include <ga/gatypes.h>




class CDELB : public GAGeneticAlgorithm  {
public:                                     // define the Identification algorithm
	GADefineIdentity("CDELB", 1100);
	static GABoolean TerminateUponMaxMin(GAGeneticAlgorithm&);

public:
	CDELB(const GAGenome&);
	virtual ~CDELB();

public:
	virtual bool	better(const GAGenome& first, const GAGenome& second);
	virtual bool	legal (const GAGenome& g);
	virtual void	reflectGenome (GAGenome& trialGenome, GAGenome& bestGenome,GAGenome& resultGenome);
	virtual void	contractGenome(GAGenome& trialGenome, GAGenome& bestGenome,GAGenome& resultGenome);
	virtual void	minimizeGenome(GAGenome& genome, double& energyEvals, double& energyEvalsPerSec);

	virtual void    minimizePopulation(GAPopulation* pop1, double &energyEvals, double& energyEvalsPerSec);
	virtual void	writeMetEnkephalinEnergyMinimizationFile(const char* fileName, double* geomData);
	virtual	void	writeMetEnkephalinEnergyMultipleMinimizationFile(const char *fileName);
	virtual	void	writeMetEnkephalinEnergyMultipleConformationInputFile(const char *fileName, GAPopulation* pop1);
	virtual void	readMetEnkephalinMultipleMinimizationFile(const char *fileName, GAPopulation* pop1);
	virtual double	readBufToETOT(char* str, int length);


	virtual double	readDataFromFile(const char* filename, double *geomData);
	virtual void	readPerformanceDataFromFile(const char* filename, double& evals, double& evalsPerSec);
	virtual void	persistConformationData(const GAPopulation& p, char* mainFileName, int generation);

public:
	virtual void initialize(unsigned int seed=0);
	virtual void step();
	CDELB & operator++() { step(); return *this; }

public:
	GAPopulation *trialPopulation;					
	
	float	mutateFactor;								// scaling mutate constant
	float	crossFactor;	                                 // crossover constant
	float	simplexFactor;                    // simplex local optimization cutoff
	
	double	energyEvals;
	double	deEnergyEvals;
	double	reflectEnergyEvals;
	double	contractEnergyEvals;
	double  energyEvalsSec;
	double	deEnergyEvalsSec;
	double	reflectEnergyEvalsSec;
	double	contractEnergyEvalsSec;

	long	deReplace;
	long	reflectReplace;
	long	contractReplace;
};

#endif 
