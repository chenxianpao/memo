#if !defined(MET_ENKEPHALIN_DE_ZGJ_ZJUT)
#define MET_ENKEPHALIN_DE_ZGJ_ZJUT


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


//#include <ga/GABaseGA.h>
#include <ga/gaRealGenome.h>
#include <ga/gatypes.h>
#include <config.h>
#include <cam.h>
#include <dataStruct.h>
#include <MyEvalData.h>
#include <EnkephalinPDB.h>


class CDifferentialEvolution : public GAGeneticAlgorithm {
public:          
	GADefineIdentity("BasicDE", 1000);
	static GABoolean TerminateUponMaxMin(GAGeneticAlgorithm& algorithm);

public:
	CDifferentialEvolution(const GAGenome&);
	virtual ~CDifferentialEvolution();

public:
	bool	better(const GAGenome& first, const GAGenome& second);
	void	persistConformationData(const GAPopulation& p, char* mainFileName, int generation);
	void	persistAllEvolutionData(const GAPopulation& p, char* allPopulationFileName);
	void	mergeConformationData(const char *mainFileName, int fileSize, int fileNumber);

public:
	double	singleEnergyEvaluate(double *pGEOM, int length);  //
	void	singleEnergyEvaluate(GAGenome& genome);
	double	singleEnergyMinimize(double *pGEOM, int length);
	double	singleEnergyMinimize(double *pGEOM, bool *spec, int length);
	double	singleEnergyMinimizeGeneUnchanged(double* pGEOM, int length);
	void	singleEnergyMinimize(GAGenome& genome);
	void	singleEnergyMinimizeGeneUnchanged(GAGenome& genome);

	void	writeSingleEnergyEvaluateSourceFile(const char*sourceFile, double *pGEOM, int length);
	void	writeSingleEnergyMinimizeSourceFile(const char* sourceFile, GAGenome& genome);
	void	writeSingleEnergyMinimizeSourceFile(const char* sourceFile, double* pGEOM, int length);
	void	writeSingleEnergyMinimizeSpecVarAngles(const char* sourceFile, bool * spec, double* pGEOM, int length);
	double	readSingleEnergyEvaluateMainOutResultFile(const char* resultFile);
	void	readSingleEnergyMinimizeResultOUTOFile(const char* resultFile, GAGenome& genome);
	void	readSingleEnergyMinimizeResultOUTOFileGeneUnchanged(const char* resultFile, GAGenome& genome);
	double	readSingleEnergyMinimizeResultETOT(const char* resultFile, double*pGEOM, int length);
	double	readSingleEnergyMinimizeResultETOTGeneUnchanged(const char* resultFile);
	void	readSingleEnergyMinimizeResultEvalsFile(const char* mainOutFile);
	void	clearSingleEnergyEvluateResultFiles();
	void	clearSingleEnergyMinimizeResultFiles();

	void	multipleEnergyEvaluate(GAPopulation& p);
	void	multipleEnergyEvaluateWithOneResidue(GAPopulation& p, char* residueName);
	void	multipleEnergyEvaluateWritePDB(GAPopulation& p, char* fileName);
	void	multipleEnergyEvaluate(double** arrayGeom, int arraySize, int dihedralAngleSize);

	void	multipleEnergyMinimize(GAPopulation& p);
	void	multipleEnergyMinimizeWithOneResidue(GAPopulation& p, char* residueName);

	void	multipleEnergyMinimizeWritePDB(GAPopulation& p, char* fileName);
	void	multipleEnergyMinimize(double** arrayGeom, int arraySize, int dihedralAngleSize);
	void	multipleEnergyMinimize(double** arrayGeom, bool *spec, int arraySize, int dihedralAngleSize);
	void	multipleEnergyMinimizeFixedGene(GAPopulation& p);

	void	writeMetEnkephalinMultipleEvaluateSourceFile(const char* sourceFile);
	void	writeMetEnkephalinMultipleEvaluateSourceFile_WritePDB(const char* sourceFile, char* pdbFileName);
	void	writeMetEnkephalinMultipleEvaluateSourceFile(const char* sourceFile, char* residueName);
	void	writeMetEnkephalinMultipleMinimizeSourceFile(const char* sourceFile);
	void	writeMetEnkephalinMultipleMinimizeSourceFile(const char* sourceFile, bool *spec);
	void	writeMetEnkephalinMultipleMinimizeSourceFile(const char* sourceFile, char* residueName);

	void	writeMetEnkephalinMultipleMinimizeSourceFile_WritePDB(const char* sourceFile, char* pdbFileName);
	void	writeMetEnkephalinMultipleInputFile(const char* inputFile, GAPopulation& p);
	void	writeMetEnkephalinMultipleInputFile(const char* inputFile, GAPopulation& p, char* residueName);
	void	writeMetEnkephalinMultipleInputFile(const char* inputFile, double ** arrayGeom, int arraySize, int variableSize);
	void	readMetEnkephalinMultipleResultOUTOFile(const char *resultFile, GAPopulation& p);
	void	readMetEnkephalinMultipleResultOUTOFile(const char *resultFile, GAPopulation& p, char *residueName);
	void	readMetEnkephalinMultipleResultOUTOFile(const char *resultFile, double **arrayGenome, int arraySize, int variableSize);
	void	readMetEnkephalinMultipleResultOUTOFileGeneUnchanged(const char *resultFile, GAPopulation& p);
	void	readMetEnkephalinEnergyEvals(const char* mainOutFile);
	void	clearEnkephalinMultipleMinimizeResultFile();
	void	clearEnkephalinMultipleEvaluateResultFile();

	void	generatePDB(double *pGEOM, int length);
	void	writeSingleEnergyEvaluatePDBSourceFile(const char*sourceFile, double *pGEOM, int length);

public:
	virtual void initialize(unsigned int seed=0);
	CDifferentialEvolution & operator++() { step(); return *this; }
	
	virtual void objectiveData(const GAEvalData& v);

	virtual void step ( );		// Basic DE without Energy Minimization 
	virtual void step1( );		// Basic DE With Energy Minimization (or Basic DE with Reduced Energy Landscape)
	virtual	void step2( );		// Basic DE with Buildup procedure
	virtual void step3( );
	virtual void step4( );


	//ÅÅ¼·²ßÂÔ
	virtual void step5( );
	float distance(const GAGenome& first, const GAGenome& second);
	GAGenome* findNearestIndividual(const GAGenome& trial, const GAPopulation& population);
    void quickSortAscendingFitness(GAGenome **c, int l, int r) ;
    void quickSortDescendingFitness(GAGenome **c, int l, int r) ;
    GAGenome** extractGenomeArrayList(GAGenome* parent, int pos, GAPopulation* pop);
	double* getPrepareArray(GAGenome** genomeArrayList, int size);
	int selectIndex(double* psum, int size);
	int removeGenomeFromArray(GAGenome**list, int lstSize, int pos);
	void CDifferentialEvolution::outPopulationInformation(GAPopulation* pop);



	GAGenome** getGenomeArrayList(GAGenome* parent, GAPopulation *pop);
	double* prepareArray(GAGenome** genomeArrayList, int size);
	int removeGenomeList(GAGenome** genomeList, int size, GAGenome* removedGenome);

	GAGenome* select(GAGenome* parent, GAPopulation *pop);
	GAGenome* select(GAGenome** list, int size);




	virtual void dihedralAngleTest( );

	virtual void save(GAPopulation* p);
	virtual void save(GAPopulation* p, char* fileName);

	long index;
protected:												
	float	mutatorFactor;
	float	crossFactor;

public:
	double	energyEvals;
	double	energyEvalsCPUTime;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RMSD analysis utility 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
private:
	double calculateRMSD(double **ref_xlist, double** mov_xlist, int n_list);
	double calculateRMSD(char* refPDBFile, char* movPDBFile, EnkephalinPDB::RMSDTYPE rmsdType, int number);





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// auxiliary analysis utility 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
public:
	void randomMinimizeRsmdAnalysis( );				// random minimize (RMSD)
	void randomMinimizeDihedralAnalysis( );			// random minimize (dihedral) 
	void variableMinimizeSmoothAnalysis( );			// smooth minimize (individual)
	void variableMinimizeSmoothAnalysis1( );		// smooth minimize (population)

	virtual void randomMinimize( );

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CAM optimization
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
public:

public:
	double CAMObjective(double x0, int index, double* pGEOM, bool *spec, int length, int type);

public:
	bool illPosed(CL& arg);
	void updateLeaf(CLLeafArray& leafArray, Cl& l);
	double CAMOptimize (double	*inputGEOM, 
						double	inputOjbective, 
						bool	*spec, 
						int		index, 
						int		length, 
						double	*outputGEOM);
	
public:
	void test( );
	void test2();




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Conformational Space Annealing Algorithm 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
public:		// parameters
	long	thresholdN;
	double	thresholdDistanceCutoff;
	double	thresholdDifference;


public:		// variables
	GAGenome* preSeed;
	double numberMinimized;
	double numberN;
	double averageDistanceOfFirstBank;
	double xCutoff;
	double currentDistanceCutoff;

public:
	double dihedralAngleDistance(double i, double j);
	double geometryDistance(double* gi, double* gj, int angleLength);
	double geometryDistance(GAGenome& gi, GAGenome& gj);
	double populationAverageDistance(const GAPopulation& trialPop);

	double computeAnnealingX(double aveDistanceCutoff, int thresholdDistanceCutoffN, double thresholdDistanceCutoff);
	double distanceCutoff(double aveDistanceCutoff, double annealingX, int n);
	void   distanceCuttingCurve(const char* filename, double aveDistanceCutoff, double annealingX);

	SeedStatus getGenomeSeedStatus(GAGenome& genome);
	void   setGenomeSeedStatus(GAGenome& genome, SeedStatus status);
	void   setPopSeedStatus(GAPopulation& triPop, SeedStatus status);

	GAGenome* findNearestCorformations(GAGenome& genome, GAPopulation& comparePop);
	GAGenome* findNearestCorformations(GAGenome& genome, GAPopulation& comparePop, int index);


	GAGenome*	  getSeedConformation(GAPopulation& bank);
	bool		  hasSeedConformation(GAPopulation& bank);
	GAPopulation* getTrialConformations(GAGenome& seedGenome, const GAPopulation& firstBank, const GAPopulation& bank);
	GAPopulation* getTrialConformationsWithoutMinimized(GAGenome& seedGenome, const GAPopulation& firstBank, const GAPopulation& bank);
	void		  updateBank(GAPopulation& trialPop, GAPopulation& bank);
	void		  updateBank2(GAPopulation& trialPop, GAPopulation& bank);


	void   conformationSpaceAnnealingOptimize();
	void   conformationSpaceAnnealingOptimize2();
	void   conformationSpaceAnnealingContinuousOptimize( );

};




#endif 

