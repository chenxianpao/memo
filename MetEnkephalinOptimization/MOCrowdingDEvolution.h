#if !defined(TMAN_MULTIMODAL_OPTIMIZATION_CROWDINGDE_SL)
#define TMAN_MULTIMODAL_OPTIMIZATION_CROWDINGDE_SL

#include <ga/GABaseGA.h>
#include <ga/gaRealGenome.h>
#include <ga/gatypes.h>

class MOCrowdingSLDEvolution : public GAGeneticAlgorithm {
public:
	GADefineIdentity("CrowdingDE-SL", 362);
public:
	static GABoolean TerminateUponMaxMin(GAGeneticAlgorithm&);
	static GABoolean TerminateUponFunctionEvaluation(GAGeneticAlgorithm& ga);
public:
	MOCrowdingSLDEvolution(const GAGenome& g);
	virtual ~MOCrowdingSLDEvolution();
public:
	virtual void initialize(unsigned int seed=0);
	virtual void step();
	MOCrowdingSLDEvolution & operator++() { step(); return *this; }
public:
	virtual bool IsBetter(const GAGenome& first, const GAGenome& second);
	float distance(const GAGenome& first, const GAGenome& second);
	GAGenome* findNearestIndividual(const GAGenome& trial, const GAPopulation& population);

public:
	void quickSortAscendingFitness (GAGenome **c, int l, int r);
	void quickSortDescendingFitness(GAGenome **c, int l, int r);

	GAGenome** extractGenomeArrayList(GAGenome* parent, int pos, GAPopulation* pop);
	double* getPrepareArray(GAGenome** genomeArrayList, int size);
	int selectIndex(double* psum, int size);
	int removeGenomeFromArray(GAGenome**list, int lstSize, int pos);
	void MOCrowdingSLDEvolution::outPopulationInformation(GAPopulation* pop);





	



	GAGenome** getGenomeArrayList(GAGenome* parent, GAPopulation *pop);
	double* prepareArray(GAGenome** genomeArrayList, int size);
	int removeGenomeList(GAGenome** genomeList, int size, GAGenome* removedGenome);



	GAGenome* select(GAGenome* parent, GAPopulation *pop);
	GAGenome* select(GAGenome** list, int size);
public:
	float	F;
	float	CR;


};

#endif 
