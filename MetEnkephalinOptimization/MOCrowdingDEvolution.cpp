// MOCrowdingDEvolution.cpp: implementation of the MOCrowdingDEvolution class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MOCrowdingDEvolution.h"
#include "math.h"

GABoolean MOCrowdingSLDEvolution::TerminateUponMaxMin(GAGeneticAlgorithm& ga) 
{
	MOCrowdingSLDEvolution& mode = dynamic_cast<MOCrowdingSLDEvolution&>(ga);
	double fmax = mode.population().max();
	double fmin = mode.population().min();
	if (fmax-fmin < 0.000001) 
		return gaTrue;     //gaTrue=1;
	else 
		return gaFalse;    //gaFalse=0;
}
GABoolean MOCrowdingSLDEvolution::TerminateUponFunctionEvaluation(GAGeneticAlgorithm& ga)
{
	MOCrowdingSLDEvolution& mode = dynamic_cast<MOCrowdingSLDEvolution&>(ga);
	if(mode.statistics().numeval<10000) 
		return gaFalse; 
	else 
		return gaTrue;
}

bool MOCrowdingSLDEvolution::IsBetter(const GAGenome& first, const GAGenome& second) 
{
	bool flag;
	if (minimaxi() == MAXIMIZE)     //MAXIMIZE=1;
	{   
		if (first.score() > second.score()) 
			flag = true; 
		else 
			flag = false;
	}
	else{
		if (first.score() < second.score()) 
			flag = true; 
		else 
			flag = false;
	}
	return flag;
}

void MOCrowdingSLDEvolution::outPopulationInformation(GAPopulation* pop){
	fstream outfile("info_sl.txt", ios::out);
	
	int popSize = pop->size();
	for(int i=0; i<popSize; i++) {
		GAGenome* g = &(pop->individual(i));
		GARealGenome * rg = dynamic_cast<GARealGenome*>(g);
		int dim = rg->size();
		for(int j=0; j<dim; j++) {
			outfile << rg->gene(j) << "\t";
		}
		outfile << endl;
	}
	outfile.close();
}

float MOCrowdingSLDEvolution::distance(const GAGenome& first, const GAGenome& second)
 {
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

GAGenome* MOCrowdingSLDEvolution::findNearestIndividual(const GAGenome& trial, const GAPopulation& population){
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


void MOCrowdingSLDEvolution::quickSortAscendingFitness(GAGenome **c, int l, int r) 
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


void MOCrowdingSLDEvolution::quickSortDescendingFitness(GAGenome **c, int l, int r) 
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


GAGenome** MOCrowdingSLDEvolution::extractGenomeArrayList(GAGenome* parent, int pos, GAPopulation* pop)
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
double* MOCrowdingSLDEvolution::getPrepareArray(GAGenome** genomeArrayList, int size) {
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

GAGenome** MOCrowdingSLDEvolution::getGenomeArrayList(GAGenome* parent, GAPopulation *pop)
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

double* MOCrowdingSLDEvolution::prepareArray(GAGenome** genomeArrayList, int size) 
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



GAGenome* MOCrowdingSLDEvolution::select(GAGenome* parent, GAPopulation *pop) 
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


GAGenome* MOCrowdingSLDEvolution::select(GAGenome** genomeList, int listSize) 
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

int MOCrowdingSLDEvolution::selectIndex(double* psum, int size) 
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
int MOCrowdingSLDEvolution::removeGenomeFromArray(GAGenome**list, int lstSize, int pos) 
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


int MOCrowdingSLDEvolution::removeGenomeList(GAGenome** genomeList, int size, GAGenome* removedGenome) 
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


MOCrowdingSLDEvolution::MOCrowdingSLDEvolution(const GAGenome& g) : GAGeneticAlgorithm(g) 
{
	F  = 0.5;	// a good initial choice used by R.Storn
	CR = 0.1;	// a good initial choice used by R.Storn
}

MOCrowdingSLDEvolution::~MOCrowdingSLDEvolution() 
{

}

void MOCrowdingSLDEvolution::initialize(unsigned int seed) 
{
	//if(TCONDITION==0) 
	//	terminator(TerminateUponMaxMin);	         // 设置终止条件
	terminator(TerminateUponFunctionEvaluation);
	GARandomSeed(seed);							// seed = 0,表示不可复现的随机数
	pop->size(100);
	pop->initialize();
	pop->evaluate(gaTrue);
	stats.reset(*pop);
}

void MOCrowdingSLDEvolution::step() {
	int popSize = pop->size();
	int genomeLength= (dynamic_cast<GARealGenome*>(&(pop->individual(0))))->length( );

	for(int i=0; i<popSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(target->clone());
		
		//selection
		stats.numsel +=3;
		GAGenome *a1, *b1, *c1;

		GAGenome** genomeList = 0;
		double *psum = 0;
		int newSize,removedIndex;
		//generate three target individual;
		genomeList = extractGenomeArrayList(target, i, pop);
		psum = getPrepareArray(genomeList, popSize-1);
		removedIndex = selectIndex(psum, popSize-1);
		a1 = genomeList[removedIndex];
		newSize = removeGenomeFromArray(genomeList, popSize-1, removedIndex);
		delete [] psum;

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
			float newgene = c->gene(k) + F*(a->gene(k)-b->gene(k));
			float gLower = trial->alleleset(k).lower();
			float gUpper = trial->alleleset(k).upper();
			if (newgene<gLower || newgene>gUpper) {
				newgene = GARandomFloat(gLower,gUpper);
			}
			trial->gene(k, newgene);
		}
		stats.nummut++;
		
		// crossover
		int j = GARandomInt(0, genomeLength-1);
		for(k=0; k<genomeLength; k++) {
			if ( (GARandomFloat() <= CR) || (k==j) ) {
				stats.numcro++;
			}
			else {
				trial->gene(k, target->gene(k));
			}
		}
		stats.numcro++;
		
		// acceptance
		stats.numeval++; // trial evaluation
		GAGenome* nearestIndividual = findNearestIndividual(*trial, *pop);
		
		if (IsBetter(*trial, *nearestIndividual)) {
			nearestIndividual->copy(*trial);
			stats.numrep++;
		}
	}
	outPopulationInformation(pop);
	pop->sort(gaTrue);
	pop->statistics(gaTrue);
	stats.update(*pop);
}

