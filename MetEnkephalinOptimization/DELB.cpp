///////////////////////////////////////////////////////////////////////////////////////////////////
// DELB.cpp: implementation of the CDELB class.
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "DELB.h"
#include <ga/garandom.h>
#include <math.h>
#include "config.h"




GABoolean CDELB::TerminateUponMaxMin(GAGeneticAlgorithm& ga) {
	float difference = ga.population().max() - ga.population().min();
	if (difference < EPSILON) 
		return gaTrue; 
	else 
		return gaFalse;
}

CDELB::CDELB(const GAGenome& c) : GAGeneticAlgorithm(c) {
	trialPopulation= 0; 
	mutateFactor = 0.5;		// a good initial choice used by R.Storn
	crossFactor  = 0.1;		// a good initial choice used by R.Storn
	simplexFactor= 0.5;

	energyEvals  = 0.0; 
	energyEvalsSec=0.0;

	energyEvals				= 0.0;
	deEnergyEvals			= 0.0;	
	reflectEnergyEvals		= 0.0;
	contractEnergyEvals		= 0.0;
	energyEvalsSec			= 0.0;
	deEnergyEvalsSec		= 0.0;
	reflectEnergyEvalsSec	= 0.0;
	contractEnergyEvalsSec	= 0.0;


	deReplace				= 0;
	reflectReplace			= 0;
	contractReplace			= 0;

}

CDELB::~CDELB() {
	delete trialPopulation;
}


bool CDELB::better(const GAGenome& firstGenome, const GAGenome& secondGenome) {
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


bool CDELB::legal(const GAGenome& genome) {
	const GARealGenome& realGenome= dynamic_cast<const GARealGenome&>(genome);
	int genomeLength = realGenome.length();

	for(int i=0; i<genomeLength; i++){
		if ( realGenome.gene(i) < realGenome.alleleset(i).lower() || 
			 realGenome.gene(i) > realGenome.alleleset(i).upper() )
			return false;
	}
	return true;
}

void CDELB::reflectGenome(GAGenome& trialGenome, GAGenome& bestGenome, GAGenome& resultGenome) {
	GARealGenome& trialRealGenome  = dynamic_cast<GARealGenome&>(trialGenome);
	GARealGenome& bestRealGenome   = dynamic_cast<GARealGenome&>(bestGenome);
	GARealGenome& resultRealGenome = dynamic_cast<GARealGenome&>(resultGenome);
	int genomeLength = trialRealGenome.length();

	int tryNumber = 0;
	do {
		float coe = GARandomFloat(-1,1);
		for(int iRef=0; iRef<genomeLength; iRef++) {
			float geneLower = resultRealGenome.alleleset(iRef).lower();
			float geneUpper = resultRealGenome.alleleset(iRef).upper();
			if(geneUpper-geneLower>0.01) {
				float gvalue = coe*trialRealGenome.gene(iRef)+(1-coe)*bestRealGenome.gene(iRef);  
				resultRealGenome.gene(iRef, gvalue);
			}
		}
		tryNumber++;
	}while(!legal(resultRealGenome) || tryNumber<10);

	if(!legal(resultRealGenome)) {
		for(int iRef=0; iRef<genomeLength; iRef++) {
			float geneLower = resultRealGenome.alleleset(iRef).lower();
			float geneUpper = resultRealGenome.alleleset(iRef).upper();

			if ( resultRealGenome.gene(iRef) < geneLower|| resultRealGenome.gene(iRef) > geneUpper) {
				float newGene = GARandomFloat(geneLower, geneUpper);
				resultRealGenome.gene(iRef, newGene);
			}
		}
		//resultGenome.copy(trialGenome); 
	}
}


void CDELB::contractGenome(GAGenome& trialGenome, GAGenome& bestGenome, GAGenome& resultGenome) {
	GARealGenome& trialRealGenome  = dynamic_cast<GARealGenome&>(trialGenome);
	GARealGenome& bestRealGenome   = dynamic_cast<GARealGenome&>(bestGenome);
	GARealGenome& resultRealGenome = dynamic_cast<GARealGenome&>(resultGenome);
	int genomeLength = trialRealGenome.length();

	int tryNumber = 0;
	do {
		float coe = GARandomFloat(-1,1);
		for(int iRef=0; iRef<genomeLength; iRef++){
			float geneLower = resultRealGenome.alleleset(iRef).lower();
			float geneUpper = resultRealGenome.alleleset(iRef).upper();
			if(geneUpper-geneLower>0.01) {
				float gvalue = coe*bestRealGenome.gene(iRef)+(1-coe)*trialRealGenome.gene(iRef);  
				resultRealGenome.gene(iRef, gvalue);
			}
			tryNumber++;
		}
	}while(!legal(resultRealGenome)|| tryNumber<10);

	if(!legal(resultRealGenome)) {
		for(int iRef=0; iRef<genomeLength; iRef++) {
			float geneLower = resultRealGenome.alleleset(iRef).lower();
			float geneUpper = resultRealGenome.alleleset(iRef).upper();
		
			if ( resultRealGenome.gene(iRef) < geneLower|| resultRealGenome.gene(iRef) > geneUpper) {
				float newGene = GARandomFloat(geneLower, geneUpper);
				resultRealGenome.gene(iRef, newGene);
			}
		}
		//resultGenome.copy(trialGenome); 
	}
}


void CDELB::minimizeGenome(GAGenome& genome, double& evals, double& evalsPerSec){
	GARealGenome& realGenome = dynamic_cast<GARealGenome&>(genome);
	int genomeLength = realGenome.length();

	double *geomData = new double[genomeLength];
	for(int i=0; i<genomeLength; i++) {
		geomData[i] = realGenome.gene(i);
	}


	writeMetEnkephalinEnergyMinimizationFile("metEnkephalinEnergyMinimizationSourceFile.inp",geomData);
	system("i.bat MINIMIZE metEnkephalinEnergyMinimizationSourceFile MET_ENERGY_MIN rr rr  1>> log 2> err");

	double genomeScore = readDataFromFile("outo.MET_ENERGY_MIN", geomData);
	readPerformanceDataFromFile("main_out.MET_ENERGY_MIN", evals, evalsPerSec);

	for(i=0; i<genomeLength; i++) {
		realGenome.gene(i, (float)geomData[i]);
	}
	realGenome.score((float)genomeScore);


	system("del outo.MET_ENERGY_MIN");
	system("del log");
	delete [] geomData;
}

void CDELB::minimizePopulation(GAPopulation* pop1, double &energyEvals, double& energyEvalsPerSec) {
	writeMetEnkephalinEnergyMultipleMinimizationFile("metEnkephalinEnergyMultipleMinimizationSourceFile.inp");
	writeMetEnkephalinEnergyMultipleConformationInputFile("outo.multipleInput", pop1);
	system("i.bat MINIMIZE metEnkephalinEnergyMultipleMinimizationSourceFile MET_MULTI_MIN multipleInput rr  1>> log 2> err");
	readMetEnkephalinMultipleMinimizationFile("outo.MET_MULTI_MIN", pop1);

	//int i = system("del log");
	//int j = system("del outo.multipleInput");
	int k = system("del outo.MET_MULTI_MIN");
	//int m = system("del main_out.MET_MULTI_MIN");
	
	//Sleep (4000);
}


void CDELB::writeMetEnkephalinEnergyMultipleMinimizationFile(const char *fileName) {
	FILE *f;
	f=fopen(fileName,"wt");

	fprintf(f,"$CNTRL\n");
	fprintf(f,"runtyp = minimize\n");
	fprintf(f,"res_code= one_letter\n");
	fprintf(f,"$END\n\n");

	fprintf(f,"$SEQ\n");
	fprintf(f,"H\n");
	fprintf(f,"YGGFM\n");
	fprintf(f,"O\n");
	fprintf(f,"$END\n\n");

	fprintf(f,"$ENERCALC\n");
	fprintf(f,"READ_CONF\n");
	fprintf(f,"$END\n\n");

	fprintf(f,"$GEOM\n\n\n\n\n\n\n\n");
	fprintf(f,"$END\n");

	fclose(f);
}

void CDELB::writeMetEnkephalinEnergyMultipleConformationInputFile(const char *fileName, GAPopulation* pop1){
	FILE *f;
	f=fopen(fileName,"wt");

	GARealGenome& realGenome = dynamic_cast<GARealGenome&>(pop1->individual(0));
	int genomeLength = realGenome.length();
	double *geomData = new double[genomeLength];

	for(int i=0; i<pop1->size(); i++) {
		GARealGenome& realGenomei = dynamic_cast<GARealGenome&>(pop1->individual(i));
		for(int j=0; j<genomeLength; j++) {
			geomData[j] = realGenomei.gene(j);
		}

		fprintf(f,"%d", i+1);
		fprintf(f,"\n");
		fprintf(f,"   1  20   6   6   5  11  11");
		fprintf(f,"\n\n");

		for(int k=0; k<6;k++) {
			fprintf(f,"%8.3f",geomData[k]);
		}
		fprintf(f,"\n");
		for(k=6; k<9;k++) {
			fprintf(f,"%8.3f",geomData[k]);
		}
		fprintf(f,"\n");
		for(k=9; k<12;k++) {
			fprintf(f,"%8.3f",geomData[k]);
		}
		fprintf(f,"\n");
		for(k=12; k<17;k++) {
			fprintf(f,"%8.3f",geomData[k]);
		}
		fprintf(f,"\n");
		for(k=17; k<24;k++) {
			fprintf(f,"%8.3f",geomData[k]);
		}
		fprintf(f,"\n\n");
	}

	delete [] geomData;
	fclose(f);
}

void CDELB::readMetEnkephalinMultipleMinimizationFile(const char *fileName, GAPopulation* pop1) {
	ifstream fin(fileName); 

	for(int i=0; i<pop->size(); i++) {
		GARealGenome& realGenome = dynamic_cast<GARealGenome&>(pop->individual(i));

		const int   LINE_LENGTH = 1000; 
		char  str[LINE_LENGTH]; 
		int   line  = 0;
		int   index = 0;
		
		while( fin.getline(str,LINE_LENGTH)) {  
			line++;
			if(line==1) {
				double ETOT = readBufToETOT(str, LINE_LENGTH);
				realGenome.score((float)ETOT);
			}
			if(line==4) {
				char substr[9]={0};
				for(int m=0; m<6;m++) {
					for(int j=0; j<8;j++) {
						substr[j] = str[m*8+j];
					}
					realGenome.gene(index++,(float)atof(substr));
				}
			}
			if(line==5) {
				char substr[9]={0};
				for(int m=0; m<3;m++) {
					for(int j=0; j<8;j++) {
						substr[j] = str[m*8+j];
					}
					realGenome.gene(index++,(float)atof(substr));
				}
			}

			if(line==6) {
				char substr[9]={0};
					for(int m=0; m<3;m++) {
						for(int j=0; j<8;j++) {
							substr[j] = str[m*8+j];
					}
					realGenome.gene(index++,(float)atof(substr));
				}
			}
		
			if(line==7) {
				char substr[9]={0};
				for(int m=0; m<5;m++) {
					for(int j=0; j<8;j++) {
						substr[j] = str[m*8+j];
					}
					realGenome.gene(index++,(float)atof(substr));
				}
			}

			if(line==8) {
				char substr[9]={0};
				for(int m=0; m<7;m++) {
					for(int j=0; j<8;j++) {
						substr[j] = str[m*8+j];
					}
					realGenome.gene(index++,(float)atof(substr));
				}
			}
			if(line==9) {
				break;
			}
		}
	}

	fin.close();
}

double CDELB::readBufToETOT(char* str, int length){
	double ret		= 0.0;
	char buf[50]	= {0};
	bool first		= true;
	bool second		= false;

	int index=0;
	for(int i=0; i<length; i++) {
		if(str[i]!=' ' && first==true) {
			if(str[i+1]==' ') {first=false; second=true;}
		} else {
			if(second && str[i] != ' ') {
				buf[index++] = str[i];
				if(str[i+1] == ' ') break;
			}
		}
	}
	
	ret =  atof(buf);
	return ret;
}



void CDELB::writeMetEnkephalinEnergyMinimizationFile(const char* filename, double* geomData) {
	FILE *f;
	f=fopen(filename,"wt");

	fprintf(f,"$CNTRL\n");
	fprintf(f,"runtyp = minimize\n");
	fprintf(f,"!PRINT_CART\n");
	fprintf(f,"!OUTFORMAT =PDB\n");
	fprintf(f,"!FILE  = metEnkephalin\n");
	fprintf(f,"res_code= one_letter\n");
	fprintf(f,"$END\n\n");

	fprintf(f,"$SEQ\n");
	fprintf(f,"H\n");
	fprintf(f,"YGGFM\n");
	fprintf(f,"O\n");
	fprintf(f,"$END\n\n");

	fprintf(f,"$GEOM\n\n");
	for(int i=0; i<6;i++) {
		fprintf(f,"%8.3f",geomData[i]);
	}
	fprintf(f,"\n");
	for(i=6; i<9;i++) {
		fprintf(f,"%8.3f",geomData[i]);
	}
	fprintf(f,"\n");
	for(i=9; i<12;i++) {
		fprintf(f,"%8.3f",geomData[i]);
	}
	fprintf(f,"\n");
	for(i=12; i<17;i++) {
		fprintf(f,"%8.3f",geomData[i]);
	}
	fprintf(f,"\n");

	for(i=17; i<24;i++) {
		fprintf(f,"%8.3f",geomData[i]);
	}
	fprintf(f,"\n");
	fprintf(f,"\n");
	fprintf(f,"$END\n");

	fclose(f);
}

double CDELB::readDataFromFile(const char* filename, double *geomData) {
	double ETOT = 100000;

	int arrayIndex =0;
	int i=0;
	int ii = 0;

	ifstream fin(filename); 
    const int LINE_LENGTH = 1000; 
    char str[LINE_LENGTH];  

	char buf[50];
	do{
		fin		>> buf;
		ii++;
		if(ii==2) {ETOT= atof(buf);}
		if(ii==5) {break;}
	}while(!fin.eof());

    while( fin.getline(str,LINE_LENGTH)) {  
		i++;
		if(i==4) {
			char substr[9]={0};
			for(int m=0; m<6;m++) {
				for(int j=0; j<8;j++) {
					substr[j] = str[m*8+j];
				}
				geomData[arrayIndex++] = atof(substr);
			}
		}

		if(i==5) {
			char substr[9]={0};
			for(int m=0; m<3;m++) {
				for(int j=0; j<8;j++) {
					substr[j] = str[m*8+j];
				}
				geomData[arrayIndex++] = atof(substr);
			}
		}

		if(i==6) {
			char substr[9]={0};
			for(int m=0; m<3;m++) {
				for(int j=0; j<8;j++) {
					substr[j] = str[m*8+j];
				}
				geomData[arrayIndex++] = atof(substr);
			}
		}
		
		if(i==7) {
			char substr[9]={0};
			for(int m=0; m<5;m++) {
				for(int j=0; j<8;j++) {
					substr[j] = str[m*8+j];
				}
				geomData[arrayIndex++] = atof(substr);
			}
		}

		if(i==8) {
			char substr[9]={0};
			for(int m=0; m<7;m++) {
				for(int j=0; j<8;j++) {
					substr[j] = str[m*8+j];
				}
				geomData[arrayIndex++] = atof(substr);
			}
		}
    }

	fin.close();

	return ETOT;
}

void CDELB::readPerformanceDataFromFile(const char* filename, double& evals, double& evalsPerSec){
	ifstream fin(filename); 
    const int LINE_LENGTH = 1000; 
    char str[LINE_LENGTH]={0};  

	int  line = -1;
	bool find = false;

	while(fin.getline(str,LINE_LENGTH) ){  
		//char *s1=" Total Number of Energy Evals.";
		//char *s2="Energy Evals. per sec =  ";
		if(find) {
			line++;
			if(line==2) {
				int  index = 24;
				char str1[20]={0};
				for(int i=0; i<20; i++) {
					str1[i] = str[index++];
				}
				evalsPerSec = atof(str1);
				//cout << str << "\t" << str1 << "\t" << evalsPerSec<<  endl;
				break;
			}
		}
		else {
			if(str[0] == ' ' && str[1] =='T' && str[2] =='o' && str[3] =='t' && str[4] =='a' && str[5]=='l' &&
				str[24]== 'E' && str[25]=='v' && str[26]=='a' && str[27]=='l' && str[28]=='s' && str[29]=='.') {
				int  index   = 32;
				char str1[20]= {0};
				for(int i=0; i<20; i++) {
					str1[i] = str[index++];
				}
				evals = atof(str1);
				line = 0; find = true;
				//cout << str << "\t" << str1 << "\t" << evals<<  endl;
			}
		}
	}
	fin.close();
}

void CDELB::persistConformationData(const GAPopulation& p, char* mainFileName, int generation) {
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
		outfile << i+1		<< "\t"
			    << fscore	<< "\t";
		for(int j = 0; j<trial.length(); j++) {
			outfile << trial.gene(j) << "\t";
		}
		outfile << "\n";
	}
	outfile.close();

	ofstream outSummaryfile(resultSummaryFile, ios::app);

	GARealGenome& bestGenome = dynamic_cast<GARealGenome&>(p.best());
	double bestScore = bestGenome.score();
	outSummaryfile		<< generation		
						<< "\t"
						<< bestScore	
						<< "\t"
						<< deReplace	
						<< "\t"
						<< reflectReplace	
						<< "\t"
						<< contractReplace	
						<< "\t"
						<< deEnergyEvals	
						<< "\t\t\t"
						<< deEnergyEvalsSec	
						<< "\t\t"
						<< reflectEnergyEvals	
						<< "\t\t\t"
						<< reflectEnergyEvalsSec	
						<< "\t\t"
						<< contractEnergyEvals	
						<< "\t\t\t"
						<< contractEnergyEvalsSec	
						<< "\t\t"
						<< energyEvals	
						<< "\t\t\t"
						<< energyEvalsSec	
						<< "\t\t";
	for(int j = 0; j<bestGenome.length(); j++) {
		outSummaryfile << bestGenome.gene(j) << "\t";
	}
	outSummaryfile << "\n";
	outSummaryfile.close();
}

void CDELB::initialize(unsigned int seed) {
	GARandomSeed(seed);
	pop->initialize( );

	double evals = 0;
	double evalsPerSec = 0;

	/*
	minimizePopulation(pop, evals, evalsPerSec);
	pop->customEvaluate();
	for(int i=0; i<pop->size(); i++) {
		GAGenome& genome = pop->individual(i);
		genome.customEvaluate();
	}
	*/


	for(int i=0; i<pop->size(); i++) {
		minimizeGenome(pop->individual(i), evals, evalsPerSec);
		energyEvals += evals;
		deEnergyEvals += evals;
		if (evalsPerSec>0.00001) {
			energyEvalsSec += evals/ evalsPerSec; deEnergyEvalsSec += evals/ evalsPerSec;
		}
		evals = evalsPerSec = 0.0;
	}


	pop->evaluate(gaFalse);

	if(trialPopulation) delete trialPopulation; trialPopulation = pop->clone( );
	stats.reset(*pop);

	if(TCONDITION==0) terminator(TerminateUponMaxMin);
}


///////////////////////////////////////////////////////////////////////////////////////////////////
/*
 1. mutateFactor: in this implementation mutator=0.5; maybe use strategy as follows: 
		if (GAFlipCoin(0.5)) mutateFactor = GARandomFloat(0.4,1);
		else mutateFactor = GARandomFloat(-1,-0.4);
*/
///////////////////////////////////////////////////////////////////////////////////////////////////
void CDELB::step() {
	double evals = 0;
	double evalsPerSec = 0;

	int	populationSize = pop->size( );
	int	genomeLength = (dynamic_cast<GARealGenome*>(&(pop->individual(0))))->length();

	GAPopulation *targetPopulation = pop->clone();			             /*temp target population*/
	GARealGenome& bestGenome = dynamic_cast<GARealGenome&>(targetPopulation->best());
	
	for (int i=0; i<populationSize; i++) {
		GARealGenome *target = dynamic_cast<GARealGenome*>(&(targetPopulation->individual(i)));
		GARealGenome *trial  = dynamic_cast<GARealGenome*>(&(trialPopulation->individual(i)));
		GARealGenome *popi	 = dynamic_cast<GARealGenome*>(&(pop->individual(i)));
	
		/******************************************* Select **************************************/
		GARealGenome *a, *b, *c;                     /*uniform random points in target population*/
		int	ai, bi, ci;               /* the indexes of three select points in target population */
		do ai = GARandomInt(0, populationSize-1); while (ai==i);
		do bi = GARandomInt(0, populationSize-1); while (bi==i || bi==ai);
		do ci = GARandomInt(0, populationSize-1); while (ci==i || ci==bi || ci==ai);
		a = dynamic_cast<GARealGenome*>(&(targetPopulation->individual(ai)));
		b = dynamic_cast<GARealGenome*>(&(targetPopulation->individual(bi)));
		c = dynamic_cast<GARealGenome*>(&(targetPopulation->individual(ci)));
		stats.numsel += 3;

		/******************************************** Mutate *************************************/
		for(int k=0; k<genomeLength; k++) {
			float newGene	= c->gene(k) + mutateFactor*(a->gene(k)-b->gene(k));
			float geneLower	= trial->alleleset(k).lower();
			float geneUpper	= trial->alleleset(k).upper();
			if (newGene < geneLower || newGene > geneUpper) 
				newGene = GARandomFloat(geneLower, geneUpper);    /*customed for specific problem*/
			trial->gene(k, newGene);
		}
		stats.nummut++;
		/******************************************* Crossover ***********************************/
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
		
		
		/******************************************* Acceptance **********************************/
		stats.numeval++;                                                         //trial evaluation
		minimizeGenome(*trial,evals, evalsPerSec); 
		energyEvals += evals; 
		deEnergyEvals += evals;
		if (evalsPerSec>0.00001) {
			energyEvalsSec += evals/evalsPerSec;
			deEnergyEvalsSec += evals/evalsPerSec;
		}
		evals = 0; evalsPerSec = 0;

		if (better(*target, *trial)) {	                          
			trial->copy(*target);	
			deReplace++;
		}else {
			stats.numrep++;
		} 


		/*
		else
		{
			if ((GARandomFloat()<simplexFactor) && (better(bestGenome, *trial))) {
				GARealGenome *reflectionGenome = dynamic_cast<GARealGenome*>(target->clone());
				reflectGenome(*trial, bestGenome, *reflectionGenome);                  //Reflection
				stats.numeval++;	                                        //Reflection evaluation
				minimizeGenome(*reflectionGenome,evals, evalsPerSec);
				energyEvals += evals; 
				reflectEnergyEvals += evals;
				if (evalsPerSec>0.00001) {
					energyEvalsSec += evals/evalsPerSec;
					reflectEnergyEvalsSec += evals/evalsPerSec;
				}
				evals = 0; evalsPerSec = 0;
				if( better(*reflectionGenome, *trial) ) {
					trial->copy(*reflectionGenome);
					reflectReplace++;
				}
				else {				                                               
					GARealGenome *contractionGenome=dynamic_cast<GARealGenome*>(target->clone());
					contractGenome(*trial, bestGenome, *contractionGenome);           //Contraction
					stats.numeval++;                                       //Contraction evaluation
					minimizeGenome(*contractionGenome,evals, evalsPerSec);
					energyEvals += evals; 
					contractEnergyEvals += evals;
					if (evalsPerSec>0.00001) {
						energyEvalsSec += evals/ evalsPerSec; 
						contractEnergyEvalsSec += evals/ evalsPerSec; 
					}
					evals = 0; evalsPerSec = 0;
					if( better(*contractionGenome, *trial) ) {
						trial->copy(*contractionGenome);
						contractReplace++;
					}
					delete contractionGenome;
				}
				delete reflectionGenome;
			}
			stats.numrep++;
		} 
		/********************************************* update ************************************/
		popi->copy(*trial);
		/********************************************* finish ************************************/
	}

	/*
	minimizePopulation(trialPopulation, evals, evalsPerSec);
	trialPopulation->customEvaluate();

	for(int k=0; k<pop->size();k++) {
		GARealGenome& genome1 = dynamic_cast<GARealGenome&>(trialPopulation->individual(k));
		GARealGenome& genome2 = dynamic_cast<GARealGenome&>(pop->individual(k));
		genome1.customEvaluate();
		genome2.customEvaluate();
		if(better(genome1, genome2)) {
			genome2.copy(genome1);
			stats.numrep++;
		} else {
			genome1.copy(genome2);
		}
	}
	*/

	delete targetPopulation;
	
	pop->sort(gaTrue);
	pop->statistics(gaTrue);
	stats.update(*pop);

	persistConformationData(*pop, "data\\ConformationData", generation());
}






