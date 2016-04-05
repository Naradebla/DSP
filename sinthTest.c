#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//400 Hz
//#define FREQUNECY 400 
//20 Ksample/sec
#define SAMPLERATE 20000
//2 V
#define AMPLITUDE 2 
static inline int periodSample(double sRate,double freq){
	double samples = (1/freq)*sRate;
	return samples;
	//return 50;
}
double* sineWave(int sRate, int freq,int ampl){
	//sRate-> sample rate (sample/sec)
	//freq-> sine wave frequency (Hz)
	//ampl-> wave amplitude (V?)
	
	double tSample = periodSample(sRate,freq);
	double* data;
	data = malloc(tSample*sizeof(double));
	double* currentT = data;
	for(int i = 0; i < tSample;i++){
		*currentT = ampl * sin(2*3.1415*i/tSample);
		currentT++;
	}
	return data;
}
double* squareWave(int sRate, int freq,int ampl){
	double tSample = periodSample(sRate,freq);
	double* data;
	data = malloc(tSample*sizeof(double));
	double* currentT = data;
	for(int i = 0; i < tSample;i++){
		if(i < tSample/2){
			*currentT = 0;
		}else{
			*currentT = ampl;
		}
		currentT++;
	}
	return data;
}
double* trgWave(int sRate, int freq,int ampl){
	double tSample = periodSample(sRate,freq);
	double* data;
	data = malloc(tSample*sizeof(double));
	double* currentT = data;
	double delta = ampl/(tSample/2);
	for(int i = 0; i < tSample;i++){
		if(i < tSample/2){
			*currentT = delta*i;
		}else{
			*currentT = ampl-(i-(tSample/2))*delta;
		}
		currentT++;
	}
	return data;
}
double* sawWave(int sRate, int freq,int ampl){
	double tSample = periodSample(sRate,freq);
	double* data;
	data = malloc(tSample*sizeof(double));
	double* currentT = data;
	double delta = ampl/tSample;
	for(int i = 0; i < tSample;i++){
		*currentT = delta*i;
		currentT++;
	}
	return data;
}
double* linearCombine(double* waveA, double* waveB,double factor,int tSample){
	double* temp;
	temp = malloc(tSample*sizeof(double));
	double factorA = factor/100;
	printf("FactorA = %f \n",factorA);
	double factorB = 1 - factorA;
	printf("FactorB = %f \n",factorB);
	double* currentT = temp;
	double* currentA = waveA;
	double* currentB = waveB;
	for(int i = 0; i < tSample; i++){
		*currentT = (*currentA * factorA) + (*currentB * factorB);
		currentT++;
		currentA++;
		currentB++;
	}
	return temp;
}
int main(int argc, char **argv)
{
	FILE *f = fopen("test.csv", "w");
	if(f == NULL){
		printf("Error");
		return -1;
	}
	int FREQUNECY;
	printf("Ciao! \n");
	printf("Sample rate: %i sample/s\n",SAMPLERATE);
	printf("Inserisci frequezta:");
	scanf("%d",&FREQUNECY);
	printf("Frequency: %i Hz\n",FREQUNECY);
	printf("Amplitude: %i V\n",AMPLITUDE);
	int samples = periodSample(SAMPLERATE,FREQUNECY);
	printf("Target sample %i\n",samples);
	printf("Enter tipe of wave (Sq,Si,Tr,St)");
	char str[10];
	scanf("%9s",str);
	double* data;
	if(strcmp(str,"Sq") == 0){	
		printf("\n You have selected: Square Wave \n");
		data = squareWave(SAMPLERATE,FREQUNECY,AMPLITUDE);
	} else if(strcmp(str,"Si") == 0){
		printf("\n You have selected: Sine Wave \n");
		data = sineWave(SAMPLERATE,FREQUNECY,AMPLITUDE);
	}else if(strcmp(str,"Tr") == 0){
		printf("\n You have selected: Triange Wave \n");
		data = trgWave(SAMPLERATE,FREQUNECY,AMPLITUDE);
	}else if(strcmp(str,"St") == 0){
		printf("\n You have selected: Saw Tooth Wave \n");
		data = sawWave(SAMPLERATE,FREQUNECY,AMPLITUDE);
	}else{
		printf("\n Invalid input \n");
		return 0;
	}
	printf("Do you want to make a combined waveform? (S/N)");
	char strB[10];
	scanf("%9s",strB);
	if(strcmp(strB,"N")==0){
		double* cSample = data;
		for(double i = 0; i < samples; i++){
			double t = (i*1000)/SAMPLERATE;
			fprintf(f,"%f,%f \n",t,*cSample);
			cSample++;
		}
		fclose(f);
	}else{
		printf("Enter tipe of second wave (Sq,Si,Tr,St)");
		char strC[10];
		scanf("%9s",strC);
		double* dataB;
		if(strcmp(strC,"Sq") == 0){	
			printf("\n You have selected: Square Wave \n");
			dataB = squareWave(SAMPLERATE,FREQUNECY,AMPLITUDE);
		} else if(strcmp(strC,"Si") == 0){
			printf("\n You have selected: Sine Wave \n");
			dataB = sineWave(SAMPLERATE,FREQUNECY,AMPLITUDE);
		}else if(strcmp(strC,"Tr") == 0){
			printf("\n You have selected: Triange Wave \n");
			dataB = trgWave(SAMPLERATE,FREQUNECY,AMPLITUDE);
		}else if(strcmp(strC,"St") == 0){
			printf("\n You have selected: Saw Tooth Wave \n");
			dataB = sawWave(SAMPLERATE,FREQUNECY,AMPLITUDE);
		}else{
			printf("\n Invalid input \n");
			return 0;
		}
		printf("Enter modulation facotr (%%)");
		int factor;
		scanf("%d",&factor);
		double* finalWave = linearCombine(data,dataB,factor,samples);
		double* cSample = finalWave;
		for(double i = 0; i < samples; i++){
			double t = (i*1000)/SAMPLERATE;
			fprintf(f,"%f,%f \n",t,*cSample);
			cSample++;
		}
		fclose(f);
		
	}


	
	
	
	return 0;
	
}