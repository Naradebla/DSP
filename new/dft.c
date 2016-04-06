#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include "utils.h"

#define FillArray(a,b,c) memset(a,c,b)

#ifndef maxSample
#define maxSample 1000
#endif
#define pi 3.1415
#define SAMPLERATE 100000

#define INPUTFILE "test.csv"
#define FFTFILE "complex.csv"
#define FILTERFILE "lp.csv"
#define DPFILE "dp.csv"
void createLowPass(double target[], int len, double cutOffFreq, double slopeFreq){
	double factor = SAMPLERATE/len;
	double cutOff = cutOffFreq/factor;
	double slope = (1000*slopeFreq/factor);
	double startF = ((slope*cutOff)-0.5)/slope;
	double endF = (1/slope)+startF;
	printf("Factor = %f\n",factor);
	printf("Normalized cutoff = %f\n",cutOff);
	printf("Normalized slope = %f\n",slope);
	printf("Start freq = %f\n",startF);
	printf("End freq = %f\n",endF);
	for(double i = 0; i < len; i++){
		int h = i;
		if(i < startF){
			target[h] = 1;
			continue;
		}
		if(i > endF){
			target[h] = 0;
			continue;
		} 
		target[h] = 1 - ((i - startF)*slope);
	}
}
void getMaxId(fftw_complex data[],int out[],int len){
	int maxRe = data[0][0];
	int maxIm = data[0][1];
	int reId = 0; 
	int imId = 0;
	for(int i = 0; i < len; i++){
		if(data[i][0] > maxRe){
			maxRe = data[i][0];
			reId = i;
		}
		if(data[i][1] > maxIm){
			maxIm = data[i][1];
			imId = i;
		}
	}
	out[0] = reId;
	out[1] = imId;
}

int shazam(double maxF, double samples[],fftw_complex results[],int len){
	int sampleFactor = SAMPLERATE/maxF;
	int targetLen = len * sampleFactor;
	int t = 0;
	double factor = SAMPLERATE/targetLen;
	fftw_complex timeD[targetLen],freqD[targetLen];
	for(int i = 0; i < sampleFactor; i++){
		for(int c = 0; c < len; c++){
			timeD[t][0] = samples[c];
			timeD[t][1] = 0;
			t++;
		}
	}
	fftw_plan p;
	p = fftw_plan_dft_1d(targetLen,timeD,freqD,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	for(int i = 0; i < len;i++){
		results[i][0] = freqD[i][0];
		results[i][1] = freqD[i][1];
	}
	fftw_destroy_plan(p);
	return factor;
}

int main(int argc, char **argv)
{
	printf("Frequency calculator\n");
	printf("Sample rate: %i\n",SAMPLERATE);
	printf("Input: %9s\n",INPUTFILE);
	
	double in[maxSample];
	double timeA[maxSample];
	int len = loadCSV(INPUTFILE,&in[0],&timeA[0]);
	
	double cleanIn[len];
	for(int i = 0; i < len; i++){
		cleanIn[i] = in[i];
	}
	fftw_complex freqDom[len];
	double factor = shazam(20000,cleanIn,freqDom,len);
	/*
	fftw_complex timeDom[len], freqDom[len];
	fftw_plan p;
	
	*/
	
	/*double f = timeA[1]-timeA[0];
	double time[len];
	for(int i = 0; i < len; i++){
		time[i] = f*i;
		if(i < lenA){
			timeDom[i][0] = cleanIn[i];
			timeDom[i][1] = 0;
		}else{
			timeDom[i][0] = cleanIn[i - lenA];
			timeDom[i][1] = 0;
		}
	}*/
	/*writeCSVComplex(timeDom, time,len, DPFILE);
	
	
	p = fftw_plan_dft_1d(len,timeDom,freqDom,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);*/

	//double factor = SAMPLERATE/len;
	
	double freqs[len];
	for(int i = 0; i < len; i++){
		freqs[i] = i * factor;
	}
	
	writeCSVComplex(freqDom, &freqs[0],len, FFTFILE);
	
	double magnitude[len];
	for(int i = 0; i < len; i++){
		magnitude[i] = sqrt(pow(freqDom[i][0],2)+pow(freqDom[i][1],2));
	}
	
	double maxPeak = 0;
	int index = 0;
	for(int i = 0;i < len; i++){
		if(magnitude[i] > maxPeak){
			maxPeak = magnitude[i];
			index = i;
		}
	}
	
	printf("Magniute peak detection, index: %i, scaled: %f\n",index,index*factor);
	
	double cof;
	double slope;
	printf("Generating low pass filter...\nInsert cut off frequency (Hz):");
	scanf("%lf",&cof);
	printf("\nInsert slope: (DG/DKHz)");
	scanf("%lf",&slope);
	printf("\n");
	
	printf("Samples: %i, max freq:%f\n",len,len*factor);
	
	double lowPass[len];
	createLowPass(lowPass,len,cof,slope);
	
	writeCSV(&lowPass[0], &freqs[0],len, FILTERFILE);
	fftw_cleanup();
	return 0;
}