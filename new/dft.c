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
#define OUTFREQFILE "ffout.csv"
#define OUTTIMEFILE "tout.csv"
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
void applyFilter(fftw_complex signal[],fftw_complex out[],double filter[],int len){
	for(int i = 0; i < len; i++){
		out[i][0] = signal[i][0] * filter[i];
		out[i][1] = signal[i][1] * filter[i];
	}
	
}
int hdfft(int sampleFactor, double samples[],fftw_complex results[],int len){
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
	/*for(int i = 0; i < len;i++){
		results[i][0] = freqD[i][0];
		results[i][1] = freqD[i][1];
	}*/
	fftw_destroy_plan(p);
	return targetLen;
}
int main(int argc, char **argv)
{
	printf("Frequency calculator\n");
	printf("Sample rate: %i\n",SAMPLERATE);
	printf("Input: %9s\n",INPUTFILE);
	
	double in[maxSample];
	double time[maxSample];
	int len = loadCSV(INPUTFILE,&in[0],&time[0]);
	
	double cleanIn[len];
	for(int i = 0; i < len; i++){
		cleanIn[i] = in[i];
	}
	fftw_complex freqDom[len*5];
	int eLen = hdfft(5,cleanIn,freqDom,len);

	double factor = SAMPLERATE/eLen;
	double freqs[eLen];
	for(int i = 0; i < eLen; i++){
		freqs[i] = i * factor;
	}
	
	writeCSVComplex(freqDom, &freqs[0],eLen, FFTFILE);
	
	double magnitude[eLen];
	for(int i = 0; i < eLen; i++){
		magnitude[i] = sqrt(pow(freqDom[i][0],2)+pow(freqDom[i][1],2));
	}
	
	double maxPeak = 0;
	int index = 0;
	for(int i = 0;i < eLen; i++){
		if(magnitude[i] > maxPeak){
			maxPeak = magnitude[i];
			index = i;
		}
	}
	
	printf("Magniute peak detection, index: %i, scaled: %f\n",index,index*factor);
	
	/*double cof;
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
	
	fftw_complex signal[eLen];
	applyFilter(freqDom, signal, lowPass, len);
	writeCSVComplex(signal,&freqs[0],len, OUTFREQFILE);*/
	
	fftw_complex outTimeDom[eLen];
	
	fftw_plan q;
	q = fftw_plan_dft_1d(eLen,outTimeDom,freqDom,FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(q);
	
	writeCSVComplex(outTimeDom,&time[0],eLen, OUTTIMEFILE);
	
	fftw_destroy_plan(q);
	fftw_cleanup();
	return 0;
}