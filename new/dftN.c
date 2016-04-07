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
#define OUTFILE "out.csv"
#define OUTBFILE "outB.csv"
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
void applyFilter(fftw_complex in[],fftw_complex out[],double filter[], int len){
	for(int i = 0; i < len; i++){
		out[i][0] = in[i][0] * filter[i];
		out[i][1] = in[i][1] * filter[i];
	}
}
int main(int argc, char **argv)
{
	//Greating
	printf("Frequency calculator\n");
	printf("Sample rate: %i\n",SAMPLERATE);
	printf("Input: %9s\n",INPUTFILE);
	//Data input
	double in[maxSample];
	double time[maxSample];
	int len = loadCSV(INPUTFILE,&in[0],&time[0]);
	//Data cleanup
	double cleanIn[len];
	for(int i = 0; i < len; i++){
		cleanIn[i] = in[i];
	}
	//Creating complex array
	fftw_complex timeDom[len];
	for(int i = 0; i < len; i++){
		timeDom[i][0] = cleanIn[i];
		timeDom[i][1] = 0;
	}
	//Computing FFT
	fftw_complex freqDom[len];
	fftw_plan p = fftw_plan_dft_1d(len,timeDom,freqDom,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	//Calculating freqs list
	double freqs[len];
	double factor = SAMPLERATE/len;
	for(int i = 0; i < len; i++){
		freqs[i] = i * factor;
	}
	//Computing reverse FFT
	fftw_complex timeDomB[len];
	fftw_plan q = fftw_plan_dft_1d(len,freqDom,timeDomB,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(q);
	fftw_destroy_plan(q);
	//Converting FFT data to real
	double out[len];
	double scaleFactor = 1./len;
	for(int i = 0; i < len; i++){
		out[i] = scaleFactor * timeDomB[i][0];
	}
	//Outputting reverse FFT data
	writeCSV(&out[0], &time[0],len, OUTFILE);

	//Get filter params
	double cutoff;
	double slope;
	printf("Insert low pass paramether: cutoff frequency (Hz):");
	scanf("%lf",&cutoff);
	printf("\n Insert low pass paramether: slope (G/KHz):");
	scanf("%lf",&slope);
	printf("\n");
	
	//Creating filter
	double lowPass[len];
	createLowPass(lowPass, len, cutoff, slope);
	
	//Appling filter
	fftw_complex freqDomB[len];
	applyFilter(freqDom,freqDomB,lowPass,len);
	
	//Complute filtered reverse FFT
	fftw_complex timeDomC[len];
	fftw_plan r = fftw_plan_dft_1d(len,freqDomB,timeDomC,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(r);
	fftw_destroy_plan(r);
	///Converting FFT data to real
	double outB[len];
	for(int i = 0; i < len; i++){
		outB[i] = scaleFactor * timeDomC[i][0];
	}
	//Outputting reverse FFT data
	writeCSV(&outB[0], &time[0],len, OUTBFILE);
	
	//Outputting FFT data
	writeCSVComplex(freqDom,&freqs[0],len,"FreqDom.csv");
	//Outputting filter data
	writeCSV(&lowPass[0],&freqs[0],len,"Filter.csv");
	//Outputting filtered FFT data
	writeCSVComplex(freqDomB,&freqs[0],len,"FreqDomB.csv");
	
	//END OF PROGRAM
	fftw_cleanup();
	return 0;

}