#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#define ZeroMemory(a,b) memset(a,0,b)
#define FillArray(a,b,c) memset(a,c,b)

#define maxSample 10000
#define pi 3.1415
#define SAMPLERATE 20000
int loadCSV(char* path,double* outA,double* outB){
	FILE *fp;
	fp = fopen(path,"r");
	char buffer;
	char dataBuffer[10];
	char dataBufferB[10];
	double outputBuffer[maxSample];
	double outputBufferB[maxSample];
	
	int i = 0;
	int s = 0;
	int p = 0;
	int c = 0;
	while(1){
		buffer = fgetc(fp);
		if(feof(fp)){
			break;
		}
		if(buffer == ','){
			s = 1;	
		}else if(buffer == '\n'){
			s = 0;
			i = 0;
			c = 0;
			outputBuffer[p] = atof(dataBuffer);
			outputBufferB[p] = atof(dataBufferB);
			ZeroMemory(dataBuffer,10);
			ZeroMemory(dataBufferB,10);
			p++;	
		}else{
			if(s == 1){
				dataBuffer[i] = buffer;
				i++;
			}else{
				dataBufferB[c] = buffer;
				c++;
			}
		}
	}
	for(int i = 0; i < p; i++){
		*outA = outputBuffer[i];
		*outB = outputBufferB[i];
		outA++;
		outB++;
	}
	return p;
	
}
void newDFT(const double inReal[],const double inImag[],double outReal[],double outImag[],int n){
	for(int k = 0; k < n; k++){
		double sumReal = 0;
		double sumImag = 0;
		for(int t = 0; t < n; t++){
			double angle = 2 * M_PI * t *k/n;
			sumReal += inReal[t] * cos(angle) + inImag[t] * sin(angle);
			sumImag += -inReal[t] * sin(angle) + inImag[t] * cos(angle);
		}
		outReal[k] = sumReal;
		outImag[k] = sumImag;
	}
}

void writeCSV(double* x, double* y,int len, char* path){
	FILE *f = fopen(path, "w");
	for(double i = 0; i < len; i++){
		fprintf(f,"%f,%f \n",*y,*x);
		x++;
		y++;
	}
	fclose(f);
}
void createLowPass(double target[], int len, double cutOff, double slope){
	double factor = SAMPLERATE/len;
	double startF = ((slope*(cutOff/factor))-0.5)/slope;
	double endF = (1/slope)+startF;
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
void applyFilter(double data[],double filter[],double out[],int len){
	for(int i = 0; i < len; i++){
		out[i] = data[i]*filter[i];
	}
}
void applyFilterComplex(fftw_complex data[],double filter[], fftw_complex out[],int len){
	for(int i = 0; i < len; i++){
		out[i][0] = data[i][0]*filter[i];
		out[i][1] = data[i][1]*filter[i];
	}
}
void calculateAmplitude(fftw_complex in[],double out[],int len){
	for(int i = 0; i < len; i++){
		out[i] = sqrt(pow(in[i][0],2)+pow(in[i][1],2));
	}
}
int main(int argc, char **argv)
{
	printf("Ciao!\n");
	printf("Sample rate: %i \n",SAMPLERATE);
	printf("Filtro LowPass, inserisci frequenza di taglio (Hz)");
	double cutoff;
	scanf("%lf",&cutoff);
	printf("\n Inserisci slope (Db/Dec)");
	double slope;
	scanf("%lf",&slope);
	printf("\n");
	double data[maxSample];
	double time[maxSample];
	int len = loadCSV("test.csv",&data[0],&time[0]);
	fftw_complex in[len], out[len];
	fftw_plan p,q;
	for(int i = 0; i < len; i++){
		in[i][0] = data[i];
		in[i][1] = 0;
	}
	p = fftw_plan_dft_1d(len,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	double meanOut[len];
	double freq[len];
	double factor = SAMPLERATE/len;
	for(int i = 0; i < len; i++){
		meanOut[i] = sqrt(pow(out[i][0],2)+pow(out[i][1],2));
		freq[i] = i * factor;
	}
	writeCSV(&meanOut[0],&freq[0],len,"DFT.csv");
	double lowPass[len];
	createLowPass(lowPass,len,cutoff,slope);
	writeCSV(&lowPass[0],&freq[0],len,"lowPass.csv");
	fftw_complex result[len];
	//applyFilter(meanOut,lowPass,result,len);
	applyFilterComplex(out,lowPass,result,len);
	double meanResult[len];
	calculateAmplitude(result,meanResult,len);
	writeCSV(&meanResult[0],&freq[0],len,"filteredFreqs.csv");
	fftw_complex filteredOut[len];
	q = fftw_plan_dft_1d(len,result,filteredOut,FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(q);
	for(int i = 0; i < len; i++){
		filteredOut[i][0] *= 1./len;
		filteredOut[i][1] *= 1./len;
	}
	for(int i = 0; i<len;i++){
		printf("recover: %3d %+9.5f %+9.5f I \n",i,filteredOut[i][0],filteredOut[i][1]);
	}
	calculateAmplitude(filteredOut,meanResult,len);
	writeCSV(&meanResult[0],&time[0],len,"out.csv");
	
	fftw_destroy_plan(q);
	fftw_destroy_plan(p);
	fftw_cleanup();
	return 0;
	
}