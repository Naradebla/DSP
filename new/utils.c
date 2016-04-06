#include "utils.h"

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
void writeCSV(double* x, double* y,int len, char* path){
	FILE *f = fopen(path, "w");
	for(double i = 0; i < len; i++){
		fprintf(f,"%f,%f \n",*y,*x);
		x++;
		y++;
	}
	fclose(f);
}
void writeCSVComplex(fftw_complex x[], double* y,int len, char* path){
	FILE *f = fopen(path, "w");
	for(int i = 0; i < len; i++){
		fprintf(f,"%f,%f,%f \n",*y,x[i][0],x[i][1]);
		y++;
	}
	fclose(f);
}