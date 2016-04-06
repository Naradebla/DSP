#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#define ZeroMemory(a,b) memset(a,0,b)
#define maxSample 1000
int loadCSV(char* path,double* outA,double* outB);
void writeCSV(double* x, double* y,int len, char* path);
void writeCSVComplex(fftw_complex x[], double* y,int len, char* path);

#endif