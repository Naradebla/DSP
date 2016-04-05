#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define N 16

int main(void){
	fftw_complex in[N], out[N], in2[N];
	fftw_plan p,q;
	for(int i = 0; i < N; i++){
		in[i][0] = cos(3*2*M_PI*i/N);
		in[i][1] = 0;
	}
	p = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	for(int i = 0; i < N; i++){
		printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][0],out[i][1]);
	}
	
}