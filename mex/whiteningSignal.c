#include <stdio.h>
#include "mex.h"
#include <math.h>

#define MAXORDER 30

#define MIN(x,y) (x<y)?x:y
#define MAX(x,y) (x>y)?x:y

/*-------------------------------------------------------------------/
/																	 /
/     signalOut = whitening(signalIn, predictorOrder, windowLength)	 /
/																	 /
/-------------------------------------------------------------------*/

/* --------- Function declarations --------- */

void whitening(double *signal, int tmp1, int tmp2, int order, double *signalOut, double *coefficients);
void autoCorrelation(double *signal1, int tmp1, int tmp2, int maxLag, long double *out);
void arCoefficients(double *signal, int order, long double *autocorr, double *coefficients);

/* --------- Main --------- */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	// ---- Inputs ----
	double *signal;
	double *predictorOrder, *windowLength;
	
	// ---- Outputs ----
	double *whiteSignal, *coefficients;
	
	// ---- Temporary variables ----
	int i, semiLen, nSamples;
	int tmp1, tmp2, tmpOrder;
	
	
	// ---- Input initialization ----
	signal = mxGetPr(prhs[0]);
	predictorOrder = mxGetPr(prhs[1]);
	windowLength = mxGetPr(prhs[2]);
	
	semiLen = (int) floor(*windowLength/2);
	nSamples=mxGetM(prhs[0]);
	
	// ---- Output initialization ----
	plhs[0] = mxCreateDoubleMatrix(nSamples,1,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(*predictorOrder,1,mxREAL);
	whiteSignal = mxGetPr(plhs[0]);
	coefficients = mxGetPr(plhs[1]);
	
	if(*predictorOrder>MAXORDER) *predictorOrder=MAXORDER;
	
		
	// ---- Filtering ----
	for(i=0;i<nSamples;i++){
		
		tmp1=MAX(0,i-semiLen);
		tmp2=MIN(nSamples-1,i+semiLen);
		
		tmpOrder = MIN(*predictorOrder,(int)(tmp2-tmp1+1)/2);
				
		whitening(signal,tmp1,tmp2,tmpOrder,&whiteSignal[i],coefficients);
		
	}

}

/* --------- Functions definition --------- */

void whitening(double *signal, int tmp1, int tmp2, int order, double *sampleOut, double *coefficients){		// AR filtering for whitening.
	
	long double autoCorr[order+1];
	int i, x;

	x=(int) (tmp2+tmp1)/2;
	
	autoCorrelation(signal, tmp1, tmp2, order+1, autoCorr);
	
	arCoefficients(signal, order, autoCorr, coefficients);
	
	*sampleOut = signal[x];
	
	for(i=1;i<=order;i++){
		
		*sampleOut += coefficients[i-1]*signal[x-i];
		
	}
			
	return;
	
}

void autoCorrelation(double *signal, int tmp1, int tmp2, int maxLag, long double *out){				// Autocorrelation function
	
	double tmpCorr;
	int i, j;
	
	for(i=0;i<=maxLag;i++){
		
		tmpCorr=0;
		
		for(j=tmp1+i;j<tmp2;j++){
			
			tmpCorr += signal[j]*signal[j-i];
			
		}
		
		out[i]=tmpCorr/(tmp2-tmp1);
		
	}
	
	return;
	
}

void arCoefficients(double *signal, int order, long double *autocorr, double *coefficients){
	
	long double a[MAXORDER][MAXORDER]={{0}};
	long double sigma[order];
	int i,j;
		
	a[0][0] = - autocorr[1]/autocorr[0];
	sigma[0] = (1 - pow(a[0][0],2.))*autocorr[0];
	
	for(i=1;i<order;i++){
		
		a[i][i] = autocorr[i+1];
		
		for(j=0;j<i;j++){
			
			a[i][i] += a[i-1][j]*autocorr[i-j];
			
		}
		
		a[i][i] /= -sigma[i-1];
		
		for(j=0;j<i;j++){
			
			a[i][j] = a[i-1][j]+a[i][i]*a[i-1][i-j-1];
			
		}
		
		sigma[i] = (1 - pow(a[i][i],2.))*sigma[i-1];
		
	}
	
	for(i=0;i<order;i++){
		
		coefficients[i] = a[order-1][i];
		
	}	
	
}
