#include <stdio.h>
#include "mex.h"

#define MIN(x,y) (x<y)?x:y
#define MAX(x,y) (x>y)?x:y

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){	
	
	double *signal, *out;
	int l;
	
	double tmpCorr;
	int idxS, idxE;
	int i, j, k;
	
	signal = mxGetPr(prhs[0]);
	l = *mxGetPr(prhs[1]);
	
	plhs[0] = mxCreateDoubleMatrix(l,1,mxREAL);
	out = mxGetPr(plhs[0]);
	
	for(i=0;i<l;i++){
		
		tmpCorr=0;
		
		for(j=i;j<l;j++){
			
			tmpCorr += signal[j]*signal[j-i];
			
		}
		
		out[i]=tmpCorr/(l-i);
		
	}
	
}