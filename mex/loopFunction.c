#include <stdio.h>
#include "mex.h"
#include <math.h>

#define ALPHA 1.0
#define NU 2.0
#define PI 3.141592653589793
#define P_NORM 0.797884560802866
#define F_FACTOR 0.535398163397448
#define MIN(x,y) (x<y)?x:y
#define MAX(x,y) (x>y)?x:y

void filterLength(double *_m, double *_env, double *_d, double *_d2, double *_idx, int _nSamples);
void envelopeEstimation(double *_m, double *_signal, double *_env, double *_idx, int _nSamples, int _lSignal);
void derivativesEstimation(double *_m, double *_signal, double *_d1, double *_d2, double *_idx, int _nSamples, int _lSignal);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	// ---- Inputs ----
    int nSamples, lSignal;

    double *m, *idx;
    double *signal, *envelope, *d1, *d2;

	// ---- Outputs ----
    double *mOut, *envOut, *d1Out, *d2Out;
	
	// ---- Temporary variables ----
    int i;
	
	// ---- Input initialization ----
	nSamples=(int) mxGetN(prhs[5]);
    lSignal=(int) mxGetN(prhs[0]);
    signal = mxGetPr(prhs[0]);
    envelope = mxGetPr(prhs[1]);
    d1 = mxGetPr(prhs[2]);
    d2 = mxGetPr(prhs[3]);
    m = mxGetPr(prhs[4]);
    idx = mxGetPr(prhs[5]);

	// ---- Output initialization ----
	plhs[0] = mxCreateDoubleMatrix(1,lSignal,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1,lSignal,mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1,lSignal,mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1,lSignal,mxREAL);
    mOut = mxGetPr(plhs[0]);
    envOut = mxGetPr(plhs[1]);
    d1Out = mxGetPr(plhs[2]);
    d2Out = mxGetPr(plhs[3]);
    
    // ---- Copy arrays ----
    for(i=0;i<lSignal;i++){
        mOut[i] = m[i];
        envOut[i] = envelope[i];
        d1Out[i] = d1[i];
        d2Out[i] = d2[i];
    }
    
    // ---- Update step ----
    filterLength(mOut, envOut, d1Out, d2Out, idx, nSamples);
    envelopeEstimation(mOut, signal, envOut, idx, nSamples, lSignal);
    derivativesEstimation(mOut, signal, d1Out, d2Out, idx, nSamples, lSignal);

}

void filterLength(double *_m, double *_env, double *_d, double *_d2, double *_idx, int _nSamples){

    int i, updateSample, mTmp;
    double A, B;
    double numerator, den;
    double d1Tmp, d2Tmp, envTmp;

    for(i=0;i<_nSamples;i++){
        updateSample = (int)_idx[i]-1;
        mTmp = _m[updateSample];
        d1Tmp = _d[updateSample];
        d2Tmp = _d2[updateSample];
        envTmp = _env[updateSample];
        A = -0.5*d1Tmp;
        B = (1.0/6) * (d2Tmp + ((ALPHA*NU - 1)*pow(d1Tmp,2.0))/(4*envTmp));

        numerator = 4 * F_FACTOR * pow(envTmp,4.0);
        den = pow((B*envTmp + (ALPHA*NU - 1)*pow(A,2.0))/2,2.0);

        mTmp = (int)round(pow(fabs(numerator/den),1/5.0));
        if(mTmp<1){ mTmp=1; }
        if(mTmp>10000){ mTmp=10000; }

        _m[updateSample] = mTmp;

    }

    return;
}


void envelopeEstimation(double *_m, double *_signal, double *_env, double *_idx, int _nSamples, int _lSignal){

    int i, j, semiLen, upperLimit, lowerLimit, winLen, updateSample;
    double wSample;

    for(i=0;i<_nSamples;i++){
        updateSample = (int) round(_idx[i])-1;
        semiLen = (int) ceil(_m[updateSample] * 0.5);
        upperLimit = MIN(_lSignal,updateSample+semiLen);
        lowerLimit = MAX(0,updateSample-semiLen);
        winLen = upperLimit - lowerLimit + 1;
        wSample = 0.0;
        for(j=lowerLimit;j<upperLimit+1;j++){
            wSample += pow(fabs(_signal[j]),NU);
        }
        wSample /= winLen;
        _env[updateSample] = pow(wSample/P_NORM, 1/(ALPHA*NU));
    }
}

void derivativesEstimation(double *_m, double *_signal, double *_d1, double *_d2, double *_idx, int _nSamples, int _lSignal){

    int i, j, a, semiLen, upperLimit, lowerLimit, winLen, updateSample;
    double est1, est2_1, est2_2, d1Factor, d2Factor;

    for(i=0;i<_nSamples;i++){
        updateSample = (int) round(_idx[i])-1;
        semiLen = (int) ceil(_m[updateSample] * 0.5);
        upperLimit = MIN(_lSignal,updateSample+semiLen);
        lowerLimit = MAX(0,updateSample-semiLen);
        winLen = upperLimit - lowerLimit + 1;
        d1Factor = 0.0;
        d2Factor = 0.0;
        est1 = 0.0;
        est2_1 = 0.0;
        est2_2 = 0.0;
        for(j=lowerLimit;j<upperLimit+1;j++){
            a = j - lowerLimit - (int) ceil(0.5*winLen); // Check
            d1Factor += pow(a,2.0);
            d2Factor += pow(a,4.0);
        }
        for(j=lowerLimit;j<upperLimit+1;j++){
            a = j - lowerLimit - (int) ceil(0.5*winLen); // Check
            est1 += a * pow(fabs(_signal[j]),1/ALPHA);
            est2_1 += pow(a,2.0)*pow(fabs(_signal[j]),1/ALPHA);
            est2_2 += (1 - pow(a,2.0)*(d1Factor/d2Factor))*pow(fabs(_signal[j]),1/ALPHA);
        }
        _d2[updateSample] = 2 * (est2_1/(d2Factor*P_NORM) - (d1Factor/(d2Factor*P_NORM)) * est2_2 /(winLen + (d1Factor/d2Factor)*d1Factor));
        _d1[updateSample] = est1/(d1Factor*P_NORM);
    }
}