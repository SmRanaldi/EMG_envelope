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
#define MAX_ITER 20

void filterLength(double *_m, double *_env, double *_d, double *_d2, double *_idx, int _nSamples);
void envelopeEstimation(double *_m, double *_signal, double *_env, double *_idx, int _nSamples);
void derivativesEstimation(double *_m, double *_signal, double *_d1, double *_d2, double *_idx, int _nSamples);
void updateM(double *_tmpM, double *_outM, double *_idx, int _nSamples);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	// ---- Inputs ----
    int nSamples; 

    double *m, *chiSq_in, *chiSq, *idx;
    double *signalMAT, *envelope, *d1, *d2;
    double *signal;

	// ---- Outputs ----
    double *mOut, *envOut, *d1Out, *d2Out;
    double *mOut_mx, *envOut_mx;
	
	// ---- Temporary variables ----
    int i,j,k,MChi,NChi,nConv;
    double tmpChi, tmpEnt;
    double *entT2, *entT1, *mTmp;
	
	// ---- Input initialization ----
	nSamples=(int) mxGetN(prhs[0]);
    signalMAT = mxGetPr(prhs[0]);
    envelope = mxGetPr(prhs[1]);
    d1 = mxGetPr(prhs[2]);
    d2 = mxGetPr(prhs[3]);
    m = mxGetPr(prhs[4]);
    chiSq_in = mxGetData(prhs[5]);
    MChi = mxGetM(prhs[5]);
    NChi = mxGetN(prhs[5]);

    idx = mxCalloc(nSamples,sizeof(double));
    chiSq = mxCalloc(MChi*NChi,sizeof(double));

    for(i=0;i<MChi;i++){
        for(j=0;j<NChi;j++){
            chiSq[i+j*MChi] = chiSq_in[i+j*MChi];
        }
    }

	// ---- Output initialization ----
	plhs[0] = mxCreateDoubleMatrix(1,nSamples,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1,nSamples,mxREAL);
    mOut_mx = mxGetPr(plhs[0]);
    envOut_mx = mxGetPr(plhs[1]);

    // ---- Tmp initialization ----
    entT1 = mxCalloc(nSamples,sizeof(double));
    entT2 = mxCalloc(nSamples,sizeof(double));
    
    // ---- Copy arrays ----
    mOut = mxCalloc(nSamples,sizeof(double));
    mTmp = mxCalloc(nSamples,sizeof(double));
    envOut = mxCalloc(nSamples,sizeof(double));
    d1Out = mxCalloc(nSamples,sizeof(double));
    d2Out = mxCalloc(nSamples,sizeof(double));
    signal = mxCalloc(nSamples,sizeof(double));
    for(i=0;i<nSamples;i++){
        mOut[i] = m[i];
        envOut[i] = envelope[i];
        d1Out[i] = d1[i];
        d2Out[i] = d2[i];
        signal[i] = fabs(signalMAT[i]);
    }
    
    // ---- Update step ----
    nConv = 0;
    i=0;
    while(i<MAX_ITER){
        for(j=0;j<nSamples;j++){
            mTmp[j] = mOut[j];
        }
        filterLength(mTmp, envOut, d1Out, d2Out, idx, nSamples);
        envelopeEstimation(mTmp, signal, envOut, idx, nSamples);
        derivativesEstimation(mTmp, signal, d1Out, d2Out, idx, nSamples);
        for(j=0;j<nSamples;j++){
            if(idx[j]==0){
                tmpEnt = 0.0;
                for(k=0;k<NChi;k++){
                    if(chiSq[(int)mOut[j]-1+k*MChi]!=0.0){
                        tmpEnt += chiSq[(int)mOut[j]-1+k*MChi]*log(chiSq[(int)mOut[j]-1+k*MChi]);
                    }
                }
                if(i>1){
                    // mexPrintf("!\n");
                    if(((tmpEnt-entT1[j]) - (entT1[j]-entT2[j]))>0.0){
                        idx[j] = 1;
                        nConv++;
                    }
                }
            entT2[j] = entT1[j];
            entT1[j] = tmpEnt;
            }
        }
        updateM(mTmp, mOut, idx, nSamples);
        if((double)nConv/nSamples > 0.95){
            mexPrintf("Convergence at iteration: %d\n",i+1);
            i = MAX_ITER;
        }
        i++;
    }

    // ---- Copy output ----
    for(i=0;i<nSamples;i++){
        mOut_mx[i] = mOut[i];
        envOut_mx[i] = envOut[i];
    }

    return;

}

void filterLength(double *_m, double *_env, double *_d, double *_d2, double *_idx, int _nSamples){

    int updateSample, mTmp;
    double A, B;
    double numerator, den;
    double d1Tmp, d2Tmp, envTmp;

    for(updateSample=0;updateSample<_nSamples;updateSample++){
        if(_idx[updateSample] == 0){
            mTmp = _m[updateSample];
            d1Tmp = _d[updateSample];
            d2Tmp = _d2[updateSample];
            envTmp = _env[updateSample];
            A = -0.5*d1Tmp;
            B = (1.0/6) * (d2Tmp + pow(d1Tmp,2.0)/(4*envTmp));

            numerator = 4 * F_FACTOR * pow(envTmp,4.0);
            den = pow((B*envTmp + pow(A,2.0))/2,2.0);

            mTmp = round(pow(fabs(numerator/den),1/5.0));
            if(mTmp<1){ mTmp=1; }
            if(mTmp>10000){ mTmp=10000; }

            _m[updateSample] = mTmp;
        }

    }

    return;
}


void envelopeEstimation(double *_m, double *_signal, double *_env, double *_idx, int _nSamples){

    int i, j, semiLen, upperLimit, lowerLimit, winLen, updateSample;
    double wSample;

    for(updateSample=0;updateSample<_nSamples;updateSample++){
        if(_idx[updateSample] == 0){
            semiLen = ceil(_m[updateSample] * 0.5);
            upperLimit = MIN(_nSamples,updateSample+semiLen);
            lowerLimit = MAX(0,updateSample-semiLen);
            winLen = upperLimit - lowerLimit + 1;
            wSample = 0.0;
            for(j=lowerLimit;j<upperLimit+1;j++){
                wSample += pow((_signal[j]),NU);
            }
            wSample /= winLen;
            _env[updateSample] = pow(wSample/P_NORM, 0.5);
        }
    }
    return;
}

void derivativesEstimation(double *_m, double *_signal, double *_d1, double *_d2, double *_idx, int _nSamples){

    int i, j, a, semiLen, upperLimit, lowerLimit, winLen, updateSample;
    double est1, est2_1, est2_2, d1Factor, d2Factor;

    for(updateSample=0;updateSample<_nSamples;updateSample++){
        if(_idx[updateSample] == 0){
            semiLen = ceil(_m[updateSample] * 0.5);
            upperLimit = MIN(_nSamples,updateSample+semiLen);
            lowerLimit = MAX(0,updateSample-semiLen);
            winLen = upperLimit - lowerLimit + 1;
            d1Factor = 0.0;
            d2Factor = 0.0;
            est1 = 0.0;
            est2_1 = 0.0;
            est2_2 = 0.0;
            for(j=lowerLimit;j<upperLimit+1;j++){
                a = j - lowerLimit - ceil(0.5*winLen); // Check
                d1Factor += pow(a,2.0);
                d2Factor += pow(a,4.0);
            }
            for(j=lowerLimit;j<upperLimit+1;j++){
                a = j - lowerLimit - ceil(0.5*winLen); // Check
                est1 += a * (_signal[j]);
                est2_1 += pow(a,2.0)*(_signal[j]);
                est2_2 += (1 - pow(a,2.0)*(d1Factor/d2Factor))*(_signal[j]);
            }
            _d2[updateSample] = 2 * (est2_1/(d2Factor*P_NORM) - (d1Factor/(d2Factor*P_NORM)) * est2_2 /(winLen + (d1Factor/d2Factor)*d1Factor));
            _d1[updateSample] = est1/(d1Factor*P_NORM);
        }
    }
    return;
}

void updateM(double *_tmpM, double *_outM, double *_idx, int _nSamples){
    int i;

    for(i=0;i<_nSamples;i++){
        if(_idx[i]==0.0){
            _outM[i] = _tmpM[i];
        }
    }
    return;
}