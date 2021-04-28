# EMG_envelope
Algorithm for the automatic estimation of the RMS envelope of the surface electromyography signal.

## Description

This package contains the MATLAB codes for implementing the algorithm described in:

  * S.Ranaldi, C. De Marchis and S. Conforto "An automatic, adaptive, information-based method for the extraction of the sEMG envelope"
https://doi.org/10.1016/j.jelekin.2018.06.001

This package contains some mex function to speed up the algorithm using C functions.

## List of functions

MATLAB codes:

* adaptiveEnvelope.m - Main function
* conditionEMG.m - Conditioning block (whitening and normalization)
* derivativesEstimation.m - Point by point first and second derivative of the envelope estimation
* envelopeEstimation.m - Point by point envelope computation
* entropyEst.m - Point by point entropy estimation for convergence
* f.m - Normalization factor
* filterLength.m - Update of the adaptive filter window lengths
* staticEstimationD.m - Initialization of the estimation of the derivatives (might be removed in the future)
* staticEstimationW.m - Initialization of the estimation of the envelope (might be removed in the future)
* whiteTest.m - Tests for the whiteness of the signal (borrowed function, source in the comments)
* whitenSignal.m - Whitening filter in MATLAB

C codes:

* posAutoCorr.c - Correlation
* whiteningSignal - Whitening filter in C

## Install

Make sure that the MATLAB C package is installed. To use the adaptiveEnvelope.m function add its path and its subfolder to MATLAB.

C functions have been compiled on Windows 10. MacOs or Linux users must recompile all the C codes in order for them to work on their systems.
