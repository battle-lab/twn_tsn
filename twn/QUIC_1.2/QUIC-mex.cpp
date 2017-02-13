// This is the MEX wrapper for QUIC.  The algorithm is in QUIC.C.

// Invocation form within Matlab or Octave:
// [X W opt time iter dGap] = QUIC(mode, ...)
// [X W opt time iter dGap] = QUIC("default", S, L, tol, msg, maxIter,
//                            X0, W0)
// [X W opt time iter dGap] = QUIC("path", S, L, path, tol, msg, maxIter,
//                            X0, W0)
// [X W opt time iter dGap] = QUIC("trace", S, L, tol, msg, maxIter,
//                                 X0, W0)
// See the README file and QUIC.m for more information.

#include <mex.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

extern "C"
void QUIC(char mode, uint32_t& p, const double* S, double* Lambda,
	  uint32_t& pathLen, const double* path, double& tol,
	  int32_t& msg, uint32_t& maxIter,
	  double* X, double* W, double* opt, double* time,
	  uint32_t* iter, double* dGap);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2) {
	mexErrMsgIdAndTxt("QUIC:arguments",
			  "Missing arguments, please specify\n"
	    "             S - the empirical covariance matrix, and\n"
	    "             L - the regularization parameter.");
    }
    long argIdx = 0;
    char mode[8];
    mode[0] = 'd';
    if (mxIsChar(prhs[0])) {
	mxGetString(prhs[0], mode, 8);
	if (strcmp(mode, "path") &&
	    strcmp(mode, "trace") &&
	    strcmp(mode, "default"))
	    mexErrMsgIdAndTxt("QUIC:arguments",
			      "Invalid mode, use: 'default', 'path' or "
			      "'trace'.");
	argIdx++;
    }
    // The empirical covariance matrix:
    if (!mxIsDouble(prhs[argIdx]))
	mexErrMsgIdAndTxt("QUIC:type",
			  "Expected a double matrix. (Arg. %d)",
			  argIdx + 1);
    const double* S = mxGetPr(prhs[argIdx]);
    uint32_t p = (uint32_t) mxGetN(prhs[argIdx]);
    if (p != uint32_t(mxGetM(prhs[argIdx]))) {
	mexErrMsgIdAndTxt("QUIC:dimensions",
			  "Expected a square empirical covariance matrix.");
    }
    argIdx++;

    // Regularization parameter matrix:
    if (!mxIsDouble(prhs[argIdx]))
	mexErrMsgIdAndTxt("QUIC:type",
			  "Expected a double matrix. (Arg. %d)",
			  argIdx + 1);
    double* Lambda;
    unsigned long LambdaAlloc = 0;
    if (mxGetN(prhs[argIdx]) == 1 && mxGetM(prhs[argIdx]) == 1) {
	Lambda = (double*) malloc(p*p*sizeof(double));
	LambdaAlloc = 1;
	double lambda = mxGetPr(prhs[argIdx])[0];
	for (unsigned long i = 0; i < p*p; i++)
	    Lambda[i] = lambda;
    } else {
	if (mxGetN(prhs[argIdx]) != p && mxGetM(prhs[argIdx]) != p) {
	    mexErrMsgIdAndTxt("QUIC:dimensions",
			      "The regularization parameter is not a scalar\n"
		"              or a matching matrix.");
	}
	Lambda = mxGetPr(prhs[argIdx]);
    }
    argIdx++;

    uint32_t pathLen = 1;
    double* path = NULL;
    if (mode[0] == 'p') {
	if (!mxIsDouble(prhs[argIdx]))
	    mexErrMsgIdAndTxt("QUIC:type",
			      "Expected a double matrix. (Arg. %d)",
		argIdx + 1);
	pathLen = (uint32_t) mxGetN(prhs[argIdx]);
	path = mxGetPr(prhs[argIdx]);
	if (pathLen <= 0) {
	    mexErrMsgIdAndTxt("QUIC:dimensions",
			      "Please specify the path scaling values.");
	}
	argIdx++;
    }

    double tol = 1e-6;
    if (nrhs > argIdx) {
	if (!mxIsDouble(prhs[argIdx]))
	    mexErrMsgIdAndTxt("QUIC:type",
			      "Expected a double scalar. (Arg. %d)",
		argIdx + 1);
	tol = mxGetScalar(prhs[argIdx]);
	if (tol < 0) {
	    mexErrMsgIdAndTxt("QUIC:tol",
			      "Negative tolerance value.");
	}
	argIdx++;
    }

    int32_t msg = 0;    
    if (nrhs > argIdx) {
	msg = (uint32_t) mxGetScalar(prhs[argIdx]);
	argIdx++;
    }
    
    // Maximum number of Newton ierations (whole matrix update):
    uint32_t maxIter = 1000;
    if (nrhs > argIdx) {
	maxIter = (uint32_t) mxGetScalar(prhs[argIdx]);
	argIdx++;
    }

    double* X0 = NULL;
    double* W0 = NULL;
    if (nrhs > argIdx) {
	if (!mxIsDouble(prhs[argIdx]))
	    mexErrMsgIdAndTxt("QUIC:type",
			      "Expected a double matrix. (Arg. %d)",
			      argIdx + 1);
	if (p != mxGetM(prhs[argIdx]) || p != mxGetN(prhs[argIdx]))
	    mexErrMsgIdAndTxt("QUIC:dimensions",
			      "Matrix dimensions should match.\n"
	         "             Maybe incorrect mode is specified?");
	X0 = mxGetPr(prhs[argIdx]);
	argIdx++;
	if (nrhs == argIdx)
	    mexErrMsgIdAndTxt("QUIC:initializations",
			      "Please specify both the initial estimate\n"
	         "             and the inverse.\n"
	         "             Maybe incorrect mode is specified?");
	if (!mxIsDouble(prhs[argIdx]))
	    mexErrMsgIdAndTxt("QUIC:type",
			      "Expected a double matrix. (Arg. %d)",
			      argIdx + 1);
	if (p != mxGetM(prhs[argIdx]) || p != mxGetN(prhs[argIdx])) {
	    mexErrMsgIdAndTxt("QUIC:dimensions",
			      "Matrix dimensions should match.\n"
	         "             Maybe incorrect mode is specified?");
	}
	W0 = mxGetPr(prhs[argIdx]);
	argIdx++;
    }

    double* X = NULL;
    double* W = NULL;
    if (mode[0] == 'p') {
	mwSize dims[] = {int32_t(p), int32_t(p), int32_t(pathLen)};
	mxArray* tmp = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	X = (double*) mxGetPr(tmp);
	if (nlhs > 0)
	    plhs[0] = tmp;
	tmp = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	W = (double*) mxGetPr(tmp);
	if (nlhs > 1)
	    plhs[1] = tmp;
    } else {
	mwSize dims[] = {int32_t(p), int32_t(p)};
	mxArray* tmp = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
	X = (double*) mxGetPr(tmp);
	if (nlhs > 0)
	    plhs[0] = tmp;
	tmp = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
	W = (double*) mxGetPr(tmp);
	if (nlhs > 1)
	    plhs[1] = tmp;	
    }
    if (X0 != NULL) {
	memcpy(X, X0, sizeof(double)*p*p);
	memcpy(W, W0, sizeof(double)*p*p);
    } else {
	memset(X, 0, sizeof(double)*p*p);
	memset(W, 0, sizeof(double)*p*p);
	for (unsigned long i = 0; i < p*p; i += (p+1)) {
	    X[i] = 1.0;
	    W[i] = 1.0;
	}
    }
    double* opt = NULL;
    double* dGap = NULL;
    double* cputime = NULL;
    uint32_t* iter = NULL;
    unsigned long optSize = 1;
    unsigned long iterSize = 1;
    if (mode[0] == 'p') {
	optSize = pathLen;
	iterSize = pathLen;
    } else if (mode[0] == 't')
	optSize = maxIter;
    if (nlhs > 2) {
	plhs[2] = mxCreateDoubleMatrix(optSize, 1, mxREAL);
	opt = mxGetPr(plhs[2]);
    }
    if (nlhs > 3) {
	plhs[3] = mxCreateDoubleMatrix(optSize, 1, mxREAL);
	cputime = mxGetPr(plhs[3]);
    }
    if (nlhs > 4) {
	mwSize dims[] = {int32_t(iterSize)};
	plhs[4] = mxCreateNumericArray(1, dims, mxINT32_CLASS, mxREAL);
	iter = (uint32_t*) mxGetData(plhs[4]);
    }
    if (nlhs > 5) {
	plhs[5] = mxCreateDoubleMatrix(optSize, 1, mxREAL);
	dGap = mxGetPr(plhs[5]);
    }
    QUIC(mode[0], p, S, Lambda, pathLen, path, tol, msg, maxIter, X, W,
	 opt, cputime, iter, dGap);
    if (LambdaAlloc)
	free(Lambda);
}
