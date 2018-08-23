//function [val] = AKVF( A2, imX)
//change to imX1imX2imX3 A2, imX
#include <mex.h>
#include <iostream>


void mexFunction(int nlhs, mxArray *plhs[],	int nrhs, const mxArray*prhs[])
{
  // jacobi(HNZ, SC, Hi, Hj) = HNZS 
   const size_t n = mxGetNumberOfElements(prhs[0]);
   if(mxGetNumberOfElements(prhs[2])!=n || mxGetNumberOfElements(prhs[3])!=n)
	   mexErrMsgTxt("bad input, inconsistent dimensions");
  
   const double * HNZ = (double*)mxGetData(prhs[0]);
   const double * Sc  = (double*)mxGetData(prhs[1]);
   const size_t* pidi = (size_t*)mxGetData(prhs[2]);
   const size_t* pidj = (size_t*)mxGetData(prhs[3]);
   
   plhs[0] = mxCreateDoubleMatrix(n,1, mxREAL);
   double * Val = (double*) mxGetData(plhs[0]);
   
   #pragma omp parallel for shared(Val, HNZ,Sc, pidi, pidj)
   for (int i = 0; i<n; i++){
	   Val[i] = HNZ[i] * Sc[pidi[i]-1] * Sc[pidj[i]-1];
   }
   
}
