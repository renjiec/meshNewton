//function [val] = AKVF( A2, imX)
//change to imX1imX2imX3 A2, imX
#include <mex.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>


using namespace Eigen;

void mexFunction(int nlhs, mxArray *plhs[],	int nrhs, const mxArray*prhs[])
{
   const size_t n = mxGetNumberOfElements(prhs[0]);
   const double *Apt = mxGetPr(prhs[0]);
   Map<const MatrixXd> A2(mxGetPr(prhs[0]), n, 1);
   Map<const MatrixXd> ImsX(mxGetPr(prhs[1]), 1, n);
   Map<const MatrixXd> ImsY(mxGetPi(prhs[1]), 1, n);
   const size_t* pidx = (size_t*)mxGetData(prhs[2]);
   
   plhs[0] = mxCreateDoubleMatrix(n,36, mxREAL);
   Map<MatrixXd> Val(mxGetPr(plhs[0]), n,36);
   
   #pragma omp parallel for shared(Val, ImsX, ImsY, pidx)
   for (int i = 0; i<n; i++){
	   size_t t1 = pidx[i*3 + 0] - 1;
	   size_t t2 = pidx[i*3 + 1] - 1;
	   size_t t3 = pidx[i*3 + 2] - 1;
	   double x12 = ImsX(t1) - ImsX(t2);
	   double x23 = ImsX(t2) - ImsX(t3);
	   double x31 = ImsX(t3) - ImsX(t1);
	   
	   double y12 = ImsY(t1) - ImsY(t2);
	   double y23 = ImsY(t2) - ImsY(t3);
	   double y31 = ImsY(t3) - ImsY(t1);

		double areas = abs ( y12*x31 - x12*y31 );
		areas = (areas* areas)/(A2(i,0)/2);
		areas = sqrt(0.5 * areas);
		x12 /=areas;
		y12 /=areas;
		x23 /=areas;
		y23 /=areas;
		x31 /=areas;
		y31 /=areas;

	   Val(i, 0) =  x23*x23 + 2*y23*y23;
	   Val(i, 1) = 2*x23*x23 + y23*y23;
	   Val(i, 2) = -x23*y23;                       
	   Val(i, 3) = x31*x31 + 2*y31*y31;
	   Val(i, 4) = 2*x31*x31 + y31*y31;
	   Val(i, 5) = -x31*y31;                       
	   Val(i, 6) = x12*x12 + 2*y12*y12;
	   Val(i, 7) = 2*x12*x12 + y12*y12;
	   Val(i, 8) = -x12*y12;                       
	   Val(i, 9) = x31*x23 + 2*y23*y31;  
	   Val(i, 10) = 2*x31*x23 + y23*y31;  
	   Val(i, 11) = -x31*y23;                        
	   Val(i, 12) = -x23*y31;                        
	   Val(i, 13) = x12*x23 + 2*y12*y23;
	   Val(i, 14) = 2*x12*x23 + y12*y23;
	   Val(i, 15) = -x12*y23;                       
	   Val(i, 16) =  -x23*y12;                       
	   Val(i, 17) = x12*x31 + 2*y12*y31;
	   Val(i, 18) = 2*x12*x31 + y12*y31;
	   Val(i, 19) = -x12*y31;
	   Val(i, 20) = -x31*y12;
	   Val(i, 21) = -x23*y23;
	   Val(i, 22) = -x31*y31;
	   Val(i, 23) = -x12*y12;
	   Val(i, 24) = x31*x23 + 2*y23*y31;
	   Val(i, 25) = 2*x31*x23 + y23*y31;
	   Val(i, 26) = -x31*y23;
	   Val(i, 27) = -x23*y31;
	   Val(i, 28) = x12*x23 + 2*y12*y23;
	   Val(i, 29) = 2*x12*x23 + y12*y23;
	   Val(i, 30) = -x12*y23;
	   Val(i, 31) = -x23*y12;
	   Val(i, 32) = x12*x31 + 2*y12*y31;
	   Val(i, 33) = 2*x12*x31 + y12*y31;
	   Val(i, 34) = -x12*y31;
	   Val(i, 35) = -x31*y12;
   }
   
}
