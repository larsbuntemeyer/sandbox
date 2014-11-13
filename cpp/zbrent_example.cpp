/****************************************************
*      Program to demonstrate the real domain       *
*               Zbrent subroutine                   *
* ------------------------------------------------- *
* Reference:  BORLAND MATHEMATICAL LIBRARY          *
*                                                   *
*                C++ version by J-P Moreau, Paris.  *
* ------------------------------------------------- *
* Example:    Find a real root of f(x)=(x+1)^5      *
*                                                   *
* SAMPLE RUN:                                       *
*                                                   *
*  Input interval (X1,X2):                          *
*                                                   *
*        X1 = -2                                    *
*        X2 =  0                                    *
*                                                   *
*  Convergence criterion: 1e-10                     *
*  Maximum number of iterations: 10                 *
*                                                   *
*  The estimated root is:                           *
*                                                   *
*        X = -1.000000                              *
*                                                   *
*  The associated Y value is Y = 0.000000           *
*                                                   *
*  The number of iterations was: 2                  *
*  The error code is: 0                             *
*                                                   *
****************************************************/
#include <stdio.h>
#include <math.h>

#define FALSE  0
#define TRUE   1

double  e,r,x1,x2,yr;
int     err,maxiter,k;


// Test function for Brent method
double AFunction(double x) {
  return (x+1)*(x+1)*(x+1)*(x+1)*(x+1);
}

// TRUE if x1*x2 negative
int RootBracketed(double x1,double x2) {
  int result;
  if ((x1 > 0 && x2 > 0) || (x1 < 0 && x2 < 0)) 
    result = FALSE;
  else
    result = TRUE;
  return result;
}

// returns the minimum of two real numbers
double Minimum(double x1,double x2) {
  double result;
  if (x1 < x2) result = x1;
  else result = x2;
  return result;
}

/****************************************************
*              Brent Method Function                *
* ------------------------------------------------- *
* The purpose is to find a real root of a real      *
* function f(x) using Brent method.                 *
*                                                   *
* INPUTS:  x1,x2     : interval of root             *
*          Tolerance : desired accuracy for root    *
*          maxIter   : maximum number of iterations *
*                                                   *
* OUTPUTS: The function returns the root value      *
*          ValueAtRoot : value of f(root)           *
*          niter    : number of done iterations     *
*          error    : =0, all OK                    *
*                   : =1, no root found in interval *
*                   : =2, no more iterations !      *
****************************************************/
double BrentRoots( double x1, double x2, double Tolerance,
                   int maxIterations,
		           double *valueAtRoot,
                   int *niter, int *error )  {

  double FPP = 1e-11, nearzero = 1e-20;

  double result, AA, BB, CC, DD, EE, FA, FB, FC, Tol1, PP, QQ, RR, SS, xm;
  int i, done;

  i = 0; done = FALSE;   error = 0;
  AA = x1;  BB = x2;  FA = AFunction(AA); FB = AFunction(BB);
  if (!(RootBracketed(FA,FB))) 
    *error = 1;
  else {
    FC = FB;
    do {
      if (!(RootBracketed(FC,FB))) {
        CC = AA; FC = FA; DD = BB - AA; EE = DD;
      }
      if (fabs(FC) < fabs(FB)) {
        AA = BB; BB = CC; CC = AA;
        FA = FB; FB = FC; FC = FA;
      }
      Tol1 = 2.0 * FPP * fabs(BB) + 0.5 * Tolerance;
      xm = 0.5 * (CC-BB);
      if ((fabs(xm) <= Tol1) || (fabs(FA) < nearzero)) {
        result = BB;
        done = TRUE;
        *valueAtRoot = AFunction(result);
      } // A root has been found
      else {
        if ((fabs(EE) >= Tol1) && (fabs(FA) > fabs(FB))) {
          SS = FB/ FA;
          if (fabs(AA - CC) < nearzero) {
            PP = 2.0 * xm * SS;
            QQ = 1.0 - SS;
          }
          else {
            QQ = FA/FC;
            RR = FB /FC;
            PP = SS * (2.0 * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1.0));
            QQ = (QQ - 1.0) * (RR - 1.0) * (SS - 1.0);
          }
          if (PP > nearzero) QQ = -QQ;
          PP = fabs(PP);
          if ((2.0 * PP) < Minimum(3.0*xm *QQ-fabs(Tol1 * QQ), fabs(EE * QQ))) {
            EE = DD;  DD = PP/QQ;
          }
          else {
            DD = xm;   EE = DD;
          }
        }
        else {
          DD = xm;
          EE = DD;
        }
        AA = BB;
        FA = FB;
        if (fabs(DD) > Tol1) 
          BB = BB + DD;
        else {
          if (xm > 0) BB = BB + fabs(Tol1);
          else BB = BB - fabs(Tol1);
        }
        FB = AFunction(BB);
        i++;
      }
	}  while ((!done) && (i < maxIterations));
    if (i >= maxIterations) *error = 2;
  }
  *niter = i;
  return result;
} // BrentRoots()


void main()  {
  printf("\n Input interval (X1,X2):\n\n");
  printf("       X1 = "); scanf("%lf",&x1);
  printf("       X2 = "); scanf("%lf",&x2);
  printf("\n Convergence criterion: "); scanf("%lf",&e);
  printf("\n Maximum number of iterations: "); scanf("%d",&maxiter);

  r = BrentRoots(x1,x2,e,maxiter,&yr,&k,&err);

  printf("\n\n The estimated root is:\n\n");
  printf("       X = %f\n\n", r);
  printf(" The associated Y value is Y = %f\n\n", yr);
  printf(" The number of iterations was: %d\n",k);  
  printf(" The error code is: %d\n\n",err);
}

// End of file zbrent.cpp
