#include "eff_shanno.h"
#include <math.h>
#include <stdio.h>

/*
int main()
{ 
  double x[4], g[4], w[22], eps, f, acc;
  acc = 1e-20;
  eps = .00001;

  // nmeth = 0 for CG 
  //       = 1 for newton 
  int nmeth = 0;
  x[0] = -3.0;
  x[1] = -1.0;
  x[2] = -3.0;
  x[3] = -1.0;

  int ifun, iter;
  conmin(4, x, &f, g, &ifun, &iter, eps, 300, w, 1, 22, acc, nmeth, &calcfg, &update_output);

  printf("final f = %e, iterations = %i, function calls = %i\n", f, iter, ifun);
  printf("final x.\n");
  for (int i = 0; i < 4; i++)
    printf("%e\n", x[i]);
  printf("final g.\n");
  for (int i = 0; i < 4; i++)
    printf("%e\n", g[i]);
}

void update_output(int num_iterations, int num_fevals, double f_val, double g_sq)
{
  printf("\t%i\t%i\t%e\t%e\n", num_iterations, num_fevals, f_val, g_sq);
}

void calcfg(double *x, double *f, double *g)
{
  double a, b;
  a = x[1] - x[0] * x[0];
  b = x[3] - x[2] * x[2];
  *f = 100 * a * a + (1.0 - x[0]) * (1.0 - x[0])
       + 90 * b * b + (1.0 - x[2]) * (1.0 - x[2])
       + 10.1 * ((x[1] - 1) * (x[1] - 1) + (x[3] - 1) * (x[3] - 1))
       + 19.8 * (x[1] - 1) * (x[3] - 1);
  g[0] = -2 * (200 * x[0] * a + 1.0 - x[0]);
  g[1] = 2 * (100 * a + 10.1 * (x[1] - 1.0) + 9.9 * (x[3] - 1.0));
  g[2] = -2 * (180 * x[2] * b + 1.0 - x[2]);
  g[3] = 2 * (90 * b + 10.1 * (x[3] - 1) + 9.9 * (x[1] - 1.0));
  return;
}
*/

double dmin(double a, double b)
{
  if (a < b) return a; else return b;
}

double dmax(double a, double b)
{
  if (a > b) return a; else return b;
}

/* Changes from fortran version:
   -- NFLAG returned as return value, not entered as argument.
   -- IDEV not included, defaults to std. output.
   -- A function to be called every IOUT iterations is included here.
   -- arrays made into C-arrays that are indexed from zero, not one.
 
C PURPOSE:    SUBROUTINE CONMIN MINIMIZES AN UNCONSTRAINED NONLINEAR
C             SCALAR VALUED FUNCTION OF A VECTOR VARIABLE X
C             EITHER BY THE BFGS VARIABLE METRIC ALGORITHM OR BY A
C             BEALE RESTARTED CONJUGATE GRADIENT ALGORITHM.
C
C USAGE:      CALL CONMIN(N,X,F,G,IFUN,ITER,EPS,NFLAG,MXFUN,W,
C             IOUT,MDIM,IDEV,ACC,NMETH)
C
C PARAMETERS: N      THE NUMBER OF VARIABLES IN THE FUNCTION TO
C                    BE MINIMIZED.
C             X      THE VECTOR CONTAINING THE CURRENT ESTIMATE TO
C                    THE MINIMIZER. ON ENTRY TO CONMIN,X MUST CONTAIN
C                    AN INITIAL ESTIMATE SUPPLIED BY THE USER.
C                    ON EXITING,X WILL HOLD THE BEST ESTIMATE TO THE
C                    MINIMIZER OBTAINED BY CONMIN. X MUST BE DOUBLE
C                    PRECISIONED AND DIMENSIONED N.
C             F      ON EXITING FROM CONMIN,F WILL CONTAIN THE LOWEST
C                    VALUE OF THE OBJECT FUNCTION OBTAINED.
C                    F IS DOUBLE PRECISIONED.
C             G      ON EXITING FROM CONMIN,G WILL CONTAIN THE
C                    ELEMENTS OF THE GRADIENT OF F EVALUATED AT THE
C                    POINT CONTAINED IN X. G MUST BE DOUBLE
C                    PRECISIONED AND DIMENSIONED N.
C             IFUN   UPON EXITING FROM CONMIN,IFUN CONTAINS THE
C                    NUMBER OF TIMES THE FUNCTION AND GRADIENT
C                    HAVE BEEN EVALUATED.
C             ITER   UPON EXITING FROM CONMIN,ITER CONTAINS THE
C                    TOTAL NUMBER OF SEARCH DIRECTIONS CALCULATED
C                    TO OBTAIN THE CURRENT ESTIMATE TO THE MINIZER.
C             EPS    EPS IS THE USER SUPPLIED CONVERGENCE PARAMETER.
C                    CONVERGENCE OCCURS WHEN THE NORM OF THE GRADIENT
C                    IS LESS THAN OR EQUAL TO EPS TIMES THE MAXIMUM
C                    OF ONE AND THE NORM OF THE VECTOR X. EPS
C                    MUST BE DOUBLE PRECISIONED.
C             NFLAG  UPON EXITING FROM CONMIN,NFLAG STATES WHICH
C                    CONDITION CAUSED THE EXIT.
C                    IF NFLAG=0, THE ALGORITHM HAS CONVERGED.
C                    IF NFLAG=1, THE MAXIMUM NUMBER OF FUNCTION
C                       EVALUATIONS HAVE BEEN USED.
C                    IF NFLAG=2, THE LINEAR SEARCH HAS FAILED TO
C                       IMPROVE THE FUNCTION VALUE. THIS IS THE
C                       USUAL EXIT IF EITHER THE FUNCTION OR THE
C                       GRADIENT IS INCORRECTLY CODED.
C                    IF NFLAG=3, THE SEARCH VECTOR WAS NOT
C                       A DESCENT DIRECTION. THIS CAN ONLY BE CAUSED
C                       BY ROUNDOFF,AND MAY SUGGEST THAT THE
C                       CONVERGENCE CRITERION IS TOO STRICT.
C             MXFUN  MXFUN IS THE USER SUPPLIED MAXIMUM NUMBER OF
C                    FUNCTION AND GRADIENT CALLS THAT CONMIN WILL
C                    BE ALLOWED TO MAKE.
C             W      W IS A VECTOR OF WORKING STORAGE.IF NMETH=0,
C                    W MUST BE DIMENSIONED 5*N+2. IF NMETH=1,
C                    W MUST BE DIMENSIONED N*(N+7)/2. IN BOTH CASES,
C                    W MUST BE DOUBLE PRECISIONED.
C             IOUT   IOUT IS A USER  SUPPLIED OUTPUT PARAMETER.
C                    IF IOUT = 0, THERE IS NO PRINTED OUTPUT FROM
C                    CONMIN. IF IOUT > 0,THE VALUE OF F AND THE
C                    NORM OF THE GRADIENT SQUARED,AS WELL AS ITER
C                    AND IFUN,ARE WRITTEN EVERY IOUT ITERATIONS.
C             MDIM   MDIM IS THE USER SUPPLIED DIMENSION OF THE
C                    VECTOR W. IF NMETH=0,MDIM=5*N+2. IF NMETH=1,
C                    MDIM=N*(N+7)/2.
C             IDEV   IDEV IS THE USER SUPPLIED NUMBER OF THE OUTPUT
C                    DEVICE ON WHICH OUTPUT IS TO BE WRITTEN WHEN
C                    IOUT>0.
C             ACC    ACC IS A USER SUPPLIED ESTIMATE OF MACHINE
C                    ACCURACY. A LINEAR SEARCH IS UNSUCCESSFULLY
C                    TERMINATED WHEN THE NORM OF THE STEP SIZE
C                    BECOMES SMALLER THAN ACC. IN PRACTICE,
C                    ACC=10.D-20 HAS PROVED SATISFACTORY. ACC IS
C                    DOUBLE PRECISIONED.
C             NMETH  NMETH IS THE USER SUPPLIED VARIABLE WHICH
C                    CHOOSES THE METHOD OF OPTIMIZATION. IF
C                    NMETH=0,A CONJUGATE GRADIENT METHOD IS
C                    USED. IF NMETH=1, THE BFGS METHOD IS USED.
C
C REMARKS:    IN ADDITION TO THE SPECIFIED VALUES IN THE ABOVE
C             ARGUMENT LIST, THE USER MUST SUPPLY A SUBROUTINE
C             CALCFG WHICH CALCULATES THE FUNCTION AND GRADIENT AT
C             X AND PLACES THEM IN F AND G(1),...,G(N) RESPECTIVELY.
C             THE SUBROUTINE MUST HAVE THE FORM:
C                    SUBROUTINE CALCFG(N,X,F,G)
C                    DOUBLE PRECISION X(N),G(N),F
C
C             AN EXAMPLE SUBROUTINE FOR THE ROSENBROCK FUNCTION IS:
C
C                    SUBROUTINE CALCFG(N,X,F,G)
C                    DOUBLE PRECISION X(N),G(N),F,T1,T2
C                    T1=X(2)-X(1)*X(1)
C                    T2=1.0-X(1)
C                    F=100.0*T1*T1+T2*T2
C                    G(1)=-400.0*T1*X(1)-2.0*T2
C                    G(2)=200.0*T1
C                    RETURN
C                    END
C

*/
int conmin(int n, double *x, double *f, double *g,
           int *ifun, int *iter, double eps, int mxfun, double *w,
           int iout, int mdim, double acc, int nmeth, void (*calcfg)(double *, double *, double *),
           void (*update_func)(int, int, double, double))
{
  double fp, fmin, alpha, at, ap, gsq, dg, dg1;
  double dp, step, dal, u1, u2, u3, u4;
  double xsq, rtst, dsqrt, dmin1, dmax1, dabs;
  int ioutk, nflag, nx, ng;
  int nry, nrd, ncons, ncons0, ncons1;
  int nxpi, ngpi, nrdpi, nrypi, ngpj;
  int ij, ii, i, j;
  int nrst;
  int rsw;

  /* INITIALIZE ITER,IFUN,NFLAG,AND IOUTK,WHICH COUNTS OUTPUT ITERATIONS. */
  *iter = 0;
  *ifun = 0;
  ioutk = 0;
  nflag = 0;

  /* SET PARAMETERS TO EXTRACT VECTORS FROM W.
  W(I) HOLDS THE SEARCH VECTOR,W(NX+I) HOLDS THE BEST CURRENT
  ESTIMATE TO THE MINIMIZER,AND W(NG+I) HOLDS THE GRADIENT
  AT THE BEST CURRENT ESTIMATE.
  */

  nx = n;
  ng = nx + n;

  /* TEST WHICH METHOD IS BEING USED.
  IF NMETH=0, W(NRY+I) HOLDS THE RESTART Y VECTOR AND
  W(NRD+I) HOLDS THE RESTART SEARCH VECTOR.
  */

  if (nmeth == 1) goto label_10;
  nry = ng + n;
  nrd = nry + n;
  ncons = 5 * n;
  ncons0 = ncons + 0;
  ncons1 = ncons + 1;
  goto label_20;

  /* IF NMETH=1,W(NCONS+I) HOLDS THE APPROXIMATE INVERSE HESSIAN. */
label_10:
  ncons = 3 * n;
  
  /* CALCULATE THE FUNCTION AND GRADIENT AT THE INITIAL
  POINT AND INITIALIZE NRST,WHICH IS USED TO DETERMINE
  WHETHER A BEALE RESTART IS BEING DONE. NRST=N MEANS THAT THIS
  ITERATION IS A RESTART ITERATION. INITIALIZE RSW,WHICH INDICATES
  THAT THE CURRENT SEARCH DIRECTION IS A GRADIENT DIRECTION.
  */
label_20:
  (*calcfg)(x, f, g);
  (*ifun)++;
  nrst=n;
  rsw = 1;

  /* CALCULATE THE INITIAL SEARCH DIRECTION , THE NORM OF X SQUARED,
  AND THE NORM OF G SQUARED. DG1 IS THE CURRENT DIRECTIONAL
  DERIVATIVE,WHILE XSQ AND GSQ ARE THE SQUARED NORMS.
  */
  dg1 = 0;
  xsq = 0;
  for (i = 0; i < n; i++)
  {
    w[i] = -g[i];
    xsq += x[i] * x[i];
    dg1 -= g[i] * g[i];
  }
  gsq = -dg1;

  /* TEST IF THE INITIAL POINT IS THE MINIMIZER. */
  if (gsq <= eps * eps * dmax(xsq, 1.0))
    return nflag;

  /* BEGIN THE MAJOR ITERATION LOOP. NCALLS IS USED TO GUARANTEE THAT
  AT LEAST TWO POINTS HAVE BEEN TRIED WHEN NMETH=0. FMIN IS THE
  CURRENT FUNCTION VALUE.
  */
label_40:
  fmin = *f;
  int ncalls = *ifun;

  /* IF OUTPUT IS DESIRED,TEST IF THIS IS THE CORRECT ITERATION
  AND IF SO, WRITE OUTPUT. 
  */

  if (iout == 0) goto label_60;
  if (ioutk != 0) goto label_50;
  (*update_func)(*iter, *ifun, fmin, gsq);

label_50:
  ioutk++;
  if (ioutk == iout) ioutk = 0;

  /* BEGIN LINEAR SEARCH. ALPHA IS THE STEPLENGTH.
  SET ALPHA TO THE NONRESTART CONJUGATE GRADIENT ALPHA.
  */
label_60:
  alpha *= dg / dg1;

  /* IF NMETH=1 OR A RESTART HAS BEEN PERFORMED, SET ALPHA=1.0. */
  if (nrst == 1 || nmeth == 1) alpha = 1.0;

  /* IF A GRADIENT DIRECTION IS USED, SET ALPHA=1.0/DSQRT(GSQ),
  WHICH SCALES THE INITIAL SEARCH VECTOR TO UNITY.
  */
  if (rsw) alpha = 1.0 / sqrt(gsq);

  /* THE LINEAR SEARCH FITS A CUBIC TO F AND DAL, THE FUNCTION AND ITS
  DERIVATIVE AT ALPHA, AND TO FP AND DP,THE FUNCTION
  AND DERIVATIVE AT THE PREVIOUS TRIAL POINT AP.
  INITIALIZE AP ,FP,AND DP.
  */
  ap = 0;
  fp = fmin;
  dp = dg1;

  /* SAVE THE CURRENT DERIVATIVE TO SCALE THE NEXT SEARCH VECTOR. */
  dg = dg1;

  /* UPDATE THE ITERATION. */
  (*iter)++;

  /* CALCULATE THE CURRENT STEPLENGTH  AND STORE THE CURRENT X AND G. */
  step = 0;
  for (i = 0; i < n; i++)
  {
    step += w[i] * w[i];
    int nxpi = nx + i;
    int ngpi = ng + i;
    w[nxpi] = x[i];
    w[ngpi] = g[i];
  }
  step = sqrt(step);

  /* BEGIN THE LINEAR SEARCH ITERATION.
  TEST FOR FAILURE OF THE LINEAR SEARCH.
  */

label_80:
  if (alpha * step > acc) goto label_90;

  /* TEST IF DIRECTION IS A GRADIENT DIRECTION. */
  if (!rsw) goto label_20;
  nflag = 2;
  return nflag;

  /* CALCULATE THE TRIAL POINT. */
label_90:
  for (i = 0; i < n; i++)
  {
    nxpi = nx + i;
    x[i] = w[nxpi] + alpha * w[i];
  }

  /* EVALUATE THE FUNCTION AT THE TRIAL POINT. */
  (*calcfg)(x, f, g);

  /* TEST IF THE MAXIMUM NUMBER OF FUNCTION CALLS HAVE BEEN USED. */
  (*ifun)++;
  if (*ifun <= mxfun) goto label_110;
  nflag = 1;
  return nflag;

  /* COMPUTE THE DERIVATIVE OF F AT ALPHA. */
label_110:
  dal = 0.0;
  for (i = 0; i < n; i++)
    dal += g[i] * w[i];

  /* TEST WHETHER THE NEW POINT HAS A NEGATIVE SLOPE BUT A HIGHER
  FUNCTION VALUE THAN ALPHA=0. IF THIS IS THE CASE,THE SEARCH
  HAS PASSED THROUGH A LOCAL MAX AND IS HEADING FOR A DISTANT LOCAL
  MINIMUM. 
  */
  if (*f > fmin && dal < 0) goto label_160;

  /* IF NOT, TEST WHETHER THE STEPLENGTH CRITERIA HAVE BEEN MET. */
  if (*f > (fmin+.0001*alpha*dg) || fabs(dal/dg) > .9) goto label_130;

  /* IF THEY HAVE BEEN MET, TEST IF TWO POINTS HAVE BEEN TRIED
  IF NMETH=0 AND IF THE TRUE LINE MINIMUM HAS NOT BEEN FOUND.
  */
  if ((*ifun-ncalls) <= 1 && fabs(dal/dg) > eps && nmeth == 0) goto label_130;
  goto label_170;

  /* A NEW POINT MUST BE TRIED. USE CUBIC INTERPOLATION TO FIND
  THE TRIAL POINT AT.
  */
label_130:
  u1 = dp + dal - 3.0 * (fp - *f) / (ap - alpha);
  u2 = u1 * u1 - dp * dal;
  if (u2 < 0) u2 = 0;
  u2 = sqrt(u2);
  at = alpha - (alpha - ap) * (dal + u2 - u1) / (dal - dp + 2 * u2);

  /* TEST WHETHER THE LINE MINIMUM HAS BEEN BRACKETED. */
  if ((dal/dp) > 0) goto label_140;

  /* THE MINIMUM HAS BEEN BRACKETED. TEST WHETHER THE TRIAL POINT LIES
  SUFFICIENTLY WITHIN THE BRACKETED INTERVAL.
  IF IT DOES NOT, CHOOSE AT AS THE MIDPOINT OF THE INTERVAL.
  */
  if (at < (1.01 * dmin(alpha, ap)) || at > 0.99 * dmax(alpha,ap))
    at = (alpha + ap) / 2.0;
  goto label_150;

  /* THE MINIMUM HAS NOT BEEN BRACKETED. TEST IF BOTH POINTS ARE
  GREATER THAN THE MINIMUM AND THE TRIAL POINT IS SUFFICIENTLY
  SMALLER THAN EITHER.
  */
label_140:
  if (dal > 0.0 && 0.0 < at && at < 0.99 * dmin(ap, alpha))
    goto label_150;

  /* TEST IF BOTH POINTS ARE LESS THAN THE MINIMUM AND THE TRIAL POINT
  IS SUFFICIENTLY LARGE. */
  if (dal <= 0.0 && at > 1.01 * dmax(ap, alpha))
    goto label_150;

  /* IF THE TRIAL POINT IS TOO SMALL,DOUBLE THE LARGEST PRIOR POINT. */
  if (dal <= 0) at = 2.0 * dmax(ap, alpha);

  /* IF THE TRIAL POINT IS TOO LARGE, HALVE THE SMALLEST PRIOR POINT. */
  if (dal > 0) at = dmin(ap, alpha) / 2.0;

  /* SET AP=ALPHA, ALPHA=AT,AND CONTINUE SEARCH. */
label_150:
  ap = alpha;
  fp = *f;
  dp = dal;
  alpha = at;
  goto label_80;

  /* A RELATIVE MAX HAS BEEN PASSED.REDUCE ALPHA AND RESTART THE SEARCH. */
label_160:
  alpha = alpha / 3.0;
  ap = 0;
  fp = fmin;
  dp = dg;
  goto label_80;

  /* THE LINE SEARCH HAS CONVERGED. TEST FOR CONVERGENCE OF THE ALGORITHM. */
label_170:
  gsq = 0.0;
  xsq = 0.0;
  for (i = 0; i < n; i++)
  {
    gsq += g[i] * g[i];
    xsq += x[i] * x[i];
  }
  if (gsq <= eps * eps * dmax(1.0, xsq)) return nflag;

  /* SEARCH CONTINUES. SET W(I)=ALPHA*W(I),THE FULL STEP VECTOR. */
  for (i = 0; i < n; i++)
    w[i] *= alpha;

  /* COMPUTE THE NEW SEARCH VECTOR. FIRST TEST WHETHER A
  CONJUGATE GRADIENT OR A VARIABLE METRIC VECTOR IS USED. 
  */
  if (nmeth == 1) goto label_330;

  /* CONJUGATE GRADIENT UPDATE SECTION.
  TEST IF A POWELL RESTART IS INDICATED.
  */
  rtst = 0;
  for (i = 0; i < n; i++)
  {
    ngpi = ng + i;
    rtst += g[i] * w[ngpi];
  }
  if (fabs(rtst/gsq) > 0.2) nrst = n;

  /* IF A RESTART IS INDICATED, SAVE THE CURRENT D AND Y
  AS THE BEALE RESTART VECTORS AND SAVE D'Y AND Y'Y
  IN W(NCONS+1) AND W(NCONS+2).
  */

  if (nrst != n) goto label_220;
  w[ncons + 0] = 0;
  w[ncons + 1] = 0;
  for (i = 0; i < n; i++)
  {
    nrdpi = nrd + i;
    nrypi = nry + i;
    ngpi = ng + i;
    w[nrypi] = g[i] - w[ngpi];
    w[nrdpi] = w[i];
    w[ncons0] += w[nrypi] * w[nrypi];
    w[ncons1] += w[i] * w[nrypi];
  }

  /* CALCULATE  THE RESTART HESSIAN TIMES THE CURRENT GRADIENT. */
label_220:
  u1 = 0.0;
  u2 = 0.0;
  for (i = 0; i < n; i++)
  {
    nrdpi = nrd + i;
    nrypi = nry + i;
    u1 -= w[nrdpi] * g[i] / w[ncons0];
    u2 += w[nrdpi] * g[i] * 2 / w[ncons1] - w[nrypi] * g[i] / w[ncons0];
  }
  u3 = w[ncons1] / w[ncons0];
  for (i = 0; i < n; i++)
  {
    nxpi = nx + i;
    nrdpi = nrd + i;
    nrypi = nry + i;
    w[nxpi] = -u3 * g[i] - u1 * w[nrypi] - u2 * w[nrdpi];
  }

  /* IF THIS IS A RESTART ITERATION,W(NX+I) CONTAINS THE NEW SEARCH VECTOR. */
  if (nrst == n) goto label_300;

  /* NOT A RESTART ITERATION. CALCULATE THE RESTART HESSIAN TIMES THE CURRENT Y. */
label_250:
  u1 = u2 = u3 = u4 = 0;
  for (i = 0; i < n; i++)
  {
    ngpi = ng + i;
    nrdpi = nrd + i;
    nrypi = nry + i;
    u1 = u1 - (g[i] - w[ngpi]) * w[nrdpi] / w[ncons0];
    u2 = u2 - (g[i] - w[ngpi]) * w[nrypi] / w[ncons0] + 2.0 * w[nrdpi] * (g[i] - w[ngpi]) / w[ncons1];
    u3 = u3 + w[i] * (g[i] - w[ngpi]);
  }
  step = 0;
  for (i = 0; i < n; i++)
  {
    ngpi = ng + i;
    nrdpi = nrd + i;
    nrypi = nry + i;
    step = (w[ncons1] / w[ncons0]) * (g[i] - w[ngpi]) + u1 * w[nrypi] + u2 * w[nrdpi];
    u4 = u4 + step * (g[i] - w[ngpi]);
    w[ngpi] = step;
  }

  /* CALCULATE THE DOUBLY UPDATED HESSIAN TIMES THE CURRENT
  GRADIENT TO OBTAIN THE SEARCH VECTOR. */
  u1 = 0.0;
  u2 = 0.0;
  for (i = 0; i < n; i++)
  {
    u1 = u1 - w[i] * g[i] / u3;
    ngpi = ng + i;
    u2 = u2 + (1.0 + u4 / u3) * w[i] * g[i] / u3 - w[ngpi] * g[i] / u3;
  }
  for (i = 0; i < n; i++)
  {
    ngpi = ng + i;
    nxpi = nx + i;
    w[nxpi] = w[nxpi] - u1 * w[ngpi] - u2 * w[i];
  }

  /* CALCULATE THE DERIVATIVE ALONG THE NEW SEARCH VECTOR. */
label_300:
  dg1 = 0;
  for (i = 0; i < n; i++)
  {
    nxpi = nx + i;
    w[i] = w[nxpi];
    dg1 = dg1 + w[i] * g[i];
  }

  /* IF THE NEW DIRECTION IS NOT A DESCENT DIRECTION, STOP. */
  if (dg1 > 0) goto label_320;

  /* UPDATE NRST TO ASSURE AT LEAST ONE RESTART EVERY N ITERATIONS. */
  if (nrst == n) nrst = 0;
  nrst++;
  rsw = 0;
  goto label_40;

label_320:
  nflag = 3;
  return nflag;

  /* A VARIABLE METRIC ALGORITM IS BEING USED. CALCULATE Y AND D'Y. */
label_330:
  u1 = 0.0;
  for (i = 0; i < n; i++)
  {
    ngpi = ng + i;
    w[ngpi] = g[i] - w[ngpi];
    u1 = u1 + w[i] * w[ngpi];
  }

  /* IF RSW=.TRUE.,SET UP THE INITIAL SCALED APPROXIMATE HESSIAN. */
  if (!rsw) goto label_380;

  /* CALCULATE Y'Y. */
  u2 = 0;
  for (i = 0; i < n; i++)
  {
    ngpi = ng + i;
    u2 = u2 + w[ngpi] * w[ngpi];
  }
  
  /* CALCULATE THE INITIAL HESSIAN AS H=(P'Y/Y'Y)*I AND THE INITIAL U2=Y'HY AND W(NX+I)=HY. */
  ij = 0;
  u3 = u1 / u2;
  for (i = 0; i < n; i++)
  {
    for (j = i; j < n; j++)
    {
      ncons0 = ncons + ij;
      w[ncons0] = 0.0;
      if (i == j) w[ncons0] = u3;
      ij++;
    }
    nxpi = nx + i;
    ngpi = ng + i;
    w[nxpi] = u3 * w[ngpi];
  }
  u2 = u3 * u2;
  goto label_430;

 /* CALCULATE W(NX+I)=HY AND U2=Y'HY. */
label_380:
  u2 = 0.0;
  for (i = 0; i < n; i++)
  {
    u3 = 0.0;
    ij = i;
    if (i == 0) goto label_400;
    ii = i;
    for (j = 0; j < ii; j++)
    {
      ngpj = ng + j;
      ncons0 = ncons + ij;
      u3 = u3 + w[ncons0] * w[ngpj];
      ij = ij + n - j - 1;
    }
label_400:
    for (j = i; j < n; j++)
    {
      ncons0 = ncons + ij;
      ngpj = ng + j;
      u3 = u3 + w[ncons0] * w[ngpj];
      ij++;
    }
    ngpi = ng + i;
    u2 = u2 + u3 * w[ngpi];
    nxpi = nx + i;
    w[nxpi] = u3;
  }

  /* CALCULATE THE UPDATED APPROXIMATE HESSIAN. */
label_430:
  u4 = 1.0 + u2 / u1;
  for (i = 0; i < n; i++)
  {
    nxpi = nx + i;
    ngpi = ng + i;
    w[ngpi] = u4 * w[i] - w[nxpi];
  }
  ij = 0;
  for (i = 0; i < n; i++)
  {
    nxpi = nx + i;
    u3 = w[i] / u1;
    u4 = w[nxpi] / u1;
    for (j = i; j < n; j++)
    {
      ncons0 = ncons + ij;
      ngpj = ng + j;
      w[ncons0] = w[ncons0] + u3 * w[ngpj] - u4 * w[j];
      ij++;
    }
  }

  /* CALCULATE THE NEW SEARCH DIRECTION W(I)=-HG AND ITS DERIVATIVE. */
  dg1 = 0.0;
  for (i = 0; i < n; i++)
  {
    u3 = 0.0;
    ij = i;
    if (i == 0) goto label_470;
    ii = i;
    for (j = 0; j < ii; j++)
    {
      ncons0 = ncons + ij;
      u3 = u3 - w[ncons0] * g[j];
      ij = ij + n - j - 1;
    }
  label_470:
    for (j = i; j < n; j++)
    {
      ncons0 = ncons + ij;
      u3 = u3 - w[ncons0] * g[j];
      ij++;
    }
    dg1 = dg1 + u3 * g[i];
    w[i] = u3;
  }

  /* TEST FOR A DOWNHILL DIRECTION. */
  if (dg1 > 0) goto label_320;
  rsw = 0;
  goto label_40;
}
