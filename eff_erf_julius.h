#ifndef ERF_H
#define ERF_H

/* ------------------------------------------------------------------------------- */
/* Added routines for calculating erf, (erf x)/x, and derivatives                  */
/* Originally from Hoops3D floating orbital code.  -- JTS 7/23/03                  */
/*                                                                                 */

double erfoverx(double x);
double erfoverx1(double x, double *df);
double erfoverx2(double x, double *df, double *df2);

double poly01(double x);
double poly02(double x);
double poly1(double x);
double poly2(double x);

#endif
