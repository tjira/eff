double dmin(double a, double b);
double dmax(double a, double b);
int conmin(int n, double *x, double *f, double *g,
           int *ifun, int *iter, double eps, int mxfun, double *w,
           int iout, int mdim, double acc, int nmeth, void (*calcfg)(double *, double *, double *),
           void (*update_func)(int, int, double, double));
