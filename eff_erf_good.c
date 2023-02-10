/* ------------------------------------------------------------------------------- */
/* Added routines for calculating erf, (erf x)/x, and derivatives                  */
/* Originally from Hoops3D floating orbital code.  -- JTS 7/23/03                  */
/*                                                                                 */
/* Chebyshev routine adapted from Numerical Recipes 5.8                            */
/* Polynomial coefficients and bilinear mapping from                               */
/* Schonfelder, J.L. Math. Comp. 32, 144, (1978), 1232-40.                         */

/* Chebechev polynomial coefficients. */

#include "eff_erf.h"
#include <math.h>

#define ERF_TERMS1 12
#define ERF_TERMS2 7
#define DERF_TERMS 13

static double E1[] = 
{
	1.483110564084803581889448079057,
	-3.01071073386594942470731046311E-1,
	  6.8994830689831566246603180718E-2,
	 -1.3916271264722187682546525687E-2,
	   2.420799522433463662891678239E-3,
	   -3.65863968584808644649382577E-4,
	     4.8620984432319048282887568E-5,
		 -5.749256558035684835054215E-6,
		   6.11324357843476469706758E-7,
		   -5.8991015312958434390846E-8,
		     5.207009092068648240455E-9,
		    -4.23297587996554326810E-10,
		      3.1881135066491749748E-11,
		      -2.236155018832684273E-12,
		        1.46732984799108492E-13,
		         -9.044001985381747E-15,
			       5.25481371547092E-16,
			   	   -2.8874261222849E-17,
					 1.504785187558E-18,
					  -7.4572892821E-20,
					    3.522563810E-21,
					    -1.58944644E-22,
					       6.864365E-24,
				   		   -2.84257E-25,
							 1.1306E-26,
							  -4.33E-28,
							    1.6E-29,
						       -1.0E-30
};

static double E2[] = 
{
     1.077977852072383151168335910348,
	 -2.6559890409148673372146500904E-2,
	  -1.487073146698099509605046333E-3,
	   -1.38040145414143859607708920E-4,
	    -1.1280303332287491498507366E-5,
		 -1.172869842743725224053739E-6,
		  -1.03476150393304615537382E-7,
		   -1.1899114085892438254447E-8,
		    -1.016222544989498640476E-9,
		    -1.37895716146965692169E-10,
		      -9.369613033737303335E-12,
			  -1.918809583959525349E-12,
			    -3.7573017201993707E-14,
			    -3.7053726026983357E-14,
				  2.627565423490371E-15,
				 -1.121322876437933E-15,
				   1.84136028922538E-16,
				   -4.9130256574886E-17,
				    1.0704455167373E-17,
				    -2.671893662405E-18,
				      6.49326867976E-19,
				     -1.65399353183E-19,
				       4.2605626604E-20,
					  -1.1255840765E-20,
					    3.025617448E-21,
					    -8.29042146E-22,
					     2.31049558E-22,
		    		     -6.5469511E-23,
						  1.8842314E-23,
						  -5.504341E-24,
						   1.630950E-24,
						   -4.89860E-25,
						    1.49054E-25,
						    -4.5922E-26,
						     1.4318E-26,
						     -4.516E-27,
						      1.440E-27,
						      -4.64E-28,
						       1.51E-28,
							   -5.0E-29,
							    1.7E-29,
							   -6.0E-30,
							    2.0E-30,
							   -1.0E-30
};

static double DE1[] = 
{
	-0.689379974848418501361491576718,
	0.295939056851161774752959335568,
	-0.087237828075228616420029484096,
	0.019959734091835509766546612696,
	-0.003740200486895490324750329974,
	0.000593337912367800463413186784,
	-0.000081560801047403878256504204,
	9.886099179971884018535968E-6,
	-1.071209234904290565745194E-6,
	1.0490945447626050322784E-7,
	-9.370959271038746709966E-9,
	7.6927263488753841874E-10,
	-5.8412335114551520146E-11,
	4.125393291736424788E-12,
	-2.72304624901729048E-13,
	1.6869717361387012E-14,
	-9.84565340276638E-16,
	5.4313471880068E-17,
	-2.840458699772E-18,
	1.4120512798E-19,
	-6.688772574E-21,
	3.0257558E-22,
	-1.3097526E-23,
	5.4352E-25,
	-2.1704E-26,
	8.32E-28,
	-5.4E-29
};

static double DE2[] = 
{
0.717710208167480928473053690384
-0.379868973985143305103199928808
0.125832094465157378967135019248
-0.030917661684228839423081992424
0.006073689914144320367855343072
-0.000996057789064916825079352632
0.000140310790466315733723475232
-0.000017328176496070286001302184
1.90540194670935746397168e-6
-1.8882873760163694937908e-7
1.703176613666840587056e-8
-1.40955218086201517976e-9
1.0776816914256065828e-10
-7.656138112778696256e-12
5.07943557413613792e-13
-3.1608615530282912e-14
1.852036572003432e-15
-1.02524641430496e-16
5.37852808112e-18
-2.68128238704e-19
1.273321788e-20
-5.77335744e-22
2.504352e-23
-1.0446e-24
4.16e-26
-2.808e-27
};

double erf(double x)
{
	int i;
	double b0, b1, b2;

	double xp;   /* transformed x times two. */
	double x2;   /* x squared. */

	if (x < 2.0)
	{
		/* erf(x) = x * y(t)    */
		/* t = 2 * (x/2)^2 - 1. */
		xp = x * x - 2;   /* first change of variables */
		b1 = 0.0; b0 = 0.0;
		for (i = ERF_TERMS1; i >= 0; i--)
		{
			b2 = b1;
			b1 = b0;
			b0 = xp * b1 - b2 + E1[i];
		}
		return x * (b0 - b2) / 2.0;
	}
	else
	{
		/* erf(x) = 1 - exp(-x^2)/x * y(t) */
		/* t = (10.5 - x^2) / (2.5 + x^2)  */
		x2 = x * x;
		xp = 2 * (10.5 - x2) / (2.5 + x2);   /* first change of variables. */
		b1 = 0.0; b0 = 0.0;
		for (i = ERF_TERMS2; i >= 0; i--)
		{
			b2 = b1;
			b1 = b0;
			b0 = xp * b1 - b2 + E2[i];
		}
		return 1.0 - exp(-x2) * (b0 - b2) / (2 * x);
	}
}
		
double erfoverx(double x)
{
	int i;
	double b0, b1, b2;

	double xp;   /* transformed x times two. */
	double x2;   /* x squared. */

	if (x < 2.0)
	{
		/* erf(x) = x * y(t)     */
		/* t = 2 * (x/2)^2 - 1.  */
		xp = x * x - 2;   /* first change of variables */
		b1 = 0.0; b0 = 0.0;
		for (i = ERF_TERMS1; i >= 0; i--)
		{
			b2 = b1;
			b1 = b0;
			b0 = xp * b1 - b2 + E1[i];
		}
		return (b0 - b2) / 2.0;
	}
	else
	{
		/* erf(x) = 1 - exp(-x^2)/x * y(t) */
		/* t = (10.5 - x^2) / (2.5 + x^2)  */
		x2 = x * x;
		xp = 2 * (10.5 - x2) / (2.5 + x2);   /* first change of variables. */
		b1 = 0.0; b0 = 0.0;
		for (i = ERF_TERMS2; i >= 0; i--)
		{
			b2 = b1;
			b1 = b0;
			b0 = xp * b1 - b2 + E2[i];
		}
		return 1.0 / x - exp(-x2) * (b0 - b2) / (2 * x2);
	}
}
		
double derfoverx_short(double x)
{
	// Only valid for x < 2 (!!)
	int i;
	double b0, b1, b2;

	double xp;   // transformed x times two.

	// D(erf(x) / x) = y'(t) * 2 x 
	// t = 2 * (x/2)^2 - 1.
	xp = x * x - 2;   // first change of variables
	b1 = 0.0; b0 = 0.0;
	for (i = DERF_TERMS; i >= 0; i--)
	{
		b2 = b1;
		b1 = b0;
		b0 = xp * b1 - b2 + DE1[i];
	}
	return x * (b0 - b2) / 2;
}

double derfoverx(double x, double f)
{
        int i;
        double b0, b1, b2;

        double xp;   // transformed x times two.

        if (x < 2)
        {
          // D(erf(x) / x) = y'(t) * 2 x 
          // t = 2 * (x/2)^2 - 1.
          xp = x * x - 2;   // first change of variables
          b1 = 0.0; b0 = 0.0;
          for (i = DERF_TERMS; i >= 0; i--)
          {
                b2 = b1;
                b1 = b0;
                b0 = xp * b1 - b2 + DE1[i];
          }
          return x * (b0 - b2) / 2;
        }
        else
          return (1.12837916709551257389615890312 * exp(-x * x) - f) / x;
}

double d2erfoverx_short(double x)
{
	// Only valid for x < 2 (!!)
	int i;
	double b0, b1, b2;

	double xp;   // transformed x times two.

	// D(erf(x) / x) = y'(t) * 2 x 
	// t = 2 * (x/2)^2 - 1.
	xp = x * x - 2;   // first change of variables
	b1 = 0.0; b0 = 0.0;
	for (i = DERF_TERMS; i >= 0; i--)
	{
		b2 = b1;
		b1 = b0;
		b0 = xp * b1 - b2 + DE1[i];
	}
	return x * (b0 - b2) / 2;
}


