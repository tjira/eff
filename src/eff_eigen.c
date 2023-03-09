#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ROTATE(a,i,j,k,l) g=a[n*i+j];h=a[n*k+l];a[n*i+j]=g-s*(h+g*tau); a[n*k+l]=h+s*(g-h*tau);

void jacobi(double *matrix, double eigenvalue[], double *eigenvector, int n)
{
	int j, iq, ip, i;
	double tresh, theta, tau, t, sm, s, h, g, c;
	
	double *b, *z;
	b = (double *) malloc(sizeof(double) * n);
	z = (double *) malloc(sizeof(double) * n);

	for (ip = 0; ip < n; ip++)
	{
		for (iq = 0; iq < n; iq++)
		{
		  eigenvector[n * ip + iq] = 0.0;
		}
		eigenvector[(n + 1) * ip] = 1.0;
	}

	for (ip = 0; ip < n; ip++)
	{
		b[ip] = eigenvalue[ip] = matrix[(n + 1) * ip];
		z[ip] = 0.0;
	}

	for (i = 1; i <= 100; i++)
	{
		sm = 0.0;
		for (ip = 0; ip < n - 1; ip++)
		{
			for (iq = ip + 1; iq < n; iq++)
			{
				sm += fabs(matrix[n * ip + iq]);
			}
		}
		if (sm == 0.0)
		{
			break;
		}
		if (i < 4)
		{
			tresh = 0.2 * sm / (n * n);
		}
		else
		{
			tresh = 0.0;
		}

		for (ip = 0; ip < n - 1; ip++)
		{
			for (iq = ip + 1; iq < n; iq++)
			{
				g = 100.0 * fabs(matrix[n * ip + iq]);
				if (i > 4 && (fabs(eigenvalue[ip])+g) == fabs(eigenvalue[ip])
					&& (fabs(eigenvalue[iq])+g) == fabs(eigenvalue[iq]))
				{
					matrix[n * ip + iq] = 0.0;
				}
				else if (fabs(matrix[n * ip + iq]) > tresh)
				{
					h = eigenvalue[iq] - eigenvalue[ip];
					if ((fabs(h) + g) == fabs(h))
					{
						t = (matrix[n * ip + iq]) / h;
					}
					else
					{
						theta = 0.5 * h / (matrix[n * ip + iq]);
						t = 1.0/(fabs(theta) + sqrt(1.0 + theta * theta));
						if (theta < 0.0) t = -t;
					}
					c = 1.0 / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * matrix[n * ip + iq];
					z[ip] -= h;
					z[iq] += h;
					eigenvalue[ip] -= h;
					eigenvalue[iq] += h;
					matrix[n * ip + iq] = 0.0;
					for (j = 0; j < ip; j++)
					{
						ROTATE(matrix,j,ip,j,iq)
					}
					for (j = ip + 1; j < iq; j++)
					{
						ROTATE(matrix,ip,j,j,iq)
					}
					for (j = iq + 1; j < n; j++)
					{
						ROTATE(matrix,ip,j,iq,j)
					}
					for (j = 0; j < n; j++)
					{
						ROTATE(eigenvector,j,ip,j,iq)
					}
				}
			}
		}
		for (ip = 0; ip < n; ip++)
		{
			b[ip] += z[ip];
			eigenvalue[ip] = b[ip];
			z[ip] = 0.0;
		}
	}

	// Sorts eigenvalues and eigenvectors.  Based on eigsrt Numerical Recipe's
	// routine.

	int k;
	for (i = 0; i < n - 1; i++)
	{
		double p = eigenvalue[k = i];
		for (j = i; j < n; j++)
			if (eigenvalue[j] >= p) p = eigenvalue[k = j];
		if (k != i)
		{
			eigenvalue[k] = eigenvalue[i];
			eigenvalue[i] = p;
			for (j = 0; j < n; j++)
			{
				p = eigenvector[n * j + i];
				eigenvector[n * j + i] = eigenvector[n * j + k];
				eigenvector[n * j + k] = p;
			}
		}
	}
}

