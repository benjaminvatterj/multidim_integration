/*solver/fixed_point.c
* Solves a fixed point equation using
* sequence acceleration.
*
* Author: Benjamin Vatter j.
* email : benjaminvatterj@gmail.com
* date: 15 August 2015
*/

#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include "fixed_point.h"

/*
* Finds the fixed point of a function f
* Translated from SciPy's fixed_point function
*/
void fixed_point(int n, ff_function *f, double *x0,
                 double xtol, int maxiter, double *out)
{
	//out = malloc(sizeof(double)*n);
	double *p = malloc(sizeof(double)*n);
	double *p1, *p2;
	double d;
	int pass;
	int i, j;

	for (i=0; i<n; i++)
		out[i] = x0[i];

	for (i=0; i<maxiter; i++)
	{
		p1 = malloc(sizeof(double)*n);
		p2 = malloc(sizeof(double)*n);
		f->function(out, f->params, p1);
		f->function(p1, f->params, p2);

		pass = 1;
		for (j=0; j<n; j++)
		{
			d = p2[j] - 2.0*p1[j] + out[j];
			if (d==0) {
				p[j] = p2[j];
			} else {
				p[j] = out[j] - pow(p1[j] - out[j], 2.0) / d;
			}

			if (out[j] == 0 && fabs(p[j]) > xtol){
				pass = 0;
			}
			if (out[j] != 0 && (p[j] - out[j])/out[j] > xtol){
				pass = 0;
			}
			out[j] = p[j];
		}
		if (pass == 1) {
			free(p1);
			free(p2);
			break;
		}
		if(i+1 >= maxiter){
			free(p1);
			free(p2);
			pass = -1;
			break;
		}
		free(p1);
		free(p2);
	}

	free(p);
	// Free garbage result in case of no convergence
	if (pass == -1) {
		free(out);
	}
}
