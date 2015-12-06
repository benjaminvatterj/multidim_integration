/*solver/integration.c
*
* Function for numerical integration
* using quadrature and Monte Carlo
* methods.
*
* Author: Benjamin Vatter.
* email: benjaminvatterj@gmail.com
* date: 15 August 2015
*/

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "integration.h"
#include "dbg.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>


typedef struct {
	int dims;
	int depth;
	double *lb;
	double *ub;
	double *real_x;
	double epsabs;
	double epsrel;
	size_t limit;
	gsl_integration_workspace * workspace;
	double * abserr;
	nquad_function *f;
} nquad_params;

void gsl_error (const char * reason, const char * file, int line, int gsl_errno)
{
  gsl_stream_printf ("NQUAD ERROR", file, line, reason);


  if (gsl_errno != GSL_EROUND) {
	  fflush (stderr);
	  fflush (stdout);
	  abort ();
  } else {
  	log_warn("NQUAD ERROR: rounding error detected. This result might be unreliable.");
  	fflush (stderr);
  	return;
  }
}

/*
SUbroutine for recursive multidimensional integration
*/
double nquad_f(double x, void * params)
{
	nquad_params *in_pars = (nquad_params *) params;
	double result;
	double abserr;
	int step = in_pars->depth;
	// Add value
	if (step>0) {
	 	in_pars->real_x[step-1] = x;
	}
	if (in_pars->dims == in_pars->depth){
		result = in_pars->f->function(in_pars->real_x, in_pars->dims, in_pars->f->params);
		return result;
	} else {
		// Consistency
		assert((step+1) <= in_pars->dims);
		// increase depth
		in_pars->depth += 1;
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(in_pars->limit);
		gsl_function F = {.function = &nquad_f, .params = in_pars};
		gsl_integration_qags(&F, in_pars->lb[step], in_pars->ub[step], in_pars->epsabs, in_pars->epsrel, in_pars->limit, w,&result, &abserr);
		//gsl_integration_qag(&F, in_pars->lb[step], in_pars->ub[step], in_pars->epsabs, in_pars->epsrel, in_pars->limit, 2, w,&result, &abserr);
		in_pars->abserr[step] = abserr;
		gsl_integration_workspace_free(w);
		in_pars->depth -= 1;
		return result;
	}
}

/*
An extension of gsl_integration_qags
to n dimensions. The first parameter is
the number or dimensions. the signature of
f must be
f(double *x, void * params)

Parameters:
n - dimensions
f - function and parameters
lb - lower bounds to integrals
ub - upper bounds to integrals
epsabs - desired absolute error in each integral (only in nested-quadrature)
epsrel - desired relative error in each integral (only in nested-quadrature)
limit - max number of iterations
result - on out contains result
abserr - on out contains estimated error
*/
void nquad(int n, nquad_function *f,  double *lb, double *ub,
           double epsabs, double epsrel, size_t limit, double * result, double *abserr)
{
	assert(n>=1);
	nquad_params *new_params = malloc(sizeof(nquad_params));
	new_params->dims = n;
	new_params->depth = 0;
	new_params->lb = lb;
	new_params->ub = ub;
	new_params->epsabs = epsabs;
	new_params->epsrel = epsrel;
	new_params->limit = limit;
	new_params->abserr = malloc(sizeof(double) * n);
	new_params->f = f;
	new_params->real_x = malloc(sizeof(double) * n);

	// Set error handle
	gsl_error_handler_t *old_handler = gsl_set_error_handler(&gsl_error);

	// integrate
	result[0] = nquad_f(0.0, new_params);
	abserr[0] = 1.0;
	int i;
	double err = 0.0;
	for (i = 0; i<n; i++){
		err = (err < new_params->abserr[i])? new_params->abserr[i]: err;
	}
	*abserr = err;

	// release real x
	free(new_params->real_x);
	free(new_params->abserr);
	free(new_params);
	// recover error handle
	gsl_set_error_handler(old_handler);
}


/* Monte Carlo integration routine
* It has the same signature as the nquad function
* but it doesn't really use the epsrel value. It is
* done so to allow interchanging those two functions
* during runtime.
*
* This is based on the Monte Carlo plain routine found
* in GSL, written by Michael Booth
*/
void mcint(int n, nquad_function *f, double *lb, double *ub, double epsabs,
           double epsrel, size_t limit, double *result, double *abserr)
{
	double vol=1.0, m = 0, q = 0;
	double *x = malloc(sizeof(double)*n);
	size_t i, j;
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup ();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	/* Sanity check */
	for(i=0; i<n; i++) {
		if(lb[i] >= ub[i]) {
			GSL_ERROR("Lower bound is higher than upper bound. Null feasible region.", GSL_EINVAL);
		}
		if(ub[i] - lb[i] > GSL_DBL_MAX) {
			GSL_ERROR("Range of integration is too large, please rescale",
                     GSL_EINVAL);
		}
		vol *= ub[i] - lb[i];
	}

	j=0;

	do
	{
		for (i=0; i < n; i++) {
			x[i] = lb[i] + gsl_rng_uniform_pos(r) * (ub[i] - lb[i]);
		}
		double fval = (f->function)(x, n, f->params);
		double d = fval - m;
		m += d / (n + 1.0);
		q += d * d * (n/(n +1.0));

		if (j < 2){ *abserr = GSL_POSINF;}
		else {
			*abserr = vol * sqrt(q / (limit * (limit - 1.0)));
		}
		j++;
	} while((*abserr > epsabs) && j < limit);

	*result = vol * m;
	if (*abserr > epsabs){
		log_warn("Monte Carlo integration max evaluation %d was insufficient"
		"to achieve the desired absolute error of %e", (int)limit, epsabs);
	}

	free(x);
	gsl_rng_free (r);
}



