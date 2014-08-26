#include <R.h>
#include <Rmath.h>

/******************/
/* random numbers */
/******************/

void F77_SUB(rndstart)(void) 
{ GetRNGstate(); }

void F77_SUB(rndend)(void) 
{ PutRNGstate(); }

double F77_SUB(ggrunif)(double *a, double *b)
{return runif(*a, *b); }

double F77_SUB(ggrbinom)(double *n, double *p)
{return rbinom(*n, *p); }

/* normal */
double F77_SUB(normrnd)(void) 
{ return norm_rand(); }

double F77_SUB(ggrnorm)(double *mu, double *sigma) 
{ return rnorm(*mu, *sigma); }

double F77_SUB(ggpnorm)(double *x, double *mu, double *sigma,
			int *lower_tail, int *give_log)
{return pnorm( *x, *mu, *sigma, *lower_tail, *give_log); }

double F77_SUB(ggdnorm)(double *x, double *mu, double *sigma, int *give_log)
{return dnorm( *x, *mu, *sigma, *give_log); }

/* gamma */
double F77_SUB(ggrgam)(double *a, double *scale)
{return rgamma(*a, *scale); }

double F77_SUB(ggqgam)(double *p, double *alpha, double *scale, 
		       int *lower_tail, int *log_p)
{ return qgamma(*p, *alpha, *scale, *lower_tail, *log_p); }

/* exponential */
double F77_SUB(ggqexp)(double *p, double *scale, int *lower_tail, int *log_p)
{return qexp(*p, *scale, *lower_tail, *log_p); }



/*********************/
/* special functions */
/*********************/

/* gamma function  */
double F77_SUB(gggammafn)(double *x) 
{ return gammafn(*x); }

/* log of gamma function */
double F77_SUB(gglgammafn)(double *x) 
{ return lgammafn(*x); }

/* Bessel K  */
double F77_SUB(ggbesselk)(double *x, double *nu, double *expo)
{return bessel_k(*x,*nu,*expo);}
 

/* log(1+x) for |x| << 1  */
double F77_SUB(gglog1p)(double *x) 
{ return log1p(*x);}




/* Poisson */
double F77_SUB(ggpois)(double *l) 
{ return rpois(*l); }


/* Log gamma function */
double F77_SUB(gglngamma)(double *x) 
{ return lgamma(*x); }

/* Sample */
/*double F77_SUB(ggsample)(double *x, double *size) 
  { return sample(*x, *size); }*/


