#include <stdio.h>

#define stepredn	0.2
#define acctol		0.0001
#define reltest		10.0

#include <R_ext/Applic.h>

#define MATHLIB_STANDALONE
#include <Rmath.h>

// hack macro definitions - can be improved
#define Rf_error printf
#define Rprintf printf

#include <malloc.h>
#define R_alloc(a,b) malloc((a)*(b))

static double * vect(int n)
{
    return (double *)R_alloc(n, sizeof(double));
}
static double ** matrix(int nrh, int nch)
{
    int   i;
    double **m;

    m = (double **) R_alloc((nrh + 1), sizeof(double *));
    for (i = 0; i <= nrh; i++)
	m[i] = (double*) R_alloc((nch + 1), sizeof(double));
    return m;
}

static double ** Lmatrix(int n)
{
    int   i;
    double **m;

    m = (double **) R_alloc(n, sizeof(double *));
    for (i = 0; i < n; i++)
	m[i] = (double *) R_alloc((i + 1), sizeof(double));
    return m;
}

/*  BFGS variable-metric method, based on Pascal code
in J.C. Nash, `Compact Numerical Methods for Computers', 2nd edition,
converted by p2c then re-crafted by B.D. Ripley */

void
vmmin(int n0, double *b, double *Fmin, optimfn fminfn, optimgr fmingr,
      int maxit, int trace, int *mask,
      double abstol, double reltol, int nREPORT, void *ex,
      int *fncount, int *grcount, int *fail)
{
    Rboolean accpoint, enough;
    double *g, *t, *X, *c, **B;
    int   count, funcount, gradcount;
    double f, gradproj;
    int   i, j, ilast, iter = 0;
    double s, steplength;
    double D1, D2;
    int   n, *l;

    if (maxit <= 0) {
	*fail = 0;
	*Fmin = fminfn(n0, b, ex);
	*fncount = *grcount = 0;
	return;
    }

    /* if (nREPORT <= 0) */
    /* 	error(_("REPORT must be > 0 (method = \"BFGS\")")); */
    l = (int *) R_alloc(n0, sizeof(int));
    n = 0;
    for (i = 0; i < n0; i++) if (mask[i]) l[n++] = i;
    g = vect(n0);
    t = vect(n);
    X = vect(n);
    c = vect(n);
    B = Lmatrix(n);
    f = fminfn(n0, b, ex);
    /* if (!R_FINITE(f)) */
    /* 	error(_("initial value in 'vmmin' is not finite")); */
    if (trace) Rprintf("initial  value %f \n", f);
    *Fmin = f;
    funcount = gradcount = 1;
    fmingr(n0, b, g, ex);
    iter++;
    ilast = gradcount;

    do {
	if (ilast == gradcount) {
	    for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) B[i][j] = 0.0;
		B[i][i] = 1.0;
	    }
	}
	for (i = 0; i < n; i++) {
	    X[i] = b[l[i]];
	    c[i] = g[l[i]];
	}
	gradproj = 0.0;
	for (i = 0; i < n; i++) {
	    s = 0.0;
	    for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
	    for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
	    t[i] = s;
	    gradproj += s * g[l[i]];
	    if (trace) Rprintf("gradproj=%f\n",gradproj);
	    if (trace) Rprintf("t[0]=%f\n",t[0]);
	}

	if (gradproj < 0.0) {	/* search direction is downhill */
	    steplength = 1.0;
	    accpoint = FALSE;
	    do {
		count = 0;
		for (i = 0; i < n; i++) {
		    b[l[i]] = X[i] + steplength * t[i];
		    if (reltest + X[i] == reltest + b[l[i]]) /* no change */
			count++;
		}
		if (count < n) {
		    f = fminfn(n0, b, ex);
		    funcount++;
		    accpoint = R_FINITE(f) &&
			(f <= *Fmin + gradproj * steplength * acctol);
		    if (!accpoint) {
			steplength *= stepredn;
		    }
		}
	    } while (!(count == n || accpoint));
	    enough = (f > abstol) &&
		fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);
	    /* stop if value if small or if relative change is low */
	    if (!enough) {
		count = n;
		*Fmin = f;
	    }
	    if (count < n) {/* making progress */
		*Fmin = f;
		fmingr(n0, b, g, ex);
		gradcount++;
		iter++;
		D1 = 0.0;
		for (i = 0; i < n; i++) {
		    t[i] = steplength * t[i];
		    c[i] = g[l[i]] - c[i];
		    D1 += t[i] * c[i];
		}
		if (trace) Rprintf("D1 = %f \n", D1);
		if (D1 > 0) {
		    D2 = 0.0;
		    for (i = 0; i < n; i++) {
			s = 0.0;
			for (j = 0; j <= i; j++)
			    s += B[i][j] * c[j];
			for (j = i + 1; j < n; j++)
			    s += B[j][i] * c[j];
			X[i] = s;
			D2 += s * c[i];
		    }
		    if (trace) {
		      Rprintf("X[0] = %f \n", X[0]);
		      Rprintf("c[0] = %f \n", c[0]);
		      Rprintf("t[0] = %f \n", t[0]);
		      Rprintf("g[0] = %f \n", g[0]);
		      Rprintf("B[0][0] = %f \n", B[0][0]);
		      Rprintf("D2 = %f \n", D2);
		    }
		    D2 = 1.0 + D2 / D1;
		    for (i = 0; i < n; i++) {
			for (j = 0; j <= i; j++)
			    B[i][j] += (D2 * t[i] * t[j]
					- X[i] * t[j] - t[i] * X[j]) / D1;
		    }
		    if (trace) {
		      Rprintf("D2' = %f \n", D2);
		      Rprintf("B[0][0] = %f \n", B[0][0]);
		    }
		} else {	/* D1 < 0 */
		    ilast = gradcount;
		}
	    } else {	/* no progress */
		if (ilast < gradcount) {
		    count = 0;
		    ilast = gradcount;
		}
	    }
	} else {		/* uphill search */
	    count = 0;
	    if (ilast == gradcount) count = n;
	    else ilast = gradcount;
	    /* Resets unless has just been reset */
	}
	if (trace && (iter % nREPORT == 0))
	    Rprintf("iter%4d value %f\n", iter, f);
	if (iter >= maxit) break;
	if (gradcount - ilast > 2 * n)
	    ilast = gradcount;	/* periodic restart */
    } while (count != n || ilast != gradcount);
    if (trace) {
	Rprintf("final  value %f \n", *Fmin);
	if (iter < maxit) Rprintf("converged\n");
	else Rprintf("stopped after %i iterations\n", iter);
    }
    *fail = (iter < maxit) ? 0 : 1;
    *fncount = funcount;
    *grcount = gradcount;
}


/* typedef double optimfn(int, double *, void *); */
double fn1(int n, double * x, void * ignore) {
  return x[0]*x[0];
}
double fn2(int n, double * x, void * ignore) {
  return x[0]*x[0]+x[1]*x[1];
}
double rosenbrock(int n, double * par, void * ignore) {
  double a=1.0, b=100.0, x=par[0], y=par[1];
  return (a-x)*(a-x)+b*(y-x*x)*(y-x*x);
}

/* typedef void optimgr(int, double *, double *, void *); */
void gr1(int n, double * x, double *gr, void * ignore) {
  gr[0] = 2.0*x[0];
}
void gr2(int n, double * x, double *gr, void * ignore) {
  gr[0] = 2.0*x[0];
  gr[1] = 2.0*x[1];
}
void rosenbrock1(int n, double * par, double *gr, void * ignore) {
  double a=1.0, b=100.0, x=par[0], y=par[1];
  gr[0] = -4.0*b*x*y+4.0*b*x*x*x+2.0*x-2.0*a;
  gr[1] = 2.0*b*y-2.0*b*x*x;
}

int main() {

  /* int n0 = 1; */
  /* double b[1] = {1.0}; */
  /* int mask[1] = {1}; */
  int n0 = 2;
  double b[2] = {-1.0,1.0};
  double Fmin = 0.0;
  int maxit = 100;
  int trace = 0;
  int mask[2] = {1,1};
  double abstol = 0.0;
  double reltol = 1.0e-8;
  int nREPORT=1;
  int ex = 0;
  int fncount = 0;
  int grcount = 0;
  int fail = 0;
    
  vmmin(n0, b, &Fmin, &rosenbrock, &rosenbrock1,
	maxit, trace, mask,
	abstol, reltol, nREPORT, (void *) &ex,
	&fncount, &grcount, & fail);

  printf("Rosenbrock test:\nf = %17.15f\nx = %17.15f\ny = %17.15f\n", Fmin, b[0], b[1]);
  
  return 0;
  }
