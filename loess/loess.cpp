#include <algorithm>
#include <cmath>
#include <cstring>
#include "throw.h"
#include "loess.h"

using namespace std;
#define TRUE true
#define FALSE false
typedef bool Rboolean;
#define imin2(a,b) std::min(a,b)
#define imax2(a,b) std::max(a,b)
#define fmax2(a,b) std::max(a,b)

#define isnan(x) _isnan(x)	// GB

#define rPsort(x,n,k) std::partial_sort(&x[0],&x[n],&x[k]);

Loess::Loess():
//	smoother_span(2.0/3.0),	// GB
	smoother_span(1.0/3.0),
	nsteps(3),
	delta_speed(-1.0),
	paranoid(true)
    {
    }

Loess::~Loess()
    {
    }

double Loess::fsquare(double x)
    {
    return x * x;
    }

double Loess::fcube(double x)
    {
    return x * x * x;
    }

void Loess::set_smoother_span(double span)
{
	smoother_span = span;
}

void Loess::lowest(const double *x, const double *y, int n, const double *xs, double *ys,
	int nleft, int nright, double *w,
	Rboolean userw, double *rw, Rboolean *ok)
{
    int nrt, j;
    double a, b, c, h, h1, h9, r, range;

    x--;
    y--;
    w--;
    rw--;

    range = x[n]-x[1];
    h = fmax2(*xs-x[nleft], x[nright]-*xs);
    h9 = 0.999*h;
    h1 = 0.001*h;

    /* sum of weights */

    a = 0.;
    j = nleft;
    while (j <= n) {

	/* compute weights */
	/* (pick up all ties on right) */

	w[j] = 0.;
	r = fabs(x[j] - *xs);
	if (r <= h9) {
	    if (r <= h1)
		w[j] = 1.;
	    else
		w[j] = fcube(1.-fcube(r/h));
	    if (userw)
		w[j] *= rw[j];
	    a += w[j];
	}
	else if (x[j] > *xs)
	    break;
	j = j+1;
    }

    /* rightmost pt (may be greater */
    /* than nright because of ties) */

    nrt = j-1;
    if (a <= 0.)
	*ok = FALSE;
    else {
	*ok = TRUE;

	/* weighted least squares */
	/* make sum of w[j] == 1 */

	for(j=nleft ; j<=nrt ; j++)
	    w[j] /= a;
	if (h > 0.) {
	    a = 0.;

	    /*  use linear fit */
	    /* weighted center of x values */

	    for(j=nleft ; j<=nrt ; j++)
		a += w[j] * x[j];
	    b = *xs - a;
	    c = 0.;
	    for(j=nleft ; j<=nrt ; j++)
		c += w[j]*fsquare(x[j]-a);
	    if (sqrt(c) > 0.001*range) {
		b /= c;

		/* points are spread out */
		/* enough to compute slope */

		for(j=nleft; j <= nrt; j++)
		    w[j] *= (b*(x[j]-a) + 1.);
	    }
	}
	*ys = 0.;
	for(j=nleft; j <= nrt; j++)
	    *ys += w[j] * y[j];
    }
}


void Loess::clowess(
	const double *x, const double *y, int n,
	     double f, int nsteps, double delta,
	     double *ys, double *rw, double *res)
{
    int i, iter, j, last, m1, m2, nleft, nright, ns;
    Rboolean ok;
    double alpha, c1, c9, cmad, cut, d1, d2, denom, r, sc;

    if (n < 2) {
	ys[0] = y[0]; return;
    }

    /* nleft, nright, last, etc. must all be shifted to get rid of these: */
    x--;
    y--;
    ys--;


    /* at least two, at most n points */
    ns = imax2(2, imin2(n, (int)(f*n + 1e-7)));
#ifdef DEBUG_lowess
    REprintf("lowess(): ns = %d\n", ns);
#endif

    /* robustness iterations */

    iter = 1;
    while (iter <= nsteps+1) {
	nleft = 1;
	nright = ns;
	last = 0;	/* index of prev estimated point */
	i = 1;		/* index of current point */

	for(;;) {
	    if (nright < n) {

		/* move nleft,  nright to right */
		/* if radius decreases */

		d1 = x[i] - x[nleft];
		d2 = x[nright+1] - x[i];

		/* if d1 <= d2 with */
		/* x[nright+1] == x[nright], */
		/* lowest fixes */

		if (d1 > d2) {

		    /* radius will not */
		    /* decrease by */
		    /* move right */

		    nleft++;
		    nright++;
		    continue;
		}
	    }

	    /* fitted value at x[i] */

	    lowest(&x[1], &y[1], n, &x[i], &ys[i],
		   nleft, nright, res, iter>1, rw, &ok);
	    if (!ok) ys[i] = y[i];

	    /* all weights zero */
	    /* copy over value (all rw==0) */

	    if (last < i-1) {
		denom = x[i]-x[last];

		/* skipped points -- interpolate */
		/* non-zero - proof? */

		for(j = last+1; j < i; j++) {
		    alpha = (x[j]-x[last])/denom;
		    ys[j] = alpha*ys[i] + (1.-alpha)*ys[last];
		}
	    }

	    /* last point actually estimated */
	    last = i;

	    /* x coord of close points */
	    cut = x[last]+delta;
	    for (i = last+1; i <= n; i++) {
		if (x[i] > cut)
		    break;
		if (x[i] == x[last]) {
		    ys[i] = ys[last];
		    last = i;
		}
	    }
	    i = imax2(last+1, i-1);
	    if (last >= n)
		break;
	}
	/* residuals */
	for(i = 0; i < n; i++)
	    res[i] = y[i+1] - ys[i+1];

	/* overall scale estimate */
	sc = 0.;
	for(i = 0; i < n; i++) sc += fabs(res[i]);
	sc /= n;

	/* compute robustness weights */
	/* except last time */

	if (iter > nsteps)
	    break;
	/* Note: The following code, biweight_{6 MAD|Ri|}
	   is also used in stl(), loess and several other places.
	   --> should provide API here (MM) */
	for(i = 0 ; i < n ; i++)
	    rw[i] = fabs(res[i]);

	/* Compute   cmad := 6 * median(rw[], n)  ---- */
	/* FIXME: We need C API in R for Median ! */
	m1 = n/2;
	/* partial sort, for m1 & m2 */
	rPsort(rw, n, m1);
	if(n % 2 == 0) {
	    m2 = n-m1-1;
	    rPsort(rw, n, m2);
	    cmad = 3.*(rw[m1]+rw[m2]);
	}
	else { /* n odd */
	    cmad = 6.*rw[m1];
	}
#ifdef DEBUG_lowess
	REprintf("   cmad = %12g\n", cmad);
#endif
	if(cmad < 1e-7 * sc) /* effectively zero */
	    break;
	c9 = 0.999*cmad;
	c1 = 0.001*cmad;
	for(i = 0 ; i < n ; i++) {
	    r = fabs(res[i]);
	    if (r <= c1)
		rw[i] = 1.;
	    else if (r <= c9)
		rw[i] = fsquare(1.-fsquare(r/cmad));
	    else
		rw[i] = 0.;
	}
	iter++;
    }
}
//see also ..R-2.11.0/src/library/stats/R/lowess.R

std::auto_ptr<std::vector<double> > Loess::lowess(
	    const double *x,
	    const double *y,
	    int32_t n
	    )
    {
    if(n<1) THROW("n =" << n << "<1");
	printf("lowess: %d\n",n);
    if(paranoid)
	{
	for(int32_t i=0;i< n;++i)
	    {
	    if(isnan(x[i])) THROW("NAN: x["<<i<<"]");
	    if(isnan(y[i])) THROW("NAN: y["<<i<<"]");
	    if(i>0)
		{
		if(x[i-1]> x[i]) THROW("Data not sorted on x");
		}
	    }
	}
    double delta=this->delta_speed;
    if(delta<0.0)
	{
	delta=.01*(x[n-1]-x[0]);
	}
    double* rw=new double[n];
    double* ys=new double[n];
    double* res=new double[n];
    memset((void*)rw,0,sizeof(double)*n);
    memset((void*)ys,0,sizeof(double)*n);
    memset((void*)res,0,sizeof(double)*n);
    clowess(x, y, n, this->smoother_span, this->nsteps, delta, ys, rw, res);


    std::vector<double>* v=new std::vector<double>(n,0.0);
    std::copy(&ys[0],&ys[n],v->begin());
    delete[] ys;
    delete[] rw;
    delete[] res;
    return std::auto_ptr<std::vector<double> >(v);
    }


__declspec(dllexport) int loess_interface(int n, double *xdata, double *ydata, double span, double *xsmooth, double *ysmooth)
{
    Loess app;
    vector<double> xs;
    vector<double> ys;

	printf("loess: span: %lf\n",span);

	app.set_smoother_span(span);
	for (int i=0; i<n; i++) {
		double x,y;
		x = xdata[i];
		y = ydata[i];
		xs.push_back(x);
		ys.push_back(y);
	}
    auto_ptr<vector<double> > y2= app.lowess(&xs.front(),&ys.front(),xs.size());
	if (xs.size() != n) {
		printf("loess: xs.size() != n: %d %d\n",xs.size(),n);
		return 1;
	}
    for (int i=0;i< (int)y2->size();++i) {
		xsmooth[i] = xs[i];
		ysmooth[i] = y2->at(i);
	}
    return 0;
}
