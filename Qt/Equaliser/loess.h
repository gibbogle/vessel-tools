
/**
 Implementation of the loess algorithm in C++
 original code R project: http://www.r-project.org/
 http://svn.r-project.org/R/trunk/src/library/stats/src/lowess.c

 TEST:

Using #RStats

    > T<-read.table('data.txt')
    > T
      V1   V2
    1  1  1.5
    2  3  2.0
    3  6  7.0
    4  9  8.0
    5 12 15.0

    > lowess(T$V1,T$V2)
    $x
    [1]  1  3  6  9 12

    $y
    [1] 1.5 2.0 7.0 8.0 8.0

The same dataset, but using this class


    $ ./a.out < data.txt
    1	1.5	1.5
    3	2	2
    6	7	7
    9	8	8
    12	15	8

 */
#ifndef LOESS_STATS_H
#define LOESS_STATS_H
#include <stdint.h>
#include <vector>
#include <memory>

//#define STANDALONE_VERSION 1

class Loess
    {
    public:
	Loess();
	~Loess();
	/** proportion of points in the plot which influence the smooth at each value */
	double smoother_span;
	/** the number of ‘robustifying’ iterations which should be performed. */
	int32_t nsteps;
	/** used to speed up computation */
	double delta_speed;
	/** perform some basic checks */
	bool paranoid;

	std::auto_ptr<std::vector<double> > lowess(const double *x, const double *y, int32_t n);
	void set_smoother_span(double span);

	private:
	static double fsquare(double x);
	static double fcube(double x);

	void lowest(const double *x, const double *y, int n, const double *xs, double *ys,
		int nleft, int nright, double *w,
		bool userw, double *rw, bool *ok);

	void clowess(const double  *x, const double *y, int n,
		     double f, int nsteps, double delta,
		     double *ys, double *rw, double *res);

    };

__declspec(dllexport) int loess_interface(int n, double *xdata, double *ydata, double span, double *xsmooth, double *ysmooth);

#endif
