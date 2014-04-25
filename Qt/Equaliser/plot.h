#ifndef PLOT_H
#define PLOT_H

#include <QtGui>
#include <QPen>
#include <QFileDialog>
#include <qwt_plot.h>
#include <qwt_painter.h>
#include <qwt_plot_canvas.h>
#include <qwt_plot_marker.h>
#include <qwt_plot_curve.h>
#include <qwt_scale_widget.h>
#include <qwt_legend.h>
#include <qwt_scale_draw.h>
#include <qwt_math.h>

using namespace std;

class Plot : public QwtPlot
{
//	Q_OBJECT

public:
    Plot(QString, QString, QWidget *parent = 0);
    ~Plot();

	void mousePressEvent (QMouseEvent *);
	void addCurve(QString);
	void removeCurve(QString);
	void removeAllCurves();
    void redraw(double *, double *, int, QString, QString);
//    void redraw2(double *, double *, double *, int);
	void redraw2(double *, double *, double *, double *, int, int);
	void draw2(double *, double *, double *, double *, int, int);
	void setYScale(double);
	double calc_yscale(double);

	QString name;
	static const int ncmax = 8;
	QwtPlotCurve *curve[ncmax];
	double xscale;
	double yscale;
	int ncurves;
	char msg[1024];
//signals:
//	void rButtonClicked(QString text, int w, int h);

};

#endif
