#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSize.h"

#include <qwt_plot_curve.h>

#include "network.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef itk::Image<unsigned char,3> ImageType_u8;

#define V(a,b,c)  p_raw[(c)*xysize+(b)*width+(a)]
#define V_scaled(a,b,c)  p_scaled[(c)*xysize+(b)*width+(a)]

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void readRawData();
    void tiffInFileSelecter();
    void tiffOutFileSelecter();
    void thresholdChanged();
    void checkReady();
    void smoother();
    void sumIntensity();
    void sumScaledIntensity();
    void determineScale();
    void axisChanged();
    void saveScaledTiff();

private:
    Ui::MainWindow *ui;

public:
    FILE *fpout;
    int readTiff(const char *);
    int createTiff(const char *, unsigned char *buffer);
    void computeRawSum();
    void computeScaledSum();
    int summer(int w, double alpha, int *n);
    void scaleImage();

    QString tiffInFileName, tiffOutFileName;
    bool is_input_tiffname;
    bool is_output_tiffname;
    bool tiff_read;
    char axis;
    double threshold;
    bool ready;
    ImageType_u8::Pointer im_u8;
    long long width, height, depth, xysize;
    int len;
    double *rawSum;
    double *smoothSum;
    double *scaledSum;
    double *scale;
    int *nSum;
    bool use_average;
    bool scaled;
    unsigned char *p_raw;
    unsigned char *p_scaled;
    QwtPlot *qp;
    QwtPlotCurve *smoothCurve, *sumCurve, *scaledCurve;
};

#endif // MAINWINDOW_H
