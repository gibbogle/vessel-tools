#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#ifdef HAVE_QT5
#include <QtWidgets/QMainWindow>
#else
#include <QtGui/QMainWindow>
#endif
#include "network.h"
#include "imageviewer.h"

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void amFileSelecter();
    void closeFileSelecter();
    void resultFileSelecter();
    void computeVolume();
    void computeArea();
    void ComputeVessels();
    int checkSlice(int, int);
    void voxelChanged();
    void checkReady();
    void doSetup();
//    void sliceChanged();
//    void averageChanged();
    void checkBox_x();
    void checkBox_y();
    void checkBox_z();
    void getRanges();
    void getCentredRanges();
    void getAveragingRanges();
    void checkRanges();
    void reader();
//    void viewImage();

private slots:
    void on_radioButton_slice_toggled(bool checked);
    void on_radioButton_centre_toggled(bool checked);
    void on_checkBox_selectblock_toggled(bool checked);
    void on_radioButton_xaxis_toggled(bool checked);
    void on_radioButton_yaxis_toggled(bool checked);
    void on_radioButton_zaxis_toggled(bool checked);

private:
    Ui::MainWindow *ui;

public:
    int getArea(int axis, int islice, int *npixels, double *area);
    int TotalVoxelCount();
    int SliceHistology(int axis, int islice, int *nvessels, int *nvesselpixels, int *nslicepixels, double *slicearea);
    int VolumeHistology(bool *use_axis, int *nvessels, int *nvesselpixels, int *ntissuepixels, double *tissuearea);
    void fillEllipse(double z0, double S1[], double S2[], double diam, double vsize[], int nv[], int rng_x[], int rng_y[], int *npixels);
    int setup(char *input_amfile, char *close_file, char *result_file, double vsize[]);
    void setBlockAxes();
    void enableRange(char, bool);
    int getCloseSize(int nvoxels[]);
    bool isSetup();
    void reset();
    int ReadAmiraFile(char *, NETWORK *);
    int ReadCloseFile(char *);
    int branching(int *nbranchpts, double *totlen, double *totvol);
    void VesselDensity(double dmin, double dmax, double *vessellength_mm, int *nbranchpts, double *tissuevolume_mm3);
    int RangeVoxelCount();

    FILE *fpout;
    char *axisstr[3];
    double voxelsize[3];
    int nvoxels[3];
    QString amFileName;
    QString closeFileName;
    QString resultFileName;
    bool is_amfile;
    bool is_closefile;
    bool is_resultfile;
    bool is_ready;
    bool is_slice;
    bool is_average;
    bool is_block;
    int range[3][2];

    ImageViewer *imageViewer;
};

#define DEBUG false

#endif // MAINWINDOW_H
