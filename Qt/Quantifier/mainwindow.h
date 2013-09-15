#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
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
    void computeVessels();
    int checkSlice(int, int);
    void voxelChanged();
    void checkReady();
    void doSetup();
    void sliceChanged();
    void averageChanged();
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
    void on_radioButton_centre_toggled(bool checked);

private:
    Ui::MainWindow *ui;

public:
    int getArea(int axis, int islice, float *area);
    int getVolume(float *volume, int *ntvoxels);
    int histology(int axis, int islice, int *np, float *area);
    int average_histology(int *np, float *area);
    void fillEllipse(float z0, float S1[], float S2[], float diam, float vsize[], int nv[], int *npixels);
    int setup(char *input_amfile, char *close_file, char *result_file, float vsize[]);
    int getCloseSize(int nvoxels[]);
    bool isSetup();
    void reset();
    int ReadAmiraFile(char *, NETWORK *);
    int ReadCloseFile(char *);
    int branching(int *nbranchpts, float *totlen, float *totvol);
    float rangeVolume();

    FILE *fpout;
    float voxelsize[3];
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
    int range_x1;
    int range_x2;
    int range_y1;
    int range_y2;
    int range_z1;
    int range_z2;
    int range[3][2];

    ImageViewer *imageViewer;
};

#endif // MAINWINDOW_H
