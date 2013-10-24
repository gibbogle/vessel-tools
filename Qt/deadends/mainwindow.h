#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSize.h"

#include "network.h"

typedef itk::Image<unsigned char,3> ImageType_u8;

#define V(a,b,c)  p_im[(c)*xysize+(b)*width+(a)]
#define V_dead(a,b,c)  p_im_dead[(c)*xysize_d+(b)*width_d+(a)]

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
    void inputAmiraFileSelecter();
    void inputTiffFileSelecter();
    void outputAmiraFileSelecter();
    void outputTiffFileSelecter();
    void voxelChanged();
    void checkReady();
    void detectDeadends();
    void evaluate();
    void getAllIntensities();
    void optionSphere();
    void saveDeadendFile();
private:
    Ui::MainWindow *ui;

public:
    FILE *fpout;
    int readAmiraFile(const char *amFile, NETWORK *net);
    int writeAmiraFile(const char *amFileOut, const char *amFileIn, NETWORK *net);
    int readTiff(const char *, int *, int *, int *);
    int createTiff(const char *, unsigned char *, int, int, int);
    int createTiffData(NETWORK *net);
    int readNetwork(NETWORK *, const char *);
    int createNetwork(NETWORK *, NETWORK *, DEADEND *, int);
    void freeNetwork(NETWORK *);
    int findDeadends(NETWORK *, DEADEND **, int *);
    int getIntensities(NETWORK *, DEADEND *, int, int *nzero);
    void getVoxels(float r, float c[], float beta, int *nv, int v[][3]);
    void getMaxIntensity(int nv, int v[][3], int *maxval);
    void getSphere();
    bool inSphere(APOINT p);
    void showEdge(int ie, NETWORK *net);
    QString inputAmiraFileName;
    QString inputTiffFileName;
    QString outputAmiraFileName;
    QString outputTiffFileName;
    bool is_am_in;
    bool is_tiff_in;
    bool is_am_out;
    bool is_tiff_out;
    bool am_read;
    bool tiff_read;
    bool ready;
    bool is_network;
    NETWORK network;
    NETWORK deadnetwork;
    int ndead;
    DEADEND *deadlist;
    ImageType_u8::Pointer im_u8;
    int width, height, depth, xysize;
    int width_d, height_d, depth_d, xysize_d;
    int margin;
    unsigned char *p_im;
    unsigned char *p_im_dead;
    float voxelsize[3];
    bool isSphere;
    float sphereRadius;
    float sphereCentre[3];
    float maxRadius;
    bool detected;
    bool evaluated;
};

#endif // MAINWINDOW_H
