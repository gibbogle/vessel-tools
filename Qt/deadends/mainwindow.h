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
    void amFileSelecter();
    void tiffFileSelecter();
    void outputFileSelecter();
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
    int createTiff(char *, unsigned char *, int, int, int);
    int readNetwork(NETWORK *, const char *);
    int createNetwork(NETWORK *, NETWORK *, DEADEND *, int);
    void freeNetwork(NETWORK *);
    int findDeadends(NETWORK *, DEADEND **, int *);
    int getIntensities(NETWORK *, DEADEND *, int);
    void getVoxels(float r, float c[], float beta, int *nv, int v[][3]);
    void getMaxIntensity(int nv, int v[][3], int *maxval);
    void getSphere();
    bool inSphere(APOINT p);
    QString amFileName;
    QString tiffFileName;
    QString outputFileName;
    bool is_am;
    bool is_tiff;
    bool is_output;
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
    unsigned char *p_im;
    float voxelsize[3];
    bool isSphere;
    float sphereRadius;
    float sphereCentre[3];
    float maxRadius;
    bool detected;
    bool evaluated;
};

#endif // MAINWINDOW_H
