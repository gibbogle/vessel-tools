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
    void am_tiffer();
    void amFileSelecter();
    void tiffFileSelecter();
    void dataChanged();
    void checkReady();
private:
    Ui::MainWindow *ui;

public:
    FILE *fpout;
    int readAmiraFile(const char *amFile, NETWORK *net);
    int writeAmiraFile(const char *amFileOut, const char *amFileIn, NETWORK *net);
    int readTiff(const char *, int *, int *, int *);
    int createTiffData();
    int createTiff(const char *, unsigned char *, int, int, int);
    int readNetwork(NETWORK *, const char *);
    int createNetwork(NETWORK *, NETWORK *, DEADEND *, int);
    void freeNetwork(NETWORK *);
    bool inSphere(APOINT p);

    QString amFileName;
    QString tiffFileName;
    bool is_am;
    bool is_tiff;
    bool am_read;
    bool tiff_read;
    bool ready;
    bool is_network;
    NETWORK network;
    DEADEND *deadlist;
    ImageType_u8::Pointer im_u8;
    int width, height, depth, xysize;
    float margin;
    unsigned char *p_im;
    float voxelsize[3];
    float sphereCentre[3], sphereRadius;
};

#endif // MAINWINDOW_H
