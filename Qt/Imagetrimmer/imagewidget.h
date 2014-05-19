/*
 * File:   imagewidget.h
 * Author: zian
 *
 * Created on 19 de septiembre de 2011, 19:34
 */

#ifndef IMAGEWIDGET_H
#define	IMAGEWIDGET_H

#include <QtGui>
#include <QWidget>
#include <QVTKWidget.h>
#include <QFile>
#include <QCoreApplication>
#include <QTextStream>

#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <itkImage.h>
#include <itkRGBPixel.h>

#include "external/itkVTKImageToImageFilter.h"
#include "external/itkImageToVTKImageFilter.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


typedef itk::RGBPixel< unsigned char > RGBPixelType;

typedef itk::Image< unsigned char, 2 > ImageType2;
typedef itk::Image< unsigned char, 3 > ImageType3;
typedef itk::Image< float, 2 > FloatImageType;
typedef itk::Image< RGBPixelType, 2 > RGBImageType;

typedef itk::VTKImageToImageFilter <ImageType2> itkConnectorType;
typedef itk::ImageToVTKImageFilter <ImageType2> vtkConnectorType;
typedef itk::ImageToVTKImageFilter <RGBImageType> RGBVtkConnectorType;
typedef itk::ImageToVTKImageFilter <FloatImageType> vtkFloatConnectorType;

class ImageWidget : public QWidget {
    Q_OBJECT

public:
    /** 
     * Constructor for this ImageWidget 
     */
//    ImageWidget(QWidget *parent = 0);
    ImageWidget(QWidget * parent, QWidget *page);

    /** 
     * Destructor for this ImageWidget 
     */
    virtual ~ImageWidget();

    /**
     * load an display an image from file
     */
    void openWithVTK();

    /**
     * load an display an image from file
     */
    void openWithITK();

    bool frameOK(int);
    void showFrame(int);

    QString getShortFileName();
    QString getFileName();
    int getDepth();
    void getInfo(int *, int *, int *, int *);
    unsigned char *getBuffer();
    void clearAll();
    void clearActors();
    void clearImages();
    void addMark(int *);
    void moveMark(int, int *);
    bool selectMark(int pos[3], int *iactor);
    void deleteMark(int);
    void setNumber(int iactor);
    void saveLine();
    void showActors();
    void showLines();
    vtkSmartPointer<vtkActor2D> makeMarkActor(double *pos);
    void convertPosToPixel(double *pos, int *pixel);
    void convertPixelToPos(int *pixel,double *pos);
    void writeLineFile(QString);
    void readLineFile(QString);

public slots:
    void tallyMark(QMouseEvent* event);

signals:
    void SaveStatus(bool);

private:

    QWidget *m_parent;
    QVTKWidget *qvtkWidget;

    /** The image displayed for this window */
    ImageType2::Pointer itkImage2;
    ImageType3::Pointer itkImage3;
    RGBImageType::Pointer rgbItkImage;

    /** The VTK image to display in this window */
    vtkSmartPointer <vtkImageData> vtkImage;

    vtkSmartPointer<vtkImageActor> imageActor;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor; 
    vtkSmartPointer<vtkCamera> camera;
    vtkSmartPointer<vtkPolyDataMapper2D> mapper;

    /** The type of the image components RGB, scalar, etc */
    std::string pixelType;

    /** The type of the image pixels */
    std::string imageType;

    /** The number of the image dimensions */
    size_t numDimensions;
    
    /** flag for a flipped image */
    bool isFlipped;

    /** flag for kind of fit to the render window */
    bool horizontalFit;

    /**
     * to display the loaded image in the QVTKwidget
     * @param image vtkImageData
     */
    void displayImage(vtkImageData *image);

    /**
     * Set itkImage converting the vtkImage to a ITK image
     */
    void setITKImageFromVTK();

    /**
     * extract some image properties as, pixel type, image type and number of dimensions
     * @param fileName path to the file 
     * @param verbose if print the standar out
     */
    void setImageProperties(std::string fileName, bool verbose);

    void getFrame(int z);

    void show2DImage();

    void getImageXY(int pos[], int xy[]);

    long long imageWidth;
    long long imageHeight;
    long long imageDepth;
    int iframe;
    QString imageFileName;

    int rect[2];
    bool LButtonDown;
    std::vector<vtkSmartPointer<vtkActor2D>> Actors;
    std::vector<vtkSmartPointer<vtkTextActor>> NumActors;
    bool selected;
    int iselect;
    bool lineSaved;
    typedef std::vector<vtkSmartPointer<vtkActor2D>> vect_actor_t;
    std::vector<vect_actor_t> vect_actor_mat;
};

#endif	/* IMAGEWIDGET_H */

