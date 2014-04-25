 //Qt
#include "mainwindow.h"
#include "plot.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QFileInfoList>
#include <QTextStream>
#include <QDir>
#include <iostream>
//ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSize.h"

using namespace std;
int loess_interface(int n, double *xdata, double *ydata, double span, double *xsmooth, double *ysmooth);

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->textEdit->setReadOnly(true);
    QString infoFile = QCoreApplication::applicationDirPath() + "/info/equaliser_info.txt";
    QFile file(infoFile);
    bool ok = file.open(QIODevice::ReadOnly | QIODevice::Text);
    if (!ok) {
        ui->textEdit->append("The information file is missing:");
        ui->textEdit->append(infoFile);
    } else {
        QTextStream in(&file);
        QString line = in.readLine();
        while (!line.isNull()) {
            ui->textEdit->append(line);
            line = in.readLine();
        }
        file.close();
        ui->textEdit->moveCursor(QTextCursor::Start);
    }

    fpout = fopen("equaliser.out","w");
    is_input_tiffname = false;
    is_output_tiffname = false;
    tiff_read = false;
    scaled = false;
    rawSum = NULL;
    scaledSum = NULL;
    smoothSum = NULL;
    scale = NULL;
    axis = 'Z';
    use_average = false;
    ui->lineEdit_threshold->setText("10");
    ui->lineEdit_smoothing->setText("0.05");
    ui->pushButtonRead->setDisabled(true);
    ui->pushButtonSum->setDisabled(true);
    ui->pushButtonSmooth->setDisabled(true);
    ui->pushButtonScale->setDisabled(true);
    ui->pushButtonSave->setDisabled(true);    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_frame");
    qp->setTitle("Summed thresholded intensity");
    sumCurve = 0;
    smoothCurve = 0;
}

MainWindow::~MainWindow()
{
//    if (p_raw) free(p_raw);
//    if (rawSum) free(rawSum);
//    if (smoothSum) free(smoothSum);
//    if (scaledSum) free(scaledSum);
//    if (nSum) free(nSum);
//    if (scale) free(scale);
    delete ui;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::tiffInFileSelecter()
{
    ui->labelResult->setText("");
    tiffInFileName = QFileDialog::getOpenFileName(this,
        tr("Input TIFF file"), ".", tr("TIFF Files (*.tif)"));
    if (tiffInFileName.compare(ui->labelInputTiffFile->text()) != 0) {
        printf("%s\n",tiffInFileName.toAscii().constData());
        printf("%s\n",ui->labelInputTiffFile->text().toAscii().constData());
    }
    ui->labelInputTiffFile->setText(tiffInFileName);
    is_input_tiffname = true;
    tiff_read = false;
    ui->pushButtonRead->setEnabled(true);
    ui->pushButtonSum->setDisabled(true);
    ui->pushButtonSmooth->setDisabled(true);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::tiffOutFileSelecter()
{
    ui->labelResult->setText("");
    tiffOutFileName = QFileDialog::getSaveFileName(this,
        tr("Output TIFF file"), ".", tr("TIFF Files (*.tif)"));
    if (tiffOutFileName.compare(ui->labelOutputTiffFile->text()) != 0) {
        printf("%s\n",tiffOutFileName.toAscii().constData());
        printf("%s\n",ui->labelOutputTiffFile->text().toAscii().constData());
    }
    ui->labelOutputTiffFile->setText(tiffOutFileName);
    is_output_tiffname = true;
    if (scaled) {
        ui->pushButtonSave->setEnabled(true);
    }
//    checkReady();
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::axisChanged()
{
    if (ui->radioButton_X->isChecked()) {
        axis = 'X';
    } else if (ui->radioButton_Y->isChecked()) {
        axis = 'Y';
    } else if (ui->radioButton_Z->isChecked()) {
        axis = 'Z';
    }
    ui->pushButtonSum->setEnabled(true);
    ui->pushButtonSmooth->setDisabled(true);
}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::readRawData()
{
    int res = readTiff(tiffInFileName.toAscii().constData());
    tiff_read = (res == 0);
    if (tiff_read) {
        ui->pushButtonRead->setDisabled(true);
    }
    ui->pushButtonSum->setEnabled(true);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::sumIntensity()
{
    double *x, *y;
    int i;

    if (smoothCurve != 0) {
        smoothCurve->detach();
    }
    if (sumCurve != 0) {
        sumCurve->detach();
    }
    if (scaledCurve != 0) {
        scaledCurve->detach();
    }
    qp->replot();
    computeRawSum();
    x = (double *)malloc(len*sizeof(double));
    y = (double *)malloc(len*sizeof(double));
    for (i=0; i<len; i++) {
        x[i] = i;
        y[i] = rawSum[i];
    }
    if (sumCurve == 0) {
        sumCurve = new QwtPlotCurve("raw");
    }
    sumCurve->setData(x, y, len);
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    QPen *pen = new QPen();
    pen->setColor(pencolor[0]);
    sumCurve->setPen(*pen);
    sumCurve->attach(qp);
    qp->replot();
    free(x);
    free(y);
    ui->pushButtonSum->setDisabled(true);
    ui->pushButtonSmooth->setEnabled(true);
    ui->pushButtonScale->setDisabled(true);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::smoother()
{
    double *x, *y, *xs, *ys;
    double span;
    int i;

    if (scaledCurve != 0) {
        scaledCurve->detach();
    }
    qp->replot();
    x = (double *)malloc(len*sizeof(double));
    y = (double *)malloc(len*sizeof(double));
    xs = (double *)malloc(len*sizeof(double));
    ys = (double *)malloc(len*sizeof(double));
    for (i=0; i<len; i++) {
        x[i] = i;
        y[i] = rawSum[i];
    }
    span = ui->lineEdit_smoothing->text().toDouble();
    loess_interface(len,x,y,span,xs,ys);

    for (i=0; i<len; i++) {
        smoothSum[i] = ys[i];
        fprintf(fpout,"%4d %8.2lf %8.2lf\n",i,xs[i],ys[i]);
    }
    if (smoothCurve == 0) {
        smoothCurve = new QwtPlotCurve("smooth");
    }
    smoothCurve->setData(xs, ys, len);
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    QPen *pen = new QPen();
    pen->setColor(pencolor[1]);
    smoothCurve->setPen(*pen);
    smoothCurve->attach(qp);
    qp->replot();
    free(x);
    free(y);
    free(xs);
    free(ys);
    ui->pushButtonScale->setEnabled(true);
    ui->pushButtonSave->setDisabled(true);    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    scaled = false;
}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::computeRawSum()
{
    int x, y, z, n;

    printf("computeRawSum\n");
    threshold = ui->lineEdit_threshold->text().toDouble();
    if (rawSum) {
        free(rawSum);
        free(nSum);
    }
    if (smoothSum) {
        free(smoothSum);
    }
    if (axis == 'X') {
        len = width;
    } else if (axis == 'Y') {
        len = height;
    } else {
        len = depth;
    }
    rawSum = (double *)malloc(len*sizeof(double));
    nSum = (int *)malloc(len*sizeof(int));
    smoothSum = (double *)malloc(len*sizeof(double));
    if (axis == 'X') {
        for (x=0; x<width; x++) {
            rawSum[x] = summer(x,1.0,&n);
            nSum[x] = n;
        }
    } else if (axis == 'Y') {
        for (y=0; y<height; y++) {
            rawSum[y] = summer(y,1.0,&n);
            nSum[y] = n;
        }
    } else if (axis == 'Z') {
        for (z=0; z<depth; z++) {
            rawSum[z] = summer(z,1.0,&n);
            nSum[z] = n;
            fprintf(fpout,"%4d %6d %10.0f\n",z,n,rawSum[z]);
        }
    }
    fflush(fpout);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::computeScaledSum()
{
    int x, y, z, n;
    unsigned char v;
    int sum;

    printf("computeScaledSum\n");
    fflush(stdout);
    threshold = ui->lineEdit_threshold->text().toDouble();
    if (scaledSum) {
        free(scaledSum);
//        free(nSum);
    }
    scaledSum = (double *)malloc(depth*sizeof(double));
//    nSum = (int *)malloc(depth*sizeof(int));
    for (z=0; z<depth; z++) {
        sum = 0;
        n = 0;
        for (x=0; x<width; x++) {
            for (y=0; y<height; y++) {
                v = V_scaled(x,y,z);
                if (v > threshold) {
                    sum += v;
                    n++;
                }
            }
        }
        scaledSum[z] = sum;
//        nSum[z] = n;
    }
}

//----------------------------------------------------------------------------------------
// For a scaling value alpha, compute the sum of scaled intensities that exceed the threshold.
// Or, sum of intensities that exceed the threshold, scaled.
// The problem is that either way the scaled image is not right.
// If the scaled intensities are compared with the threshold, after scaling more pixels
// are included in the sum, and therefore the sum is increased by more than the scale factor,
// but the average could go either way.
// If the unscaled intensities are compared with the threshold, then those above scaled,
// after scaling many pixels will be removed from the sum, and the resulting scaled image
// will not have the expected sum.
//----------------------------------------------------------------------------------------
int MainWindow::summer(int w, double alpha, int *n)
{
    int x, y, z, k;
    int sum = 0;
    unsigned char v;
    int val;

    k = 0;
    if (axis == 'X') {
        for (y=0; y<height; y++) {
            for (z=0; z<depth; z++) {
                v = V(w,y,z);
                val = int(alpha*v + 0.5);
                if (val > threshold) {
                    sum += val;
                    k++;
                }
            }
        }
    } else if (axis == 'Y') {
        for (x=0; x<width; x++) {
            for (z=0; z<depth; z++) {
                v = V(x,w,z);
                val = int(alpha*v + 0.5);
                if (val > threshold) {
                    sum += val;
                    k++;
                }
            }
        }
    } else if (axis == 'Z') {
        for (x=0; x<width; x++) {
            for (y=0; y<height; y++) {
                v = V(x,y,w);
                val = int(alpha*v + 0.5);
                if (val > threshold) {
//                  if (v > threshold) {
//                    val = int(alpha*v + 0.5);
                    sum += val;
                    k++;
                }
            }
        }
    }
    *n = k;
    return sum;
}

//----------------------------------------------------------------------------------------
// This applies only when axis = 'Z'
// The aim is to determine the scaling values scale[z] such that:
//   Sum(a*V(x,y,z) < T) = L = smoothSum[z]
// Define E(a) = Sum(a*V(x,y,z)) - smoothSum[z]
// and find a s.t. E(a) = 0 by Euler's method:
// a(n+1) = a(n) - E(a(n))/E'(a(n))
//----------------------------------------------------------------------------------------
void MainWindow::determineScale()
{
    int z, n;
    double a0, a1, a2, a, E0, E1, E2, anew, da;

    printf("determineScale\n");
    fflush(stdout);
    scale = (double *)malloc(depth*sizeof(double));
    for (z=0; z<depth; z++) {
        if (rawSum[z] == 0) {
            scale[z] = 1.0;
            continue;
        }
        E0 = rawSum[z] - smoothSum[z];
        a0 = smoothSum[z]/rawSum[z];
        a1 = a0;
//        printf("%4d %8.0f %8.0f %8.4f\n",z,rawSum[z],smoothSum[z],a1);
//        fflush(stdout);
        E1 = (summer(z,a1,&n) - smoothSum[z]);
        if (E1 > 0)
            da = -0.01;
        else
            da = 0.01;
        for (;;) {
            a2 = a1 + da;
            E2 = (summer(z,a2,&n) - smoothSum[z]);
            if (E2*da > 0) {    // crossed
                a = a1 - E1*(a2-a1)/(E2-E1);
                break;
            }
            a1 = a2;
            E1 = E2;
        }
        scale[z] = a;
        printf("z,scale: %4d %8.0f %8.0f %8.4f\n",z,rawSum[z],smoothSum[z],scale[z]);
        fflush(stdout);
        fprintf(fpout,"z,raw,smoothed,scale: %4d %8.0f %8.0f %8.4f %8.4f\n",z,rawSum[z],smoothSum[z],a0,scale[z]);
    }
//    fclose(fpout);
    ui->pushButtonScale->setDisabled(true);
    scaleImage();     //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    printf("did scaleImage\n");
    fflush(stdout);
    sumScaledIntensity(); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::scaleImage()
{
    int x, y, z, val;
    unsigned char v;

    printf("scaleImage\n");
    fflush(stdout);
//    p_scaled = (unsigned char *)malloc(width*height*depth*sizeof(unsigned char));
    p_scaled = p_raw;
    for (x=0; x<width; x++) {
        for (y=0; y<height; y++) {
            for (z=0; z<depth; z++) {
                v = V(x,y,z);
                val = int(scale[z]*v + 0.5);
                val = MIN(val,255);
                if (val > threshold) {
//                if (v > threshold) {
//                    val = int(scale[z]*v + 0.5);
                    V_scaled(x,y,z) = (unsigned char)val;
                } else {
                    V_scaled(x,y,z) = v;    // TRY THIS
                }
            }
        }
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::sumScaledIntensity()
{
    double *x, *y;
    int i;

    printf("sumScaledIntensity\n");
    fflush(stdout);
//    if (smoothCurve != 0) {
//        smoothCurve->detach();
//    }
    if (scaledCurve != 0) {
        scaledCurve->detach();
    }
    qp->replot();
    computeScaledSum();
    printf("did computeScaledSum\n");
    fflush(stdout);
    x = (double *)malloc(len*sizeof(double));
    y = (double *)malloc(len*sizeof(double));
    for (i=0; i<len; i++) {
        x[i] = i;
        y[i] = scaledSum[i];
    }
    if (scaledCurve == 0) {
        scaledCurve = new QwtPlotCurve("scaled");
    }
    scaledCurve->setData(x, y, len);
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    QPen *pen = new QPen();
    pen->setColor(pencolor[3]);
    scaledCurve->setPen(*pen);
    scaledCurve->attach(qp);
    qp->replot();
    free(x);
    free(y);
    scaled = true;
    if (is_output_tiffname) {
        ui->pushButtonSave->setEnabled(true);
    }
//    ui->pushButtonSum->setDisabled(true);
//    ui->pushButtonSmooth->setEnabled(true);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::saveScaledTiff()
{
    createTiff(tiffOutFileName.toAscii().constData(),p_scaled);
//    createTiff(tiffOutFileName.toAscii().constData(),p_raw);    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::thresholdChanged()
{
    ui->pushButtonSum->setEnabled(true);
    ui->pushButtonSmooth->setDisabled(true);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::checkReady()
{
//    bool voxelOK;
//    voxelsize[0] = ui->lineEdit_xvoxel->text().toFloat();
//    voxelsize[1] = ui->lineEdit_yvoxel->text().toFloat();
//    voxelsize[2] = ui->lineEdit_zvoxel->text().toFloat();
//    margin = ui->lineEditMargin->text().toFloat();

//    if (voxelsize[0] > 0 && voxelsize[1] > 0 && voxelsize[2] > 0)
//        voxelOK = true;
//    else
//        voxelOK = false;
//    if (voxelOK && margin > 0 && is_am && is_tiff) {
//        ready = true;
//        ui->pushButtonGo->setEnabled(true);
//    } else {
//        ready = false;
//        ui->pushButtonGo->setEnabled(false);
//    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//    ui->labelResult->setText("");
//    if (!ready) {
//        ui->labelResult->setText("Not ready: select files");
//        return;
//    }
//    if (!am_read) {
//        resultstr = "reading network file...";
//        ui->labelResult->setText(resultstr);
//        QCoreApplication::processEvents();
//        err = readNetwork(&network, amFileName.toAscii().constData());
//        if (err != 0) {
//            resultstr = "FAILED: read_network";
//            ui->labelResult->setText(resultstr);
//            return;
//        }
//    }
//    resultstr = "creating image file...";
//    ui->labelResult->setText(resultstr);
//    QCoreApplication::processEvents();
//    err = createTiffData(&network);
//    if (err != 0) {
//        resultstr = "FAILED: createTiffData";
//        ui->labelResult->setText(resultstr);
//        return;
//    }
//    err = createTiff(tiffFileName.toAscii().constData(),p_im,width,height,depth);
//    if (err != 0) {
//        resultstr = "FAILED: createTiff";
//        ui->labelResult->setText(resultstr);
//        return;
//    }
//    if (p_im) free(p_im);
//    checkReady();
//    resultstr = "SUCCESS";
//    ui->labelResult->setText(resultstr);
//    return;
//}


