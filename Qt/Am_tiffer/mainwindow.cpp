 //Qt
#include "mainwindow.h"
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

//#include "itkImageSeriesReader.h"
//#include "itkImageToVTKImageFilter.h"
/*
//VTK
#include "vtkImageViewer2.h"
#include "vtkResliceImageViewer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
*/
using namespace std;

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->textEdit->setReadOnly(true);
    QString infoFile = QCoreApplication::applicationDirPath() + "/info/am_tiffer_info.txt";
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

    fpout = fopen("am_tiffer.out","w");
    network.ne = 0;
    network.nv = 0;
    network.np = 0;
    is_am = false;
    is_tiff = false;
    am_read = false;
    voxelsize[0] = voxelsize[1] = voxelsize[2] = 2;
    checkReady();
}

MainWindow::~MainWindow()
{
    delete ui;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::amFileSelecter()
{
    ui->labelResult->setText("");
    amFileName = QFileDialog::getOpenFileName(this,
        tr("Input Amira file"), ".", tr("Amira Files (*.am)"));
    if (amFileName.compare(ui->labelAmFile->text()) != 0) {
        am_read = false;
    }
    ui->labelAmFile->setText(amFileName);
    is_am = true;
    checkReady();
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::tiffFileSelecter()
{
    ui->labelResult->setText("");
    tiffFileName = QFileDialog::getSaveFileName(this,
        tr("Output TIFF file"), ".", tr("TIFF Files (*.tif)"));
    if (tiffFileName.compare(ui->labelTiffFile->text()) != 0) {
        printf("%s\n",tiffFileName.toAscii().constData());
        printf("%s\n",ui->labelTiffFile->text().toAscii().constData());
    }
    ui->labelTiffFile->setText(tiffFileName);
    is_tiff = true;
    checkReady();
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::dataChanged()
{
    checkReady();
}

/*
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::optionSphere()
{
    if (ui->checkBoxSphere->isChecked()){
        isSphere = true;
        ui->lineEditRadius->setEnabled(true);
        ui->lineEditX->setEnabled(true);
        ui->lineEditY->setEnabled(true);
        ui->lineEditZ->setEnabled(true);
    } else {
        isSphere = false;
        ui->lineEditRadius->setEnabled(false);
        ui->lineEditX->setEnabled(false);
        ui->lineEditY->setEnabled(false);
        ui->lineEditZ->setEnabled(false);
    }
    detected = false;
    checkReady();
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::getSphere()
{
    sphereCentre[0] = ui->lineEditX->text().toFloat();
    sphereCentre[1] = ui->lineEditY->text().toFloat();
    sphereCentre[2] = ui->lineEditZ->text().toFloat();
    sphereRadius = ui->lineEditRadius->text().toFloat();
}
*/

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::checkReady()
{
    bool voxelOK;
    voxelsize[0] = ui->lineEdit_xvoxel->text().toFloat();
    voxelsize[1] = ui->lineEdit_yvoxel->text().toFloat();
    voxelsize[2] = ui->lineEdit_zvoxel->text().toFloat();
    margin = ui->lineEditMargin->text().toFloat();

    if (voxelsize[0] > 0 && voxelsize[1] > 0 && voxelsize[2] > 0)
        voxelOK = true;
    else
        voxelOK = false;
    if (voxelOK && margin > 0 && is_am && is_tiff) {
        ready = true;
        ui->pushButtonGo->setEnabled(true);
    } else {
        ready = false;
        ui->pushButtonGo->setEnabled(false);
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
int MainWindow::createTiffData()
{
    int ie, ip, x0, y0, z0, ix, iy, iz, nx, ny, nz, x, y, z;
    float r, r2, dx, dy, dz, wx, wy, wz, d2;
    EDGE edge;
    APOINT p;

    printf("ne: %d\n",network.ne);
    fprintf(fpout,"ne: %d\n",network.ne);
    // First determine the required buffer size to hold the voxels
    wx = 0;
    wy = 0;
    wz = 0;
    for (ie=0; ie<network.ne; ie++) {
        edge = network.edgeList[ie];
        for (ip=0; ip<edge.npts; ip++) {
//            printf("%d %d %d %d\n",ie,edge.npts,ip,edge.pt[ip]);
//            fprintf(fpout,"%d %d %d %d\n",ie,edge.npts,ip,edge.pt[ip]);
//            fflush(fpout);
            p = network.point[edge.pt[ip]];
//            printf("%6.1f %6.1f %6.1f  %6.2f\n",p.x,p.y,p.z,p.d);
//            fprintf(fpout,"%6.1f %6.1f %6.1f  %6.2f\n",p.x,p.y,p.z,p.d);
//            fflush(fpout);
            wx = MAX(wx,(p.x + p.d/2));
            wy = MAX(wy,(p.y + p.d/2));
            wz = MAX(wz,(p.z + p.d/2));
        }
    }
    width = (int)((wx+margin)/voxelsize[0]+10);
    height = (int)((wy+margin)/voxelsize[1]+10);
    depth = (int)((wz+margin)/voxelsize[2]+10);
    xysize = width*height;
    printf("width, height, depth: %d %d %d\n",width, height, depth);
    fprintf(fpout,"width, height, depth: %d %d %d\n",width, height, depth);
    fflush(fpout);
    p_im = (unsigned char *)malloc(width*height*depth*sizeof(unsigned char));
    memset(p_im,0,width*height*depth);
    for (ie=0; ie<network.ne; ie++) {
        edge = network.edgeList[ie];
        for (ip=0; ip<edge.npts; ip++) {
//            fprintf(fpout,"ie, npts,ip: %d %d %d\n",ie,edge.npts,ip);
//            fflush(fpout);
            p = network.point[edge.pt[ip]];
            x0 = p.x/voxelsize[0];   // voxel nearest to the point
            y0 = p.y/voxelsize[1];
            z0 = p.z/voxelsize[2];
            r = p.d/2 + margin;
//            fprintf(fpout,"point: %d  %6.1f %6.1f %6.1f  %d %d %d  %6.2f\n",ip,p.x,p.y,p.z,x0,y0,z0,r);
//            fflush(fpout);
            r2 = r*r;
            nx = r/voxelsize[0] + 1;
            ny = r/voxelsize[1] + 1;
            nz = r/voxelsize[2] + 1;
//            printf("nx,ny,nz,x1,y1,z1: %d %d %d  %d %d %d\n",nx,ny,nz,x1,y1,z1);
            for (ix = -nx; ix <= nx; ix++) {
                dx = ix*voxelsize[0];
                x = x0+ix;
                if (x < 0 || x >= width) continue;
                for (iy = -ny; iy<=ny; iy++) {
                    dy = iy*voxelsize[1];
                    y = y0+iy;
                    if (y < 0 || y >= height) continue;
                    for (iz = -nz; iz<=nz; iz++) {
                        dz = iz*voxelsize[2];
                        z = z0+iz;
                        if (z < 0 || z >= depth) continue;
                        d2 = dx*dx+dy*dy+dz*dz;
                        if (d2 < r2) {
                            V(x0+ix,y0+iy,z0+iz) = 255;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::am_tiffer()
{
    int err;
    QString resultstr, numstr;

    ui->labelResult->setText("");
    if (!ready) {
        ui->labelResult->setText("Not ready: select files");
        return;
    }
    if (!am_read) {
        resultstr = "reading network file...";
        ui->labelResult->setText(resultstr);
        QCoreApplication::processEvents();
        err = readNetwork(&network, amFileName.toAscii().constData());
        if (err != 0) {
            resultstr = "FAILED: read_network";
            ui->labelResult->setText(resultstr);
            return;
        }
    }
    resultstr = "creating image file...";
    ui->labelResult->setText(resultstr);
    QCoreApplication::processEvents();
    err = createTiffData();
    if (err != 0) {
        resultstr = "FAILED: createTiffData";
        ui->labelResult->setText(resultstr);
        return;
    }
    err = createTiff(tiffFileName.toAscii().constData(),p_im,width,height,depth);
    if (err != 0) {
        resultstr = "FAILED: createTiff";
        ui->labelResult->setText(resultstr);
        return;
    }
    checkReady();
    resultstr = "SUCCESS";
    ui->labelResult->setText(resultstr);
    return;
}

