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
    QString infoFile = QCoreApplication::applicationDirPath() + "/info/deadends_info.txt";
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

    fpout = fopen("dead.out","w");
    network.ne = 0;
    network.nv = 0;
    network.np = 0;
    deadlist = NULL;
    ndead = 0;
    isSphere = false;
    is_am = false;
    is_tiff = false;
    is_output = false;
    am_read = false;
    tiff_read = false;
    detected = false;
    evaluated = false;
    voxelsize[0] = voxelsize[1] = voxelsize[2] = 2;
    checkReady();

//      image_view = vtkImageViewer2::New();
//      image_view = vtkSmartPointer<vtkResliceImageViewer>::New();    
//Qt Signal slot
//      connect(ui->actionDICOM_Sequence,SIGNAL(triggered()),this, SLOT(DICOMseq()));
//      connect(ui->actionExit,SIGNAL(triggered()),this, SLOT(close()));
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
    tiffFileName = QFileDialog::getOpenFileName(this,
        tr("Input TIFF file"), ".", tr("TIFF Files (*.tif)"));
    if (tiffFileName.compare(ui->labelTiffFile->text()) != 0) {
        tiff_read = false;
        printf("%s\n",tiffFileName.toAscii().constData());
        printf("%s\n",ui->labelTiffFile->text().toAscii().constData());
    }
    ui->labelTiffFile->setText(tiffFileName);
    is_tiff = true;
    //Testing
//    int err = readTiff(tiffFileName.toAscii().constData(),&width,&height,&depth);
//    if (err != 0) {
//        return;
//    }

    checkReady();
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::outputFileSelecter()
{
    ui->labelResult->setText("");
    outputFileName = QFileDialog::getSaveFileName(this,
        tr("Output Amira file"), ".", tr("Amira Files (*.am)"));
    ui->labelOutputFile->setText(outputFileName);
    is_output = true;
    checkReady();
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::voxelChanged()
{
    checkReady();
}

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

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::checkReady()
{
    bool voxelOK;
    maxRadius = ui->lineEditRlimit->text().toFloat();
    voxelsize[0] = ui->lineEdit_xvoxel->text().toFloat();
    voxelsize[1] = ui->lineEdit_yvoxel->text().toFloat();
    voxelsize[2] = ui->lineEdit_zvoxel->text().toFloat();
    if (voxelsize[0] > 0 && voxelsize[1] > 0 && voxelsize[2] > 0)
        voxelOK = true;
    else
        voxelOK = false;
    if (voxelOK && is_am && is_tiff && is_output) {
        ready = true;
        ui->pushButtonDeadends->setEnabled(true);
        ui->pushButtonEvaluateAll->setEnabled(true);
        if (detected) {
            ui->pushButtonEvaluate->setEnabled(true);
            if (evaluated)
                ui->pushButtonSaveFile->setEnabled(true);
        } else {
            ui->pushButtonEvaluate->setEnabled(false);
            ui->pushButtonSaveFile->setEnabled(false);
        }
    } else {
        ready = false;
        ui->pushButtonDeadends->setEnabled(false);
        ui->pushButtonEvaluateAll->setEnabled(false);
        ui->pushButtonEvaluate->setEnabled(false);
        ui->pushButtonSaveFile->setEnabled(false);
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::detectDeadends()
{
    int err;
    QString resultstr, numstr;

    ui->labelResult->setText("");
    if (!ready) {
        ui->labelResult->setText("Not ready: select files");
        return;
    }
    if (isSphere) {
        getSphere();
    }
    ui->lineEditNdead->setText("");
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
    if (!tiff_read) {
        resultstr = "reading image file...";
        ui->labelResult->setText(resultstr);
        QCoreApplication::processEvents();
        err = readTiff(tiffFileName.toAscii().constData(),&width,&height,&depth);
        if (err != 0) {
            resultstr = "FAILED: read_tiff";
            ui->labelResult->setText(resultstr);
            return;
        }
        xysize = width*height;
    }
    resultstr = "working...";
    ui->labelResult->setText(resultstr);
    QCoreApplication::processEvents();
    err = findDeadends(&network, &deadlist, &ndead);
    if (err != 0) {
        resultstr = "FAILED: find_deadends";
        ui->labelResult->setText(resultstr);
        return;
    }
    printf("ndead: %d\n",ndead);
    fprintf(fpout,"ndead: %d\n",ndead);
    numstr.setNum(ndead);
    ui->lineEditNdead->setText(numstr);
    detected = true;
    checkReady();
    resultstr = "SUCCESS";
    ui->labelResult->setText(resultstr);
    return;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::evaluate()
{
    QString resultstr;
    resultstr = "working...";
    ui->labelResult->setText(resultstr);
    QCoreApplication::processEvents();
    getIntensities(&network, deadlist, ndead);
    evaluated = true;
    checkReady();
    resultstr = "SUCCESS";
    ui->labelResult->setText(resultstr);
}

//----------------------------------------------------------------------------------------
// For each dead end, evaluate the image intensity of the connected vessel (edge).
// What makes this a bit tricky is the fact that the vessel walls are stained.  This
// suggests that it might be more appropriate to determine the peak intensity in the
// vicinity of the vessel, rather than the mean intensity.
// A possible approach:
// (1) look at each point along the edge
// (2) based on the vessel diameter at that point, find all voxels in a specified range of
//     distance from the point, e.g. from 0.5*r to 1.5*r
// (3) determine the maximum intensity of this set of voxels
// (4) record the average of the maximum intensities of all points
//----------------------------------------------------------------------------------------
int MainWindow::getIntensities(NETWORK *net, DEADEND *deadlist, int ndead)
{
    int id, ie, npts, ip, count;
    int nv, v[100000][3], maxval;
    float c[3], r, totval, sum;
    float beta = 0.5;
    bool toobig;
    EDGE edge;
    APOINT p;

    printf("getIntensities: %d\n",ndead);
    fprintf(fpout,"getIntensities: %d\n",ndead);
    count = 0;
    sum = 0;
    for (id=0; id<ndead; id++) {
        ie = deadlist[id].ie;
        edge = net->edgeList[ie];
        npts = edge.npts;
        toobig = false;
        totval = 0;
        for (ip=0; ip<npts; ip++) {
            p = net->point[edge.pt[ip]];
//            printf("%d  %d  %d  %d  %6.1f %6.1f %6.1f %6.2f\n",id,ie,npts,ip,p.x,p.y,p.z,p.d);
//            fprintf(fpout,"%d  %d  %d  %d  %6.1f %6.1f %6.1f %6.2f\n",id,ie,npts,ip,p.x,p.y,p.z,p.d);
            r = p.d/2;
            c[0] = p.x;
            c[1] = p.y;
            c[2] = p.z;
            getVoxels(r,c,beta,&nv,v);
//            printf("id, npts, ip, nv: %d %d %d %d\n",id,npts,ip,nv);
//            fprintf(fpout,"id, npts, ip, nv: %d %d %d %d\n",id,npts,ip,nv);
//            fflush(fpout);
            if (nv == 0) {
//                printf("getVoxels: no voxels: c, r: %6.1f %6.1f %6.1f  %6.2f\n",c[0],c[1],c[2],r);
//                fprintf(fpout,"getVoxels: no voxels: c, r: %6.1f %6.1f %6.1f  %6.2f\n",c[0],c[1],c[2],r);
//                exit(1);
                toobig = true;
//                fprintf(fpout,"Too big\n");
//                fflush(fpout);
                break;
            }
            getMaxIntensity(nv,v,&maxval);
//            printf("maxval: %d\n",maxval);
//            fprintf(fpout,"maxval: %d\n",maxval);
//            fflush(fpout);
            totval += maxval;
//            pv[0] = p.x/voxelsize[0] + 1;
//            pv[1] = p.y/voxelsize[1] + 1;
//            pv[2] = p.z/voxelsize[2] + 1;
        }
        if (toobig) {
            deadlist[id].intensity = 0;
            continue;
        }
        deadlist[id].intensity = totval/npts;
        sum += deadlist[id].intensity;
        count++;
//        printf("intensity: %6d %6.1f\n",id,deadlist[id].intensity);
        fprintf(fpout,"intensity: %6d %6.1f\n",id,deadlist[id].intensity);
        fflush(fpout);
    }
    printf("Average dead-end intensity: %6.1f\n",sum/count);
    fprintf(fpout,"Average dead-end intensity: %6.1f\n",sum/count);
    return 0;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::getAllIntensities()
{
    int id, ie, npts, ip, count, err;
    int nv, v[100000][3], maxval;
    int *vert;
    float c[3], r, totval, intensity, sum;
    float beta = 0.5;
    bool toobig;
    NETWORK net;
    QString resultstr;
    EDGE edge;
    APOINT p;

    printf("getAllIntensities:\n");
    fprintf(fpout,"getAllIntensities:\n");
    if (isSphere) {
        getSphere();
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
    if (!tiff_read) {
        resultstr = "reading image file...";
        ui->labelResult->setText(resultstr);
        QCoreApplication::processEvents();
        err = readTiff(tiffFileName.toAscii().constData(),&width,&height,&depth);
        if (err != 0) {
            resultstr = "FAILED: read_tiff";
            ui->labelResult->setText(resultstr);
            return;
        }
        xysize = width*height;
    }
    resultstr = "working...";
    ui->labelResult->setText(resultstr);
    QCoreApplication::processEvents();
    net = network;
    count = 0;
    sum = 0;
    for (ie=0; ie<net.ne; ie++) {
        edge = net.edgeList[ie];
        if (isSphere) {
            vert = edge.vert;
            if (!inSphere(net.vertex[vert[0]].point) && !inSphere(net.vertex[vert[1]].point)) continue;
        }
        npts = edge.npts;
        toobig = false;
        totval = 0;
        for (ip=0; ip<npts; ip++) {
            p = net.point[edge.pt[ip]];
//            printf("%d  %d  %d  %d  %6.1f %6.1f %6.1f %6.2f\n",id,ie,npts,ip,p.x,p.y,p.z,p.d);
//            fprintf(fpout,"%d  %d  %d  %d  %6.1f %6.1f %6.1f %6.2f\n",id,ie,npts,ip,p.x,p.y,p.z,p.d);
            r = p.d/2;
            c[0] = p.x;
            c[1] = p.y;
            c[2] = p.z;
            getVoxels(r,c,beta,&nv,v);
//            printf("ie, npts, ip, nv: %d %d %d %d\n",ie,npts,ip,nv);
//            fprintf(fpout,"ie, npts, ip, nv: %d %d %d %d\n",ie,npts,ip,nv);
//            fflush(fpout);
            if (nv == 0) {
                toobig = true;
                break;
            }
            getMaxIntensity(nv,v,&maxval);
//            printf("maxval: %d %d\n",ie,maxval);
//            fprintf(fpout,"maxval: %d %d\n",ie,maxval);
//            fflush(fpout);
            totval += maxval;
        }
        if (toobig) {
//            printf("Too big\n");
            continue;
        }
        intensity = totval/npts;
        sum += intensity;
        count++;
        printf("intensity: %6d %6.1f\n",ie,intensity);
        fprintf(fpout,"intensity: %6d %6.1f\n",ie,intensity);
        fflush(fpout);
    }
    printf("Average vessel intensity: %6.1f\n",sum/count);
    fprintf(fpout,"Average vessel intensity: %6.1f\n",sum/count);
    resultstr = "SUCCESS";
    ui->labelResult->setText(resultstr);
    return;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::getVoxels(float r, float c[], float beta, int *nv, int v[][3])
{
    int i, nto[3], ix, iy, iz, k, u[3];
    float dx, dy, dz, d2, rfr2, rto2;

//    rlimit = ui->lineEditRlimit->text().toFloat();

    if (r > maxRadius) {
        *nv = 0;
        return;
    }
    rfr2 = (1-beta)*(1-beta)*r*r;
    rto2 = (1+beta)*(1+beta)*r*r;
    for (i=0; i<3; i++) {
        nto[i] = (1 + beta)*r/voxelsize[i] + 1;
    }
//    fprintf(fpout,"r, beta, rfr2, rto2: %f %f %f %f  %d %d %d\n",r,beta,rfr2,rto2,nto[0],nto[1],nto[2]);
    k = 0;
    for (ix=-nto[0]; ix<=nto[0]; ix++) {
        dx = ix*voxelsize[0];
        for (iy=-nto[1]; iy<=nto[1]; iy++) {
            dy = iy*voxelsize[1];
            for (iz=-nto[2]; iz<=nto[2]; iz++) {
                dz = iz*voxelsize[2];
                d2 = dx*dx+dy*dy+dz*dz;
                if (d2<rfr2 || d2 > rto2) continue;
                u[0] = (c[0] + dx)/voxelsize[0];
                u[1] = (c[1] + dy)/voxelsize[1];
                u[2] = (c[2] + dz)/voxelsize[2];
//                fprintf(fpout,"%4d  %6d %6d %6d\n",k,u[0],u[1],u[2]);
                for (i=0; i<3; i++)
                    v[k][i] = u[i];
                k++;
            }
        }
    }
    if (k == 0) {
        fprintf(fpout,"r, beta, rfr2, rto2: %f %f %f %f\n",r,beta,rfr2,rto2);
        fprintf(fpout,"nto: %d %d %d\n",nto[0],nto[1],nto[2]);
    }
    *nv = k;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::getMaxIntensity(int nv, int v[][3], int *maxval)
{
    int k, x, y, z;
    float val;

    *maxval = 0;
    for (k=0; k<nv; k++) {
        x = v[k][0];
        y = v[k][1];
        z = v[k][2];
        if (x < 0 || y < 0 || z < 0) {
            printf("bad index < 0:  %d %d %d\n",x,y,z);
            fprintf(fpout,"bad index < 0:  %d %d %d\n",x,y,z);
            continue;
        }
        if (x >= width || y >= height || z >= depth) {
            printf("bad index too big:  %d %d %d\n",x,y,z);
            fprintf(fpout,"bad index too big:  %d %d %d\n",x,y,z);
            continue;
        }
        val = V(x,y,z);
        if (val > *maxval)
            *maxval = val;
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
void MainWindow::saveDeadendFile()
{
    createNetwork(&network, &deadnetwork, deadlist, ndead);
    writeAmiraFile(outputFileName.toAscii().constData(), amFileName.toAscii().constData(), &deadnetwork);
}



// Need to add:
// (1) Input for voxelsize
// (2) Optional specification of a spherical region
// (3) Ability to compute average intensity (or distribution) over the whole network.

/*
int MainWindow::DICOMseq()
{
      typedef signed short InputPixelType;//Pixel Type
      const unsigned int InputDimension = 3;//Dimension of image
      typedef itk::Image< InputPixelType, InputDimension > InputImageType;//Image Type

      typedef itk::ImageSeriesReader< InputImageType > ReaderType;//Reader of Image Type
      ReaderType::Pointer reader = ReaderType::New();

      //Replacement of name generator of ITK
      QDir dir("Dir");
      dir=QFileDialog::getExistingDirectory(0,"Select Folder: ");
      QFileInfoList list = dir.entryInfoList(QDir::Dirs| QDir::Files | QDir::NoDotAndDotDot);

            std::vector<std::string> names;
         foreach(QFileInfo finfo, list)
                  {
                        std::string str=dir.path().toStdString().c_str();
                        str=str+"/";
                        names.push_back(str+finfo.fileName().toStdString().c_str());
                  }

      reader->SetFileNames( names );
      //Exceptional handling
      try
      {
            reader->Update();
      }
      catch (itk::ExceptionObject & e)
      {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e << std::endl;
            return EXIT_FAILURE;
      }


     
      You can impliment your Filter here and connect its output to connector below.
     


      //connector to convert ITK image data to VTK image data
      typedef itk::ImageToVTKImageFilter<InputImageType> ConnectorType;
      ConnectorType::Pointer connector= ConnectorType::New();
      connector->SetInput( reader->GetOutput() );//Set ITK reader Output to connector you can replace it with filter
     
      //Exceptional handling
      try
      {
            connector->Update();
      }
      catch (itk::ExceptionObject & e)
      {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e << std::endl;
            return EXIT_FAILURE;
      }
      //deep copy connector's output to an image else connector will go out of scope
      //and vanish it will cause error while changing slice
      vtkImageData * image = vtkImageData::New();
      image->DeepCopy(connector->GetOutput());
     
      //set VTK Viewer to QVTKWidget in Qt's UI
      ui->qVTK1->SetRenderWindow(image_view->GetRenderWindow());
      image_view->SetupInteractor(ui->qVTK1->GetRenderWindow()->GetInteractor());
      //Set input image to VTK viewer
      image_view->SetInput(image);
      image_view->SetSlice(image_view->GetSliceMax()/2);
      image_view->GetRenderer()->ResetCamera();
      image_view->Render();


    ui->qVTK1->update();

      return EXIT_SUCCESS;
}

*/
