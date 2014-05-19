#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QFile>
#include <QTextStream>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
	ui->textEdit->setReadOnly(true);
    QString infoFile = QCoreApplication::applicationDirPath() + "/info/quantifier_info.txt";
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
    axisstr[0] = "X";
    axisstr[1] = "Y";
    axisstr[2] = "Z";

    imageViewer = NULL;
    fpout = NULL;
    is_ready = false;
    is_amfile = false;
    is_closefile = false;
    is_resultfile = false;
    ui->pushButton_vessels->setEnabled(false);
    ui->groupBox_slice->setEnabled(false);
    ui->checkBox_selectblock->setEnabled(false);
    is_block = ui->checkBox_selectblock->isChecked();
    checkReady();
    reset();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::amFileSelecter()
{
	ui->labelResult->setText("");
    amFileName = QFileDialog::getOpenFileName(this,
        tr("Input Amira file"), ".", tr("Amira Files (*.am)"));
    ui->labelAmFile->setText(amFileName);
    is_amfile = true;
    checkReady();
    reset();
}

void MainWindow::closeFileSelecter()
{
	ui->labelResult->setText("");
    closeFileName = QFileDialog::getOpenFileName(this,
        tr("Input close image binary file"), ".", tr("Bin Files (*.bin)"));
    ui->labelCloseFile->setText(closeFileName);
    is_closefile = true;
    checkReady();
    reset();
}

void MainWindow::resultFileSelecter()
{
    ui->labelResult->setText("");
    resultFileName = QFileDialog::getSaveFileName(this,
        tr("Result file"), ".", tr("Result Files (*.out)"));
    ui->labelResultFile->setText(resultFileName);
    is_resultfile = true;
    checkReady();
    reset();
}

void MainWindow::voxelChanged()
{
//    reset();
    checkReady();
    if (isSetup()) {
        computeVolume();
    }
}

void MainWindow::on_radioButton_slice_toggled(bool checked)
{
    is_slice = checked;
    is_average = !checked;
    ui->lineEdit_intercept->setEnabled(is_slice);
    ui->groupBox_centre->setEnabled(is_block && ui->radioButton_centre->isChecked());
    ui->groupBox_range->setEnabled(is_block && ui->radioButton_range->isChecked());
    ui->groupBox_slice->setEnabled(is_slice);
    ui->groupBox_average->setEnabled(is_average);
    ui->checkBox_selectblock->setChecked(false);
    ui->groupBoxRanges->setEnabled(false);
    if (is_block) {
        setBlockAxes();
    }
}

void MainWindow::on_radioButton_centre_toggled(bool checked)
{
    ui->groupBox_centre->setEnabled(checked);
    ui->groupBox_range->setEnabled(!checked);
    if (is_block) {
        setBlockAxes();
    }
}

void MainWindow::on_checkBox_selectblock_toggled(bool checked)
{
    is_block = checked;
    ui->groupBox_centre->setEnabled(is_block && ui->radioButton_centre->isChecked());
    ui->groupBox_range->setEnabled(is_block && ui->radioButton_range->isChecked());
    if (is_block) {
        setBlockAxes();
        ui->groupBoxRanges->setEnabled(true);
    }
}

void MainWindow::on_radioButton_xaxis_toggled(bool checked)
{
    if (is_block) {
        setBlockAxes();
    }
}

void MainWindow::on_radioButton_yaxis_toggled(bool checked)
{
    if (is_block) {
        setBlockAxes();
    }
}

void MainWindow::on_radioButton_zaxis_toggled(bool checked)
{
    if (is_block) {
        setBlockAxes();
    }
}

void MainWindow :: setBlockAxes(){
    // Adjust status of axis range selection when there is a change to the slice or average axes
    if (!is_block) return;
    enableRange('X',true);
    enableRange('Y',true);
    enableRange('Z',true);
    if (is_slice) {
        if (ui->radioButton_xaxis->isChecked()) {
            enableRange('X',false);
        }
        else if (ui->radioButton_yaxis->isChecked()) {
            enableRange('Y',false);
        }
        else if (ui->radioButton_zaxis->isChecked()) {
            enableRange('Z',false);
        }
    }
}

void MainWindow :: enableRange(char axischar, bool enable)
{
    if (ui->radioButton_range->isChecked()) {
        if (axischar == 'X') {
            ui->checkBox_xfull->setEnabled(enable);
            ui->lineEdit_x1->setEnabled(enable && !ui->checkBox_xfull->isChecked());
            ui->lineEdit_x2->setEnabled(enable && !ui->checkBox_xfull->isChecked());
        }
        else if (axischar == 'Y') {
            ui->checkBox_yfull->setEnabled(enable);
            ui->lineEdit_y1->setEnabled(enable && !ui->checkBox_yfull->isChecked());
            ui->lineEdit_y2->setEnabled(enable && !ui->checkBox_yfull->isChecked());
        }
        if (axischar == 'Z') {
            ui->checkBox_zfull->setEnabled(enable);
            ui->lineEdit_z1->setEnabled(enable && !ui->checkBox_zfull->isChecked());
            ui->lineEdit_z2->setEnabled(enable && !ui->checkBox_zfull->isChecked());
        }

    } else {
        if (axischar == 'X') {
            ui->lineEdit_xcentre->setEnabled(enable);
            ui->lineEdit_xsize->setEnabled(enable);
        }
        else if (axischar == 'Y') {
            ui->lineEdit_ycentre->setEnabled(enable);
            ui->lineEdit_ysize->setEnabled(enable);
        }
        else if (axischar == 'Z') {
            ui->lineEdit_zcentre->setEnabled(enable);
            ui->lineEdit_zsize->setEnabled(enable);
        }
    }
}

void MainWindow::checkBox_x()
{
    if (ui->checkBox_xfull->isChecked()) {
        ui->lineEdit_x1->setEnabled(false);
        ui->lineEdit_x2->setEnabled(false);
    } else {
        ui->lineEdit_x1->setEnabled(true);
        ui->lineEdit_x2->setEnabled(true);
    }
}

void MainWindow::checkBox_y()
{
    if (ui->checkBox_yfull->isChecked()) {
        ui->lineEdit_y1->setEnabled(false);
        ui->lineEdit_y2->setEnabled(false);
    } else {
        ui->lineEdit_y1->setEnabled(true);
        ui->lineEdit_y2->setEnabled(true);
    }
}

void MainWindow::checkBox_z()
{
    if (ui->checkBox_zfull->isChecked()) {
        ui->lineEdit_z1->setEnabled(false);
        ui->lineEdit_z2->setEnabled(false);
    } else {
        ui->lineEdit_z1->setEnabled(true);
        ui->lineEdit_z2->setEnabled(true);
    }
}

void MainWindow::checkReady()
{
    bool voxelOK;
    voxelsize[0] = ui->lineEdit_xvoxel->text().toDouble();
    voxelsize[1] = ui->lineEdit_yvoxel->text().toDouble();
    voxelsize[2] = ui->lineEdit_zvoxel->text().toDouble();
    if (voxelsize[0] > 0 && voxelsize[1] > 0 && voxelsize[2] > 0)
        voxelOK = true;
    else
        voxelOK = false;
    if ((voxelOK || isSetup()) && is_amfile && is_closefile && is_resultfile) {
        is_ready = true;
        ui->pushButton_setup->setEnabled(true);
    } else {
        is_ready = false;
        ui->pushButton_setup->setEnabled(false);
        ui->groupBox_slice->setEnabled(false);
        ui->pushButton_vessels->setEnabled(false);
        ui->checkBox_selectblock->setEnabled(false);
    }
}

void MainWindow::reader()
{
    doSetup();
    if (isSetup()) {
        ui->groupBox_slice->setEnabled(true);
        ui->pushButton_vessels->setEnabled(true);
        ui->checkBox_selectblock->setEnabled(true);
    } else {
        is_ready = false;
        ui->pushButton_setup->setEnabled(false);
        ui->groupBox_slice->setEnabled(false);
        ui->pushButton_vessels->setEnabled(false);
        ui->checkBox_selectblock->setEnabled(false);
    }
}

void MainWindow::computeVolume()
{
    int err;
    int ntvoxels;
    double volume;
    QString volumestr, countstr;
    bool voxelOK;

    if (!is_ready) {
        return;
    }
    if (!isSetup()) {
        return;
    }
    voxelsize[0] = ui->lineEdit_xvoxel->text().toDouble();
    voxelsize[1] = ui->lineEdit_yvoxel->text().toDouble();
    voxelsize[2] = ui->lineEdit_zvoxel->text().toDouble();
    if (voxelsize[0] > 0 && voxelsize[1] > 0 && voxelsize[2] > 0)
        voxelOK = true;
    else
        voxelOK = false;
    if (!voxelOK) return;
    ntvoxels = TotalVoxelCount();
    countstr = QString::number(ntvoxels);
    ui->lineEdit_ntvoxels->setText(countstr);
    fprintf(fpout,"Total voxel count: %d\n",ntvoxels);
    volume = ntvoxels*voxelsize[0]*voxelsize[1]*voxelsize[2];
    volume = 1.0e-9*volume;   // convert um3 -> mm3
    volumestr = QString::number(volume,'f',3);
    ui->lineEdit_volume->setText(volumestr);
    fprintf(fpout,"Total volume (mm3): %f\n",(float)volume);
}

void MainWindow::computeArea()
{
    int err;
    double area;
    int axis, islice, npixels;
    QString areastr;

    if (!is_ready) {
        return;
    }
    if (ui->radioButton_xaxis->isChecked())
        axis = 0;
    else if (ui->radioButton_yaxis->isChecked())
        axis = 1;
    else if (ui->radioButton_zaxis->isChecked())
        axis = 2;
    islice = ui->lineEdit_intercept->text().toInt();
    if (ui->radioButton_unitsmicrons->isChecked())
        islice = (int)(islice/voxelsize[axis]);
    err = checkSlice(axis,islice);
    if (err == 0)
        err = getArea(axis,islice,&npixels,&area);
    else
        area = 0;
    area = 1.0e-6*area;   // convert um2 -> mm2
    areastr = QString::number(area,'f',3);
    ui->lineEdit_area->setText(areastr);
    ui->lineEdit_count->setText("");
    ui->lineEdit_density->setText("");
}

void MainWindow::getCentredRanges()
{
    double x0, y0, z0, xw, yw, zw, xr, yr, zr, xfac, yfac, zfac;
    QString x0str, y0str, z0str, xsizestr, ysizestr, zsizestr;

    x0str = ui->lineEdit_xcentre->text();
    y0str = ui->lineEdit_ycentre->text();
    z0str = ui->lineEdit_zcentre->text();
    xsizestr = ui->lineEdit_xsize->text();
    ysizestr = ui->lineEdit_ysize->text();
    zsizestr = ui->lineEdit_zsize->text();
    x0 = x0str.toDouble();
    y0 = y0str.toDouble();
    z0 = z0str.toDouble();
    xw = xsizestr.toDouble();
    yw = ysizestr.toDouble();
    zw = zsizestr.toDouble();
    xr = xw/2;
    yr = yw/2;
    zr = zw/2;
    if (ui->radioButton_unitsvoxels->isChecked()) {
        xfac = 1;
        yfac = 1;
        zfac = 1;
    } else {
        xfac = 1/voxelsize[0];
        yfac = 1/voxelsize[1];
        zfac = 1/voxelsize[2];
    }
    range[0][0] = (int)(xfac*(x0 - xr));
    range[0][1] = (int)(xfac*(x0 + xr));
    range[1][0] = (int)(yfac*(y0 - yr));
    range[1][1] = (int)(yfac*(y0 + yr));
    range[2][0] = (int)(zfac*(z0 - zr));
    range[2][1] = (int)(zfac*(z0 + zr));
}

//---------------------------------------------------------
// Note: range[][] is 1-based!
//---------------------------------------------------------
void MainWindow::getAveragingRanges() {
    double x1, y1, z1, x2, y2, z2, xfac, yfac, zfac;
    QString x1str, y1str, z1str, x2str, y2str, z2str;

    if (ui->radioButton_unitsvoxels->isChecked()) {
        xfac = 1;
        yfac = 1;
        zfac = 1;
    } else {
        xfac = 1/voxelsize[0];
        yfac = 1/voxelsize[1];
        zfac = 1/voxelsize[2];
    }
    if (DEBUG) fprintf(fpout,"Voxel factors: %f %f %f\n",xfac,yfac,zfac);
    if (!ui->checkBox_xfull->isChecked()) {
        x1str = ui->lineEdit_x1->text();
        x2str = ui->lineEdit_x2->text();
        x1 = x1str.toDouble();
        x2 = x2str.toDouble();
        range[0][0] = (int)(xfac*x1);
        range[0][1] = (int)(xfac*x2);
    } else {
        range[0][0] = 1;
        range[0][1] = nvoxels[0];
    }
    if (!ui->checkBox_yfull->isChecked()) {
        y1str = ui->lineEdit_y1->text();
        y2str = ui->lineEdit_y2->text();
        y1 = y1str.toDouble();
        y2 = y2str.toDouble();
        range[1][0] = (int)(yfac*y1);
        range[1][1] = (int)(yfac*y2);
    } else {
        range[1][0] = 1;
        range[1][1] = nvoxels[1];
    }
    if (!ui->checkBox_zfull->isChecked()) {
        z1str = ui->lineEdit_z1->text();
        z2str = ui->lineEdit_z2->text();
        z1 = z1str.toDouble();
        z2 = z2str.toDouble();
        range[2][0] = (int)(zfac*z1);
        range[2][1] = (int)(zfac*z2);
    } else {
        range[2][0] = 1;
        range[2][1] = nvoxels[2];
    }
}

void MainWindow::getRanges()
{
    int i;
    if (ui->radioButton_range->isChecked()) {
        if (DEBUG) fprintf(fpout,"getAveragingRanges\n");
        getAveragingRanges();
    } else {
        if (DEBUG) fprintf(fpout,"getCentredRanges\n");
        getCentredRanges();
    }
    checkRanges();
    fprintf(fpout,"Block ranges:\n");
    if (is_slice) {
        if (!ui->radioButton_xaxis->isChecked()) {
            i = 0;
            fprintf(fpout,"axis: %s  %d %d\n",axisstr[i],range[i][0],range[i][1]);
        }
        if (!ui->radioButton_yaxis->isChecked()) {
            i = 1;
            fprintf(fpout,"axis: %s  %d %d\n",axisstr[i],range[i][0],range[i][1]);
        }
        if (!ui->radioButton_zaxis->isChecked()) {
            i = 2;
            fprintf(fpout,"axis: %s  %d %d\n",axisstr[i],range[i][0],range[i][1]);
        }
    } else {
        for (i=0; i<3; i++)
            fprintf(fpout,"axis: %s  %d %d\n",axisstr[i],range[i][0],range[i][1]);
        }
}

void MainWindow::checkRanges()
{
    int axis;

    for (axis=0; axis<3; axis++) {
        range[axis][0] = MAX(range[axis][0],1);
        range[axis][1] = MIN(range[axis][1],nvoxels[axis]);
    }
}

void MainWindow::ComputeVessels()
{
    int err;
    double tissuearea[3], darea, areafraction;
    double dmin, dmax, vessellength_mm, tissuevolume_mm3, mm_mm3;
    int axis, islice, nvessels[3], nvesselpixels[3], ntissuepixels[3], density, MVD, w, h, nbranchpts;
    int nslices;
    bool use_axis[3];
    QString areastr, vesselpixstr, tissuepixstr, countstr, densitystr, MVDstr, fractionstr, mm_mm3str;

    if (DEBUG) fprintf(fpout,"ComputeVessels\n");
    fflush(fpout);
    this->setCursor( QCursor( Qt::WaitCursor ) );

    if (!is_ready) {
        return;
    }
    if (is_slice) {
        ui->lineEdit_mm_mm3->clear();
        use_axis[0] = true;
        use_axis[1] = false;
        use_axis[2] = false;
        if (ui->radioButton_xaxis->isChecked()) {
            axis = 0;
            w = nvoxels[1];
            h = nvoxels[2];
            darea = voxelsize[1]*voxelsize[2];
        } else if (ui->radioButton_yaxis->isChecked()) {
            axis = 1;
            w = nvoxels[2];
            h = nvoxels[0];
            darea = voxelsize[0]*voxelsize[2];
        } else if (ui->radioButton_zaxis->isChecked()) {
            axis = 2;
            w = nvoxels[0];
            h = nvoxels[1];
            darea = voxelsize[0]*voxelsize[1];
        }
        islice = ui->lineEdit_intercept->text().toInt();
        fprintf(fpout,"\n------------------------------------------------------------------------------\n");
        fprintf(fpout,"Computing histology for a slice: axis: %s  islice: %d\n",axisstr[axis],islice);
        fflush(fpout);
        if (is_block) {
            getRanges();
        }
        fprintf(fpout,"------------------------------------------------------------------------------\n");
        fflush(fpout);
        if (ui->radioButton_unitsmicrons->isChecked())
            islice = (int)(islice/voxelsize[axis]);
        err = checkSlice(axis,islice);
        if (err == 0) {
            if (imageViewer) delete imageViewer;
            imageViewer = new ImageViewer(w,h);
            err = SliceHistology(axis, islice, nvessels, nvesselpixels, ntissuepixels, tissuearea);
            imageViewer->paintLabel();
            imageViewer->show();
        } else {
            nvessels[0] = 0;
            nvesselpixels[0] = 0;
            tissuearea[0] = 0;
            nvessels[0] = 0;
        }
    } else {
        fprintf(fpout,"\n------------------------------------------------------------------------------\n");
        fprintf(fpout,"Computing volume average histology\n");
        fflush(fpout);
        use_axis[0] = ui->checkBox_xaxis->isChecked();
        use_axis[1] = ui->checkBox_yaxis->isChecked();
        use_axis[2] = ui->checkBox_zaxis->isChecked();
        ui->lineEdit_area_x->clear();
        ui->lineEdit_areapixels_x->clear();
        ui->lineEdit_count_x->clear();
        ui->lineEdit_density_x->clear();
        ui->lineEdit_MVD_x->clear();
        ui->lineEdit_vesselpixels_x->clear();
        ui->lineEdit_fraction_x->clear();
        ui->lineEdit_area_y->clear();
        ui->lineEdit_areapixels_y->clear();
        ui->lineEdit_count_y->clear();
        ui->lineEdit_density_y->clear();
        ui->lineEdit_MVD_y->clear();
        ui->lineEdit_vesselpixels_y->clear();
        ui->lineEdit_fraction_y->clear();
        ui->lineEdit_area_z->clear();
        ui->lineEdit_areapixels_z->clear();
        ui->lineEdit_count_z->clear();
        ui->lineEdit_density_z->clear();
        ui->lineEdit_MVD_z->clear();
        ui->lineEdit_vesselpixels_z->clear();
        ui->lineEdit_fraction_z->clear();
        getRanges();
        fprintf(fpout,"------------------------------------------------------------------------------\n");
        fflush(fpout);
        err = VolumeHistology(use_axis, nvessels, nvesselpixels, ntissuepixels, tissuearea);
        dmin = ui->lineEdit_diam_min->text().toDouble();
        dmax = ui->lineEdit_diam_max->text().toDouble();
        if (dmax == 0) {
            dmax = 1.0e10;
        }
        VesselDensity(dmin, dmax, &vessellength_mm, &nbranchpts, &tissuevolume_mm3);
        mm_mm3 = vessellength_mm/tissuevolume_mm3;
        mm_mm3str = QString::number(mm_mm3,'f',1);
        ui->lineEdit_mm_mm3->setText(mm_mm3str);
//        err = branching(&nbranchpts, &totlen, &totvol);
    }
    for (axis=0; axis<3; axis++) {
        if (!use_axis[axis]) continue;
        areafraction = (double)nvesselpixels[axis]/ntissuepixels[axis];
        tissuearea[axis] = 1.0e-6*tissuearea[axis];   // convert um2 -> mm2
        areastr = QString::number(tissuearea[axis],'f',3);
        tissuepixstr = QString::number(ntissuepixels[axis]);
        countstr = QString::number(nvessels[axis]);
        if (tissuearea[axis] > 0) {
            density = (int)(nvessels[axis]/tissuearea[axis] + 0.5);
            MVD = (int)((nvessels[axis]*0.74/tissuearea[axis]) + 0.5);
        } else {
            density = 0;
            MVD = 0;
        }
        densitystr = QString::number(density);
        MVDstr = QString::number(MVD);
        vesselpixstr = QString::number(nvesselpixels[axis]);
        fractionstr = QString::number(100*areafraction,'f',2);
        if (is_slice) {
            ui->lineEdit_area->setText(areastr);
            ui->lineEdit_areapixels->setText(tissuepixstr);
            ui->lineEdit_count->setText(countstr);
            ui->lineEdit_density->setText(densitystr);
            ui->lineEdit_MVD->setText(MVDstr);
            ui->lineEdit_vesselpixels->setText(vesselpixstr);
            ui->lineEdit_fraction->setText(fractionstr);
            fprintf(fpout,"Slice area (mm2): %f\n",(float)tissuearea[0]);
            fprintf(fpout,"Slice pixels    : %8d\n",ntissuepixels[0]);
            fprintf(fpout,"Vessel count:     %6d\n",nvessels[0]);
            fprintf(fpout,"Count/area:       %6d\n",density);
            fprintf(fpout,"MVD:              %6d\n",MVD);
            fprintf(fpout,"Vessel pixels:    %6d\n",nvesselpixels[0]);
            fprintf(fpout,"Area percentage:  %8.3f\n",100*areafraction);
        } else {
            fprintf(fpout,"\nAverages for slices normal to %s axis\n",axisstr[axis]);
            if (axis == 0) {
                ui->lineEdit_area_x->setText(areastr);
                ui->lineEdit_areapixels_x->setText(tissuepixstr);
                ui->lineEdit_count_x->setText(countstr);
                ui->lineEdit_density_x->setText(densitystr);
                ui->lineEdit_MVD_x->setText(MVDstr);
                ui->lineEdit_vesselpixels_x->setText(vesselpixstr);
                ui->lineEdit_fraction_x->setText(fractionstr);
            }
            if (axis == 1) {
                ui->lineEdit_area_y->setText(areastr);
                ui->lineEdit_areapixels_y->setText(tissuepixstr);
                ui->lineEdit_count_y->setText(countstr);
                ui->lineEdit_density_y->setText(densitystr);
                ui->lineEdit_MVD_y->setText(MVDstr);
                ui->lineEdit_vesselpixels_y->setText(vesselpixstr);
                ui->lineEdit_fraction_y->setText(fractionstr);
            }
            if (axis == 2) {
                ui->lineEdit_area_z->setText(areastr);
                ui->lineEdit_areapixels_z->setText(tissuepixstr);
                ui->lineEdit_count_z->setText(countstr);
                ui->lineEdit_density_z->setText(densitystr);
                ui->lineEdit_MVD_z->setText(MVDstr);
                ui->lineEdit_vesselpixels_z->setText(vesselpixstr);
                ui->lineEdit_fraction_z->setText(fractionstr);
            }
            nslices = range[axis][1] - range[axis][0] + 1;
            fprintf(fpout,"Number of slices: %d\n",nslices);
            fprintf(fpout,"Slice area (mm2): %f\n",(float)tissuearea[axis]);
            fprintf(fpout,"Slice pixels:     %8d\n",ntissuepixels[axis]);
            fprintf(fpout,"Vessel count:     %6d\n",nvessels[axis]);
            fprintf(fpout,"Count/area:       %6d\n",density);
            fprintf(fpout,"MVD:              %6d\n",MVD);
            fprintf(fpout,"Vessel pixels:    %6d\n",nvesselpixels[axis]);
            fprintf(fpout,"Area percentage:  %8.3f\n",100*areafraction);
        }
    }
    fflush(fpout);
    this->setCursor( QCursor( Qt::ArrowCursor ) );
}

int MainWindow::checkSlice(int axis, int islice)
{
    if (islice < 1 || islice > nvoxels[axis]) {
        // need a message to the user
        return 1;
    } else {
        return 0;
    }
}

void MainWindow::doSetup()
{
    int err;
    int um[3];
    QString numstr, resultstr;
    char input_amfile[1024], close_file[1024], result_file[1024];

    this->setCursor( QCursor( Qt::WaitCursor ) );
    ui->labelResult->setText("Doing setup()");
    strcpy(input_amfile, amFileName.toAscii().constData());
    strcpy(close_file, closeFileName.toAscii().constData());
    strcpy(result_file, resultFileName.toAscii().constData());
    if (fpout) fclose(fpout);
    fpout = fopen(result_file,"w");
    if (!isSetup()) {
        err = setup(input_amfile,close_file,result_file,voxelsize);
        resultstr = QString::number(err);
//        ui->labelResult->setText(resultstr);
        if (err != 0) {
            this->setCursor( QCursor( Qt::ArrowCursor ) );
            return;
        }
    }
    getCloseSize(nvoxels);
    numstr = QString::number(nvoxels[0]);
    ui->lineEdit_nx->setText(numstr);
    numstr = QString::number(nvoxels[1]);
    ui->lineEdit_ny->setText(numstr);
    numstr = QString::number(nvoxels[2]);
    ui->lineEdit_nz->setText(numstr);
    um[0] = (int)(nvoxels[0]*voxelsize[0] + 0.5);
    numstr = QString::number(um[0]);
    ui->lineEdit_umx->setText(numstr);
    um[1] = (int)(nvoxels[1]*voxelsize[1] + 0.5);
    numstr = QString::number(um[1]);
    ui->lineEdit_umy->setText(numstr);
    um[2] = (int)(nvoxels[2]*voxelsize[2] + 0.5);
    numstr = QString::number(um[2]);
    ui->lineEdit_umz->setText(numstr);
    fprintf(fpout,"\nVoxel size: %6.2f x %6.2f x %6.2f\n",voxelsize[0],voxelsize[1], voxelsize[2]);
    fprintf(fpout,"Image size (voxels): %6d x %6d x %6d\n",nvoxels[0],nvoxels[1], nvoxels[2]);
    fprintf(fpout,"Image size (um): %6d x %6d x %6d\n",um[0],um[1],um[2]);
    computeVolume();
    this->setCursor( QCursor( Qt::ArrowCursor ) );
    //	res = system(cmdstr);
//	if (res == 0)
//		resultstr = "SUCCESS";
//	else if (res == 1)
//		resultstr = "FAILED: wrong number of arguments";
//	else if (res == 2)
//		resultstr = "FAILED: Read error on input file";
//	else if (res == 3)
//		resultstr = "FAILED: Write error on output file";
//	else if (res == 4)
//		resultstr = "FAILED: out of memory";
//	ui->labelResult->setText(resultstr);
}

