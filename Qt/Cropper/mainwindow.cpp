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
	QString infoFile = QCoreApplication::applicationDirPath() + "/info/cropper_info.txt";
	QFile file(infoFile);
	bool ok = file.open(QIODevice::ReadOnly | QIODevice::Text);
	if (!ok) {
		ui->textEdit->append("The information file is missing:");
		ui->textEdit->append(infoFile);
		return;
	}
	QTextStream in(&file);
	QString line = in.readLine();
	while (!line.isNull()) {
		ui->textEdit->append(line);
		line = in.readLine();
	}
	file.close();
	ui->textEdit->moveCursor(QTextCursor::Start);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::inputFileSelecter()
{
	ui->labelResult->setText("");
	inputFileName = QFileDialog::getOpenFileName(this,
		tr("Input TIFF file"), ".", tr("TIFF Files (*.tif)"));
	ui->labelInputFile->setText(inputFileName);
}

void MainWindow::outputFileSelecter()
{
	ui->labelResult->setText("");
	outputFileName = QFileDialog::getSaveFileName(this,
		tr("Output TIFF file"), ".", tr("TIFF Files (*.tif)"));
	ui->labelOutputFile->setText(outputFileName);
}

void MainWindow::on_radioButton_block_toggled(bool checked)
{
    ui->groupBox_blockchoice->setEnabled(checked);
    bool use_range = ui->radioButton_range->isChecked();
    ui->groupBox_range->setEnabled(checked&&use_range);
    ui->groupBox_centre->setEnabled(checked&&!use_range);
    ui->groupBox_sphere->setEnabled(!checked);
}

void MainWindow::checkBox_x()
{
	if (ui->checkBox_x->isChecked()) {
		ui->lineEdit_x1->setEnabled(false);
		ui->lineEdit_x2->setEnabled(false);
	} else {
		ui->lineEdit_x1->setEnabled(true);
		ui->lineEdit_x2->setEnabled(true);
	}
}

void MainWindow::checkBox_y()
{
	if (ui->checkBox_y->isChecked()) {
		ui->lineEdit_y1->setEnabled(false);
		ui->lineEdit_y2->setEnabled(false);
	} else {
		ui->lineEdit_y1->setEnabled(true);
		ui->lineEdit_y2->setEnabled(true);
	}
}

void MainWindow::checkBox_z()
{
	if (ui->checkBox_z->isChecked()) {
		ui->lineEdit_z1->setEnabled(false);
		ui->lineEdit_z2->setEnabled(false);
	} else {
		ui->lineEdit_z1->setEnabled(true);
		ui->lineEdit_z2->setEnabled(true);
	}
}

void MainWindow::cropper()
{
	int res;
	QString qstr, resultstr;
	char cmdstr[256];

    if (ui->radioButton_range->isChecked()) {
        if (ui->checkBox_x->isChecked()) {
            x1str = "-1";
            x2str = "-1";
        } else {
            x1str = ui->lineEdit_x1->text();
            x2str = ui->lineEdit_x2->text();
        }
        if (ui->checkBox_y->isChecked()) {
            y1str = "-1";
            y2str = "-1";
        } else {
            y1str = ui->lineEdit_y1->text();
            y2str = ui->lineEdit_y2->text();
        }
        if (ui->checkBox_z->isChecked()) {
            z1str = "-1";
            z2str = "-1";
        } else {
            z1str = ui->lineEdit_z1->text();
            z2str = ui->lineEdit_z2->text();
        }
    } else {
        makeRanges();
    }

	qstr = QCoreApplication::applicationDirPath() + "/exec/crop ";
	qstr += inputFileName;
	qstr += " ";
	qstr += outputFileName;
	qstr += " ";
    if (ui->radioButton_block->isChecked()) {
        qstr += x1str;
        qstr += " ";
        qstr += x2str;
        qstr += " ";
        qstr += y1str;
        qstr += " ";
        qstr += y2str;
        qstr += " ";
        qstr += z1str;
        qstr += " ";
        qstr += z2str;
    } else {
        qstr += ui->lineEdit_x0->text();
        qstr += " ";
        qstr += ui->lineEdit_y0->text();
        qstr += " ";
        qstr += ui->lineEdit_z0->text();
        qstr += " ";
        qstr += ui->lineEdit_R->text();
    }

	if (qstr.size()>(int)sizeof(cmdstr)-1) {
		printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
		resultstr = "FAILED: cmdstr not big enough for the command";
		ui->labelResult->setText(resultstr);
		return;
	}
	strcpy(cmdstr, qstr.toLocal8Bit().constData());

	ui->labelResult->setText(cmdstr);
	res = system(cmdstr);
	if (res == 0)
		resultstr = "SUCCESS";
	else if (res == 1)
        resultstr = "FAILED: Wrong number of arguments";
	else if (res == 2)
		resultstr = "FAILED: Read error on input file";
    else if (res == 3)
        resultstr = "FAILED: Sphere too close to image boundary";
    else if (res == 4)
		resultstr = "FAILED: Write error on output file";
	ui->labelResult->setText(resultstr);
}


void MainWindow::makeRanges(void)
{
	int x0, x1, x2, y0, y1, y2, z0, z1, z2, xw, yw, zw, xr, yr, zr;

    x0str = ui->lineEdit_xcentre->text();
    y0str = ui->lineEdit_ycentre->text();
    z0str = ui->lineEdit_zcentre->text();
    xsizestr = ui->lineEdit_xsize->text();
    ysizestr = ui->lineEdit_ysize->text();
    zsizestr = ui->lineEdit_zsize->text();
    x0 = x0str.toInt();
    y0 = y0str.toInt();
    z0 = z0str.toInt();
	xw = xsizestr.toInt();
	yw = ysizestr.toInt();
	zw = zsizestr.toInt();
	xr = xw/2;
	yr = yw/2;
	zr = zw/2;
    x1 = x0 - xr;
    x2 = x0 + xr;
    y1 = y0 - yr;
    y2 = y0 + yr;
    z1 = z0 - zr;
    z2 = z0 + zr;
    x1str = QString::number(x1);
    y1str = QString::number(y1);
    z1str = QString::number(z1);
    x2str = QString::number(x2);
    y2str = QString::number(y2);
    z2str = QString::number(z2);
    ui->lineEdit_x1->setText(x1str);
    ui->lineEdit_y1->setText(y1str);
    ui->lineEdit_z1->setText(z1str);
    ui->lineEdit_x2->setText(x2str);
    ui->lineEdit_y2->setText(y2str);
    ui->lineEdit_z2->setText(z2str);
    ui->checkBox_x->setChecked(false);
    ui->checkBox_y->setChecked(false);
    ui->checkBox_z->setChecked(false);
}

void MainWindow::on_radioButton_centre_toggled(bool checked)
{
    if (checked) {
        ui->groupBox_centre->setEnabled(true);
        ui->groupBox_range->setEnabled(false);
    } else {
        ui->groupBox_centre->setEnabled(false);
        ui->groupBox_range->setEnabled(true);
    }
}
