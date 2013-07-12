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
        QString infoFile = QCoreApplication::applicationDirPath() + "/info/vdistancer_info.txt";
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
    ui->labelResult->setText("");
    ui->labelResultImage->setText("");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::inFileSelecter()
{
    ui->labelInFile->setText("");
	infileName = QFileDialog::getOpenFileName(this,
                tr("Input file"), ".", tr("Amira Files (*.am)"));
	ui->labelInFile->setText(infileName);
    ui->labelResult->setText("");
    ui->labelResultImage->setText("");
}

void MainWindow::outFileSelecter()
{
    ui->labelOutFile->setText("");
    outfileName = QFileDialog::getSaveFileName(this,
                tr("Output file"), ".", tr("Text Files (*.out)"));
    ui->labelOutFile->setText(outfileName);
    ui->labelResult->setText("");
    ui->labelResultImage->setText("");
}

void MainWindow::tiffFileSelecter()
{
    ui->labelTiffFile->setText("");
    tifffileName = QFileDialog::getSaveFileName(this,
                tr("Tiff file"), ".", tr("Tiff Files (*.tif)"));
    ui->labelTiffFile->setText(tifffileName);
    ui->labelResult->setText("");
    ui->labelResultImage->setText("");
}

void MainWindow::randomOption()
{
    if (ui->checkBoxRandom->isChecked()) {
        ui->lineEditGrid_dx->setEnabled(false);
    } else {
        ui->lineEditGrid_dx->setEnabled(true);
    }
    ui->labelResult->setText("");
    ui->labelResultImage->setText("");
}

void MainWindow::sphereOption()
{
    if (ui->checkBoxSphere->isChecked()) {
        ui->lineEditX0->setEnabled(true);
        ui->lineEditY0->setEnabled(true);
        ui->lineEditZ0->setEnabled(true);
        ui->lineEditRadius->setEnabled(true);
    } else {
        ui->lineEditX0->setEnabled(false);
        ui->lineEditY0->setEnabled(false);
        ui->lineEditZ0->setEnabled(false);
        ui->lineEditRadius->setEnabled(false);
    }
    ui->labelResult->setText("");
    ui->labelResultImage->setText("");
}

void MainWindow::imageOption()
{
    if (ui->checkBoxImage->isChecked()) {
        ui->lineEditThreshold->setEnabled(true);
        ui->pushButtonTiffFile->setEnabled(true);
    } else {
        ui->lineEditThreshold->setEnabled(false);
        ui->pushButtonTiffFile->setEnabled(false);
    }
    ui->labelResult->setText("");
    ui->labelResultImage->setText("");
}

void MainWindow::distancer()
{
	int res;
	QString qstr, resultstr;
	char cmdstr[256];
    QString tempfileName = "imagedata.bin";

    qstr = QCoreApplication::applicationDirPath() + "/exec/vdistance ";
	qstr += infileName;
	qstr += " ";
    qstr += outfileName;
    qstr += " ";
    if (!ui->checkBoxRandom->isChecked()) {
        qstr += ui->lineEditGrid_dx->text();
    } else {
        qstr += "0";
    }
    qstr += " ";
    qstr += ui->lineEdit_ncpu->text();
    qstr += " 0";
    qstr += " ";

    if (ui->checkBoxImage->isChecked()) {
        qstr += ui->lineEditThreshold->text();
        qstr += " ";
        qstr += tempfileName;
    } else {
        qstr += "0 dummy";
    }

    if (ui->checkBoxSphere->isChecked()) {
        qstr += " ";
        qstr += ui->lineEditX0->text();
        qstr += " ";
        qstr += ui->lineEditY0->text();
        qstr += " ";
        qstr += ui->lineEditZ0->text();
        qstr += " ";
        qstr += ui->lineEditRadius->text();
    }

    if (qstr.size()>(int)sizeof(cmdstr)-1) {
		printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
		resultstr = "FAILED: cmdstr not big enough for the command";
		ui->labelResult->setText(resultstr);
		return;
	}
	strcpy(cmdstr, qstr.toAscii().constData());

	res = system(cmdstr);
	if (res == 0)
		resultstr = "SUCCESS";
	else if (res == 1)
        resultstr = "FAILED: get_command";
    else if (res == 2)
        resultstr = "FAILED: program name";
    else if (res == 3)
        resultstr = "FAILED: wrong number of arguments";
    else if (res == 4)
        resultstr = "FAILED: allocate error: increase grid spacing";
    else if (res == 5)
        resultstr = "FAILED: indx out of range";
    else
        resultstr = "WTF?";

    ui->labelResult->setText(resultstr);

    if (ui->checkBoxImage->isChecked()) {
        qstr = QCoreApplication::applicationDirPath() + "/exec/maketiff ";
        qstr += tempfileName;
        qstr += " ";
        qstr += tifffileName;
        if (qstr.size()>(int)sizeof(cmdstr)-1) {
            printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
            resultstr = "FAILED: cmdstr not big enough for the command";
            ui->labelResultImage->setText(resultstr);
            return;
        }
        strcpy(cmdstr, qstr.toAscii().constData());

        res = system(cmdstr);
        if (res == 0)
            resultstr = "Image SUCCESS";
        else if (res == 1)
            resultstr = "Image FAILED: wrong number of arguments";
        else if (res == 2)
            resultstr = "Image FAILED: reading image data file";
        else if (res == 3)
            resultstr = "Image FAILED: writing tiff file";
        ui->labelResultImage->setText(resultstr);
    }
}


