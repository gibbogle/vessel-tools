#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::inFileSelecter()
{
	ui->labelResult->setText("");
	infileName = QFileDialog::getOpenFileName(this,
		tr("Input file"), ".", tr("Image Files (*.tif *.img)"));
	ui->labelInFile->setText(infileName);
}

void MainWindow::outFileSelecter()
{
	ui->labelResult->setText("");
	QFileDialog dialog(this);
	dialog.setFileMode(QFileDialog::AnyFile);
	outfileName = dialog.getSaveFileName(this,
						tr("Output file"), ".", tr("TIFF Files (*.tif)"));
//	QString outfileName = QFileDialog::getOpenFileName(this,
//		tr("Output file"), "/users/gib/ln_structure", tr("TIFF Files (*.tif)"));
	ui->labelOutFile->setText(outfileName);
}

void MainWindow::compresser()
{
	int res;
	QString qstr, resultstr;
	char cmdstr[256];

	qstr = QCoreApplication::applicationDirPath() + "/exec/compress ";
	qstr += infileName;
	qstr += " ";
	qstr += outfileName;
	if (qstr.size()>(int)sizeof(cmdstr)-1) {
		printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
		resultstr = "FAILED: cmdstr not big enough for the command";
		ui->labelResult->setText(resultstr);
		return;
	}
	strcpy(cmdstr, qstr.toLocal8Bit().constData());

	res = system(cmdstr);
	if (res == 0)
		resultstr = "SUCCESS";
	else if (res == 1)
		resultstr = "FAILED: input and/or output file not specified";
	else if (res == 2)
		resultstr = "FAILED: Read error on input file";
	else if (res == 3)
		resultstr = "FAILED: Write error on output file";
	ui->labelResult->setText(resultstr);
}

void MainWindow::uncompresser()
{
	int res;
	QString qstr, resultstr;
	char cmdstr[256];

	qstr = QCoreApplication::applicationDirPath() + "/exec/uncompress ";
	qstr += infileName;
	qstr += " ";
	qstr += outfileName;
	if (qstr.size()>(int)sizeof(cmdstr)-1) {
		printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
		resultstr = "FAILED: cmdstr not big enough for the command";
		ui->labelResult->setText(resultstr);
		return;
	}
	strcpy(cmdstr, qstr.toLocal8Bit().constData());

	res = system(cmdstr);
	if (res == 0)
		resultstr = "SUCCESS";
	else if (res == 1)
		resultstr = "FAILED: input and/or output file not specified";
	else if (res == 2)
		resultstr = "FAILED: Read error on input file";
	else if (res == 3)
		resultstr = "FAILED: Write error on output file";
	ui->labelResult->setText(resultstr);
}
