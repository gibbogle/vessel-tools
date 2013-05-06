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
	QString infoFile = QCoreApplication::applicationDirPath() + "/info/topology_info.txt";
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

void MainWindow::objectFileSelecter()
{
	ui->labelResult->setText("");
	objectFileName = QFileDialog::getOpenFileName(this,
		tr("Input object file"), ".", tr("TIFF Files (*.tif)"));
	ui->labelObjectFile->setText(objectFileName);
}

void MainWindow::skelFileSelecter()
{
	ui->labelResult->setText("");
	skelFileName = QFileDialog::getOpenFileName(this,
		tr("Input skeleton file"), ".", tr("TIFF Files (*.tif)"));
	ui->labelSkelFile->setText(skelFileName);
}

void MainWindow::outputFileSelecter()
{
	ui->labelResult->setText("");
	outputFileName = QFileDialog::getSaveFileName(this,
		tr("Output file base name"), ".", tr("Output Files (*.out)"));
	ui->labelOutputFile->setText(outputFileName);
}

void MainWindow::diamCheckBox()
{
	if (ui->checkBoxDiam->isChecked()) {
		ui->lineEditDiam->setEnabled(true);
	} else {
		ui->lineEditDiam->setDisabled(true);
	}
}

void MainWindow::topology()
{
	int res;
	QString qstr, resultstr;
        char cmdstr[1024];

	// Strip ".out" off outputFileName
//	QFileInfo fi(outputFileName);
//	outputBaseName = fi.absolutePath() + "/" + fi.baseName();
	qstr = QCoreApplication::applicationDirPath() + "/exec/topo ";
	qstr += skelFileName;
	qstr += " ";
	qstr += objectFileName;
	qstr += " ";
//	qstr += outputBaseName;
	qstr += outputFileName;
	qstr += " ";
    qstr += ui->lineEditVoxelSize_xy->text();
	qstr += " ";
    qstr += ui->lineEditVoxelSize_z->text();
    qstr += " ";
    qstr += ui->lineEditPrune->text();
	qstr += " ";
	if (ui->checkBoxDiam->isChecked()) {
		qstr += ui->lineEditDiam->text();
	} else {
		qstr += "0";
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
		resultstr = "FAILED: wrong number of arguments";
	else if (res == 2)
		resultstr = "FAILED: read error on skeleton file";
	else if (res == 3)
		resultstr = "FAILED: read error on segmented vessel file";
	else if (res == 4)
		resultstr = "FAILED: skeleton and vessel files differ in size";
	else if (res == 5)
		resultstr = "FAILED: error in TraceSkeleton()";
	else if (res == 6)
		resultstr = "FAILED: error in GetDiameters()";
	else if (res == 7)
		resultstr = "FAILED: error in simplify()";
	else if (res == 8)
		resultstr = "FAILED: error in squeezer()";
	else if (res == 9)
		resultstr = "FAILED: error in CreateDistributions()";
	else if (res == 10)
		resultstr = "FAILED: error in WriteCmguiData()";
	else if (res == 11)
		resultstr = "FAILED: error in WriteAmiraFile()";
	else
		resultstr = "FAILED: some other error";
	ui->labelResult->setText(resultstr);
}


