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
	QString infoFile = QCoreApplication::applicationDirPath() + "/info/pruner_info.txt";
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
		tr("Input object file"), ".", tr("SpatialGraph Files (*.am)"));
	ui->labelInputFile->setText(inputFileName);
}

void MainWindow::outputFileSelecter()
{
	ui->labelResult->setText("");
	outputFileName = QFileDialog::getSaveFileName(this,
		tr("Output file base name"), ".", tr("SpatialGraph Files (*.am)"));
	ui->labelOutputFile->setText(outputFileName);
}

void MainWindow::on_checkBox_ratio_toggled(bool checked)
{
    if (checked) {
        ui->labelRatioLimit->setText("Twig length/diameter ratio limit");
    } else {
        ui->labelRatioLimit->setText("Twig length limit (um)");
    }
}


void MainWindow::pruner()
{
	int res;
	QString qstr, resultstr;
    char cmdstr[1024];

	qstr = QCoreApplication::applicationDirPath() + "/exec/prune ";
	qstr += inputFileName;
	qstr += " ";
	qstr += outputFileName;
	qstr += " ";
    if (ui->checkBox_ratio->isChecked())
        qstr += "1 ";
    else
        qstr += "0 ";
    qstr += ui->lineEditRatioLimit->text();
	qstr += " ";
    qstr += ui->lineNprunecycles->text();
    qstr += " ";
    if (ui->checkBoxPrune->isChecked())
        qstr += "1 ";
	else
        qstr += "0 ";
	if (ui->checkBoxCreateCmgui->isChecked())
        qstr += "1 ";
	else
        qstr += "0 ";
    qstr += ui->lineEdit_ddiam->text();
    qstr += " ";
    qstr += ui->lineEdit_dlen->text();
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
        resultstr = "FAILED: error reading input am file";
	else if (res == 3)
        resultstr = "FAILED: error in adjoinEdges";
    else if (res == 4)
        resultstr = "FAILED: error in deloop";
    else if (res == 5)
        resultstr = "FAILED: error in pruner";
    else if (res == 6)
        resultstr = "FAILED: error in checkEdgeEndPts";
    else if (res == 7)
        resultstr = "FAILED: error in squeezer";
    else if (res == 8)
        resultstr = "FAILED: error in WriteAmiraFile";
    else if (res == 9)
        resultstr = "FAILED: error in CreateDistributions";
    else if (res == 10)
        resultstr = "FAILED: error in WriteCmguiData";
    else
        resultstr = "FAILED: unknown error";
    ui->labelResult->setText(resultstr);
}


