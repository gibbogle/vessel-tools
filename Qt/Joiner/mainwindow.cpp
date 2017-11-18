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
    QString infoFile = QCoreApplication::applicationDirPath() + "/info/cutter_info.txt";
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

void MainWindow::inputLEFileSelecter()
{
	ui->labelResult->setText("");
    inputLEFileName = QFileDialog::getOpenFileName(this,
		tr("Input TIFF file"), ".", tr("TIFF Files (*.tif)"));
    ui->labelInputLEFile->setText(inputLEFileName);
}

void MainWindow::inputGTFileSelecter()
{
    ui->labelResult->setText("");
    inputGTFileName = QFileDialog::getOpenFileName(this,
        tr("Input TIFF file"), ".", tr("TIFF Files (*.tif)"));
    ui->labelInputGTFile->setText(inputGTFileName);
}

void MainWindow::outputFileSelecter()
{
    ui->labelResult->setText("");
    outputFileName = QFileDialog::getSaveFileName(this,
        tr("Input TIFF file"), ".", tr("TIFF Files (*.tif)"));
    ui->labelOutputFile->setText(outputFileName);
}

void MainWindow::joiner()
{
	int res;
    QString qstr, axis, resultstr;
	char cmdstr[2048];

    if (ui->radioButton_X->isChecked()) {
        axis = " x ";
    } else if (ui->radioButton_Y->isChecked()) {
        axis = " y ";
    } else if (ui->radioButton_Z->isChecked()) {
        axis = " z ";
    }

    qstr = QCoreApplication::applicationDirPath() + "/exec/join ";
    qstr += inputLEFileName;
    qstr += " ";
    qstr += inputGTFileName;
    qstr += " ";
    qstr += outputFileName;
    qstr += axis;
    if (ui->checkBoxCompress->isChecked())
        qstr += " C";
	else
        qstr += " U";

	if (qstr.size()>(int)sizeof(cmdstr)-1) {
		printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
		resultstr = "FAILED: cmdstr not big enough for the command";
		ui->labelResult->setText(resultstr);
		return;
	}
	strcpy(cmdstr, qstr.toLocal8Bit().constData());

//	ui->label_cmd->setText(cmdstr);
	res = system(cmdstr);
	if (res == 0)
		resultstr = "SUCCESS";
	else if (res == 1)
		resultstr = "FAILED: wrong number of arguments";
	else if (res == 2)
        resultstr = "FAILED: Bad axis choice";
	else if (res == 3)
        resultstr = "FAILED: Read error on input file";
	else if (res == 4)
        resultstr = "FAILED: Bad cut position";
    else if (res == 3)
        resultstr = "FAILED: Write error on output file";
    ui->labelResult->setText(resultstr);
}
