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
	QString infoFile = QCoreApplication::applicationDirPath() + "/info/connecter_info.txt";
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

void MainWindow::inFileSelecter()
{
	ui->labelResult->setText("");
	infileName = QFileDialog::getOpenFileName(this,
		tr("Input file"), ".", tr("TIFF Files (*.tif)"));
	ui->labelInFile->setText(infileName);
}

void MainWindow::connecter()
{
	int res;
	QString qstr, resultstr;
	char cmdstr[256];

	qstr = QCoreApplication::applicationDirPath() + "/exec/connect ";
	qstr += infileName;
	qstr += " ";
	qstr += ui->lineEditBaseName->text();
	qstr += " ";
	qstr += ui->lineEditMinVoxels->text();
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
        resultstr = "FAILED: Wrong number of arguments";
	else if (res == 2)
		resultstr = "FAILED: Read error on input file";
	else if (res == 3)
		resultstr = "FAILED: Write error on output file";
	else if (res == 4)
        resultstr = "FAILED: Out of memory: allocation of lablist";
    else if (res == 5)
        resultstr = "FAILED: Out of memory: allocation of label failed";
    else if (res == 6)
        resultstr = "FAILED: Size of list[] in tracer (MAXLABELS) exceeded";
    else if (res == 7)
        resultstr = "FAILED: Too many objects: MAXOBJECTS exceeded";
    else if (res == 8)
        resultstr = "FAILED: Too many pairs: MAXLABLIST exceeded";
    else if (res == 9)
        resultstr = "FAILED: MAXLABELS exceeded in connecter";
    else if (res == 10)
        resultstr = "FAILED: Max short size (max16) exceeded";
    else if (res == 11)
        resultstr = "FAILED: Size of explist (MAXEXP) exceeded";
    ui->labelResult->setText(resultstr);
}


