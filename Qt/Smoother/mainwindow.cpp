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
	QString infoFile = QCoreApplication::applicationDirPath() + "/info/smoother_info.txt";
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
    infileName = QFileDialog::getOpenFileName(this, tr("Input file"), ".", tr("TIFF Files (*.tif)"));
	ui->labelInFile->setText(infileName);
}

void MainWindow::outFileSelecter()
{
	ui->labelResult->setText("");
//	QFileDialog dialog(this);
//	dialog.setFileMode(QFileDialog::AnyFile);
//	outfileName = dialog.getSaveFileName(this,tr("Output file"), ".", tr("TIFF Files (*.tif)"));
    outfileName = QFileDialog::getSaveFileName(this, tr("Output file"), ".", tr("TIFF Files (*.tif)"));
    ui->labelOutFile->setText(outfileName);
}

void MainWindow::smoother()
{
	int res;
	QString qstr, resultstr;
	char cmdstr[256];

	qstr = QCoreApplication::applicationDirPath() + "/exec/smooth ";
	qstr += infileName;
	qstr += " ";
	qstr += outfileName;
	qstr += " ";
	qstr += ui->lineEditSize->text();
	qstr += " ";
	qstr += ui->lineEditNCPU->text();
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
		resultstr = "FAILED: Read error on input file";
	else if (res == 3)
		resultstr = "FAILED: Write error on output file";
	ui->labelResult->setText(resultstr);
}


