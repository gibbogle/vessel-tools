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
        QString infoFile = QCoreApplication::applicationDirPath() + "/info/reducer_info.txt";
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
    infileName = QFileDialog::getOpenFileName(this, tr("Input file"), ".", tr("Amira Files (*.am)"));
	ui->labelInFile->setText(infileName);
}

void MainWindow::outFileSelecter()
{
	ui->labelResult->setText("");
    outfileName = QFileDialog::getSaveFileName(this, tr("Output file"), ".", tr("Amira Files (*.am)"));
    ui->labelOutFile->setText(outfileName);
}

void MainWindow::reducer()
{
	int res;
	QString qstr, resultstr;
    char cmdstr[2048];

    qstr = QCoreApplication::applicationDirPath() + "/exec/reduce ";
	qstr += infileName;
	qstr += " ";
	qstr += outfileName;
	qstr += " ";
    qstr += ui->lineEdit_len_diam_limit->text();
    qstr += " ";
    qstr += ui->lineEdit_len_limit->text();
    qstr += " ";
    if (ui->checkBoxCmgui->isChecked())             // cmgui_flag
        qstr += "1 ";
    else
        qstr += "0 ";

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
        resultstr = "FAILED: Read error on input AM file";
    else if (res == 3)
        resultstr = "FAILED: ReduceNetwork error";
    else if (res == 4)
        resultstr = "FAILED: Write error on output AM file";
    else if (res == 5)
        resultstr = "FAILED: Write error on output CMGUI files";
    else if (res == 6)
        resultstr = "FAILED: Error generating distributions";
    else
        resultstr = "WTF?";

    ui->labelResult->setText(resultstr);
}


