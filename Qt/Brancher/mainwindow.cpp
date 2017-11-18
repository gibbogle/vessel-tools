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
        QString infoFile = QCoreApplication::applicationDirPath() + "/info/brancher_info.txt";
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
    infileName = QFileDialog::getOpenFileName(this,tr("Input file"), ".", tr("Amira Files (*.am)"));
	ui->labelInFile->setText(infileName);
}

void MainWindow::outFileSelecter()
{
	ui->labelResult->setText("");
	QFileDialog dialog(this);
	dialog.setFileMode(QFileDialog::AnyFile);
    outfileName = dialog.getSaveFileName(this,tr("Output file"), ".", tr("Output Files (*.out)"));
	ui->labelOutFile->setText(outfileName);
}


void MainWindow::brancher()
{
	int res;
	QString qstr, resultstr;
    char cmdstr[2048];

    qstr = QCoreApplication::applicationDirPath() + "/exec/short ";
	qstr += infileName;
	qstr += " ";
	qstr += outfileName;
	qstr += " ";
    if (ui->lineEditArteryVertex->text().compare("") == 0)
        qstr += "0";
    else
        qstr += ui->lineEditArteryVertex->text();
    qstr += " ";
    qstr += ui->lineEdit_aR->text();
    qstr += " ";
    qstr += ui->lineEdit_aG->text();
    qstr += " ";
    qstr += ui->lineEdit_aB->text();
    qstr += " ";
    if (ui->lineEditVeinVertex->text().compare("") == 0)
        qstr += "0";
    else
     qstr += ui->lineEditVeinVertex->text();
    qstr += " ";
    qstr += ui->lineEdit_vR->text();
    qstr += " ";
    qstr += ui->lineEdit_vG->text();
    qstr += " ";
    qstr += ui->lineEdit_vB->text();
    qstr += " ";
    if (ui->checkBox_use_branch_count->isChecked()) {
        qstr += "1 ";
        qstr += ui->lineEditRatio->text();
    } else {
        qstr += "0 0";
    }

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
		resultstr = "FAILED: wrong number of arguments";
	else if (res == 2)
        resultstr = "FAILED: Read error on input AM file";
    else if (res == 3)
        resultstr = "FAILED: CreateVertexLinks error";
    else if (res == 4)
        resultstr = "FAILED: Quantify error";
    else if (res == 5)
        resultstr = "FAILED: Write error on output CMGUI files";
    else
        resultstr = "WTF?";

    ui->labelResult->setText(resultstr);
}


