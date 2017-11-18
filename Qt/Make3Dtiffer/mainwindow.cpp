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
    QString infoFile = QCoreApplication::applicationDirPath() + "/info/make3Dtiffer_info.txt";
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

void MainWindow::directorySelecter()
{
    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::Directory);
	ui->labelResult->setText("");
    directory = dialog.getExistingDirectory();
    ui->labelDirectory->setText(directory);
}

void MainWindow::outputFileSelecter()
{
	ui->labelResult->setText("");
	outputFileName = QFileDialog::getSaveFileName(this,
		tr("Output TIFF file"), ".", tr("TIFF Files (*.tif)"));
	ui->labelOutputFile->setText(outputFileName);
}

void MainWindow::make3Dtiffer()
{
	int res;
    QString qstr, resultstr, templatestr, bitstr;
    char cmdstr[2048];

    templatestr = directory + "/" + ui->lineEditTemplate->text();
    if (ui->radioButton_8bits->isChecked())
        bitstr = " 8 ";
    else
        bitstr = " 16 ";
    qstr = QCoreApplication::applicationDirPath() + "/exec/make3Dtiff ";
    qstr += templatestr;
    qstr += bitstr;
    qstr += ui->lineEdit_firstnum->text();
    qstr += " ";
    qstr += ui->lineEdit_lastnum->text();
    if (ui->checkBoxCompress->isChecked())
        qstr += " 1 ";
    else
        qstr += " 0 ";
    qstr += outputFileName;

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
		resultstr = "FAILED: wrong number of arguments";
	else if (res == 2)
		resultstr = "FAILED: Read error on input file";
	else if (res == 3)
		resultstr = "FAILED: Write error on output file";
    else
        resultstr = "FAILED: some other error";
    ui->labelResult->setText(resultstr);
//    ui->labelResult->setText(qstr);
}

