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
        QString infoFile = QCoreApplication::applicationDirPath() + "/info/translater_info.txt";
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
                tr("Output file name"), ".", tr("Output Files (*.out)"));
	ui->labelOutputFile->setText(outputFileName);
}

void MainWindow::on_radioButton_len_limit_toggled(bool checked)
{
    ui->lineEdit_len_limit->setEnabled(checked);
}

void MainWindow::on_radioButton_len_diam_limit_toggled(bool checked)
{
    ui->lineEdit_len_diam_limit->setEnabled(checked);
}

void MainWindow::translater()
{
	int res;
    QString limitmodestr, limitvaluestr, qstr, resultstr;
    char cmdstr[512];

    if (ui->radioButton_no_limit->isChecked()) {
        limitmodestr = "0";
        limitvaluestr = "0";
    } else if (ui->radioButton_len_limit->isChecked()) {
        limitmodestr = "1";
        limitvaluestr = ui->lineEdit_len_limit->text();
    } else {
        limitmodestr = "2";
        limitvaluestr = ui->lineEdit_len_diam_limit->text();
    }
    qstr = QCoreApplication::applicationDirPath() + "/exec/translate ";
	qstr += inputFileName;
	qstr += " ";
	qstr += outputFileName;
	qstr += " ";
    qstr += limitmodestr;
    qstr += " ";
    qstr += limitvaluestr;
    qstr += " ";
    qstr += ui->lineEdit_ddiam->text();
    qstr += " ";
    qstr += ui->lineEdit_dlen->text();
    qstr += " ";
    if (ui->checkBoxCreateCmgui->isChecked())
        qstr += "1";
    else
        qstr += "0";
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
        resultstr = "FAILED: Write error on Amira file";
    else if (res == 4)
        resultstr = "FAILED: Write error on CMGUI files";
    else if (res == 5)
        resultstr = "FAILED: EdgeDimensions error";
    else if (res == 6)
        resultstr = "FAILED: CreateDistributions error";
    else
        resultstr = "WTF?";
    ui->labelResult->setText(resultstr);
}


