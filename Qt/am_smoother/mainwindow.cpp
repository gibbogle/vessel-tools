#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
	ui->textEdit->setReadOnly(true);
        QString infoFile = QCoreApplication::applicationDirPath() + "/info/am_smoother_info.txt";
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
                tr("Input file"), ".", tr("Amira Files (*.am)"));
	ui->labelInFile->setText(infileName);
}

void MainWindow::outFileSelecter()
{
	ui->labelResult->setText("");
	QFileDialog dialog(this);
	dialog.setFileMode(QFileDialog::AnyFile);
        outfileName = dialog.getSaveFileName(this,tr("Amira file"), ".", tr("Amira Files (*.am)"));
	ui->labelOutFile->setText(outfileName);
}

void MainWindow::am_smoother()
{
	int res;
	QString qstr, resultstr;
	char cmdstr[256];

    double fraction = ui->lineEditFraction->text().toDouble();
    double diameter = ui->lineEditDiameter->text().toDouble();
    if (fraction == 0 && diameter == 0) {
        QMessageBox msgBox;
        msgBox.setText("Since both faction and diameter are zero, no changes will be made to the network");
        msgBox.setInformativeText("Do you wish to continue?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::No);
        int ret = msgBox.exec();
        if (ret == QMessageBox::No) return;
    }
    qstr = QCoreApplication::applicationDirPath() + "/exec/am_smooth ";
	qstr += infileName;
	qstr += " ";
	qstr += outfileName;
	qstr += " ";
    qstr += ui->lineEditFraction->text();
    qstr += " ";
    qstr += ui->lineEditDiameter->text();
    qstr += " ";
    if (ui->checkBoxCmgui->isChecked())
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
            resultstr = "FAILED: Write error on output AM file";
    else if (res == 3)
            resultstr = "FAILED: Write error on output CMGUI files";
    ui->labelResult->setText(resultstr);
}


