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
        QString infoFile = QCoreApplication::applicationDirPath() + "/info/am_joiner_info.txt";
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

void MainWindow::outFileSelecter()
{
	ui->labelResult->setText("");
	QFileDialog dialog(this);
	dialog.setFileMode(QFileDialog::AnyFile);
        outfileName = dialog.getSaveFileName(this,tr("Amira file"), ".", tr("Amira Files (*.am)"));
	ui->labelOutFile->setText(outfileName);
}

void MainWindow::am_joiner()
{
	int res;
	QString qstr, resultstr;
	char cmdstr[256];

    qstr = QCoreApplication::applicationDirPath() + "/exec/am_join ";
	qstr += outfileName;
    qstr += " ";
    qstr += ui->lineEditBasename->text();
    qstr += " ";
    if (ui->checkBox_x->isChecked())
        qstr += "1 ";
    else
        qstr += "0 ";
    if (ui->checkBox_y->isChecked())
        qstr += "1 ";
    else
        qstr += "0 ";
    if (ui->checkBox_z->isChecked())
        qstr += "1 ";
    else
        qstr += "0 ";
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
            resultstr = "FAILED: AddNetwork error";
    else if (res == 4)
            resultstr = "FAILED: Write error on output AM file";
    else if (res == 5)
            resultstr = "FAILED: Write error on output CMGUI files";
    ui->labelResult->setText(resultstr);
}


