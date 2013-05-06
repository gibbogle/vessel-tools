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
        QString infoFile = QCoreApplication::applicationDirPath() + "/info/distancer_info.txt";
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
    outfileName = QFileDialog::getSaveFileName(this,
                tr("Output file"), ".", tr("Text Files (*.out)"));
    ui->labelOutFile->setText(outfileName);
}

void MainWindow::sphereOption()
{
    if (ui->checkBoxSphere->isChecked()) {
        ui->lineEditX0->setEnabled(true);
        ui->lineEditY0->setEnabled(true);
        ui->lineEditZ0->setEnabled(true);
        ui->lineEditRadius->setEnabled(true);
    } else {
        ui->lineEditX0->setEnabled(false);
        ui->lineEditY0->setEnabled(false);
        ui->lineEditZ0->setEnabled(false);
        ui->lineEditRadius->setEnabled(false);
    }
}

void MainWindow::distancer()
{
	int res;
	QString qstr, resultstr;
	char cmdstr[256];

    qstr = QCoreApplication::applicationDirPath() + "/exec/proximity ";
	qstr += infileName;
	qstr += " ";
    qstr += outfileName;
    qstr += " ";
    qstr += ui->lineEditGrid_dx->text();

    if (ui->checkBoxSphere->isChecked()) {
        qstr += " ";
        qstr += ui->lineEditX0->text();
        qstr += " ";
        qstr += ui->lineEditY0->text();
        qstr += " ";
        qstr += ui->lineEditZ0->text();
        qstr += " ";
        qstr += ui->lineEditRadius->text();
    }

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
        resultstr = "FAILED: get_command";
    else if (res == 2)
        resultstr = "FAILED: program name";
    else if (res == 3)
        resultstr = "FAILED: wrong number of arguments";
    else if (res == 4)
        resultstr = "FAILED: allocate error: increase grid spacing";
    else if (res == 5)
        resultstr = "FAILED: indx out of range";
    else
        resultstr = "WTF?";

    ui->labelResult->setText(resultstr);
}


