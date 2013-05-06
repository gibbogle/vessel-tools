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

void MainWindow::cboxBdryFileToggled()
{
    if (ui->checkBoxBdryFile->isChecked()) {
        ui->pushButtonBdryFile->setText("Old lymphatic bdry file");
        ui->pushButtonLymphaticFile->setDisabled(true);
    } else {
        ui->pushButtonBdryFile->setText("New lymphatic bdry file");
        ui->pushButtonLymphaticFile->setEnabled(true);
    }
}

void MainWindow::inAmiraFileSelecter()
{
	ui->labelResult->setText("");
    inAmiraFileName = QFileDialog::getOpenFileName(this,
                tr("Input network file"), ".", tr("Amira Files (*.am)"));
    ui->labelInAmiraFile->setText(inAmiraFileName);
}

void MainWindow::lymphaticFileSelecter()
{
    ui->labelResult->setText("");
    lymphaticFileName = QFileDialog::getOpenFileName(this,
                tr("Input lymphatic file"), ".", tr("Tiff Files (*.tif)"));
    ui->labelLymphaticFile->setText(lymphaticFileName);
}

void MainWindow::outFileSelecter()
{
    ui->labelResult->setText("");
    outFileName = QFileDialog::getSaveFileName(this,
                tr("Output distribution file"), ".", tr("Text Files (*.out)"));
    ui->labelOutFile->setText(outFileName);
}

void MainWindow::bdryFileSelecter()
{
    if (ui->checkBoxBdryFile->isChecked()) {
        inBdryFileSelecter();
        ui->labelBdryFile->setText(inBdryFileName);
        bdryFileName = inBdryFileName;
    } else {
        outBdryFileSelecter();
        ui->labelBdryFile->setText(outBdryFileName);
        bdryFileName = outBdryFileName;
    }
}

void MainWindow::inBdryFileSelecter()
{
    ui->labelResult->setText("");
    inBdryFileName = QFileDialog::getOpenFileName(this,
                tr("Input lymphatic boundary file"), ".", tr("Text Files (*.dat)"));
}

void MainWindow::outBdryFileSelecter()
{
    ui->labelResult->setText("");
    outBdryFileName = QFileDialog::getSaveFileName(this,
                tr("Output lymphatic boundary file"), ".", tr("Text Files (*.dat)"));
}

/*
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
*/
void MainWindow::L_distancer()
{
	int res;
	QString qstr, resultstr;
    char cmdstr[2048];

/*
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
*/

// First create the lymphatic boundary data file (if needed)
    if (!ui->checkBoxBdryFile->isChecked()) {
        qstr = QCoreApplication::applicationDirPath() + "/exec/boundary ";
        qstr += lymphaticFileName;
        qstr += " ";
        qstr += bdryFileName;

        if (qstr.size()>(int)sizeof(cmdstr)-1) {
            printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
            resultstr = "FAILED: cmdstr not big enough for the boundary command";
            ui->labelResult->setText(resultstr);
            return;
        }
        strcpy(cmdstr, qstr.toAscii().constData());

        res = system(cmdstr);
        if (res == 0)
            resultstr = "Boundary data extracted, computing distribution ...";
        else if (res == 1)
            resultstr = "FAILED: boundary: get_command";
        else if (res == 2)
            resultstr = "FAILED: boundary: program name";
        else if (res == 3)
            resultstr = "FAILED: boundary: wrong number of arguments";
        else if (res == 4)
            resultstr = "FAILED: boundary: allocate error: increase grid spacing";
        else if (res == 5)
            resultstr = "FAILED: boundary: indx out of range";
        else
            resultstr = " WTF?  boundary";

        ui->labelResult->setText(resultstr);
        if (res != 0) return;
    }
    qstr = QCoreApplication::applicationDirPath() + "/exec/lymphatic ";
    qstr += inAmiraFileName;
    qstr += " ";
    qstr += bdryFileName;
    qstr += " ";
    qstr += outFileName;

    if (qstr.size()>(int)sizeof(cmdstr)-1) {
        printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
        resultstr = "FAILED: cmdstr not big enough for the lymphatic command";
        ui->labelResult->setText(resultstr);
        return;
    }
    strcpy(cmdstr, qstr.toAscii().constData());

    res = system(cmdstr);
    if (res == 0)
        resultstr = "SUCCESS: lymphatic";
    else if (res == 1)
        resultstr = "FAILED: lymphatic: get_command";
    else if (res == 2)
        resultstr = "FAILED: lymphatic: program name";
    else if (res == 3)
        resultstr = "FAILED: lymphatic: wrong number of arguments";
    else if (res == 4)
        resultstr = "FAILED: lymphatic: allocate error: increase grid spacing";
    else if (res == 5)
        resultstr = "FAILED: lymphatic: indx out of range";
    else
        resultstr = "WTF?  lymphatic";

    ui->labelResult->setText(resultstr);
}


