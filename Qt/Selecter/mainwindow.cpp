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
    connect(ui->buttonGroup_selection, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(on_radioButton_diameter_changed()));

	ui->textEdit->setReadOnly(true);
        QString infoFile = QCoreApplication::applicationDirPath() + "/info/selecter_info.txt";
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
//	QFileDialog dialog(this);
//	dialog.setFileMode(QFileDialog::AnyFile);
//        outfileName = dialog.getSaveFileName(this,tr("Output file"), ".", tr("Amira Files (*.am)"));
    outfileName = QFileDialog::getSaveFileName(this, tr("Output file"), ".", tr("Amira Files (*.am)"));
    ui->labelOutFile->setText(outfileName);
}

void MainWindow::on_radioButton_diameter_changed()
{
    bool is_diam = ui->radioButton_diameter->isChecked();
    ui->groupBox_diameter->setEnabled(is_diam);
    ui->groupBox_criterion->setEnabled(is_diam);
    ui->groupBox_length->setEnabled(!is_diam);
    ui->checkBoxConnect->setEnabled(is_diam);
    if (!is_diam) {
        ui->checkBoxConnect->setChecked(false);
    }
}

void MainWindow::selecter()
{
	int res;
	QString qstr, resultstr;
    char cmdstr[2048];

    qstr = QCoreApplication::applicationDirPath() + "/exec/select ";
	qstr += infileName;
	qstr += " ";
	qstr += outfileName;
	qstr += " ";
    if (ui->radioButton_diameter->isChecked()) {      // mode_flag
        qstr += "1 ";
        qstr += ui->lineEditDmin->text();
        qstr += " ";
        qstr += ui->lineEditDmax->text();
        qstr += " ";
    } else {
        qstr += "0 ";
        qstr += ui->lineEditLmin->text();
        qstr += " ";
        qstr += ui->lineEditLmax->text();
        qstr += " ";
    }
    if (ui->checkBoxConnect->isChecked())           // connect_flag
        qstr += "1 ";
    else
        qstr += "0 ";
    if (ui->radioButton_average->isChecked())       // diam_flag
        qstr += "0 ";
    else if (ui->radioButton_majority->isChecked())
        qstr += "1 ";
    else
        qstr += "2 ";
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
        resultstr = "FAILED: CheckNetwork error";
    else if (res == 4)
        resultstr = "FAILED: CreateDiamSelectNet error";
    else if (res == 5)
        resultstr = "FAILED: CreateLenSelectNet error";
    else if (res == 6)
        resultstr = "FAILED: CreateLargestConnectedNet error";
    else if (res == 7)
        resultstr = "FAILED: Write error on output AM file";
    else if (res == 8)
        resultstr = "FAILED: Write error on output CMGUI files";
    else if (res == 9)
        resultstr = "FAILED: Error in computing network statistics";
    else
        resultstr = "WTF?";

    ui->labelResult->setText(resultstr);
}


