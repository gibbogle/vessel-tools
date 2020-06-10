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
    QString infoFile = QCoreApplication::applicationDirPath() + "/info/conduit_analyser_info.txt";
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
    connect(ui->buttonGroup_mode, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(mode_radioButtonChanged(QAbstractButton*)));
    mode_radioButtonChanged(ui->radioButton_clean);
    connect(ui->buttonGroup_jump, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(jump_radioButtonChanged(QAbstractButton*)));
    jumpstr = "0";
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::inputFileSelecter()
{
	ui->labelResult->setText("");
	inputFileName = QFileDialog::getOpenFileName(this,
        tr("Input Amira file"), ".", tr("Amira Files (*.am)"));
	ui->labelInputFile->setText(inputFileName);
}

void MainWindow::outputFileSelecter()
{
	ui->labelResult->setText("");
	outputFileName = QFileDialog::getSaveFileName(this,
        tr("Output file"), ".", tr("Output Files (*.out)"));
	ui->labelOutputFile->setText(outputFileName);
}

//void MainWindow::join_checkbox(){
//    if (ui->checkBox_join->isChecked()) {
//        ui->groupBox_traverse->setDisabled(true);
//        ui->groupBox_join->setEnabled(true);
//        ui->lineEdit_sfactor->setDisabled(true);
//    } else {
//        ui->groupBox_traverse->setEnabled(true);
//        ui->groupBox_join->setDisabled(true);
//        ui->lineEdit_sfactor->setEnabled(true);
//    }
//}

void MainWindow::savepaths_checkbox(){
    if (ui->checkBox_savepaths->isChecked()) {
        ui->lineEdit_npaths->setEnabled(true);
    } else {
        ui->lineEdit_npaths->setDisabled(true);
    }
}

void MainWindow::on_radioButton_len_limit_toggled(bool checked)
{
    ui->lineEdit_len_limit->setEnabled(checked);
}

void MainWindow::on_radioButton_len_diam_limit_toggled(bool checked)
{
    ui->lineEdit_len_diam_limit->setEnabled(checked);
}

//------------------------------------------------------------------------------------------------------
void MainWindow::mode_radioButtonChanged(QAbstractButton *b)
{
//    QString wtag = b->objectName();
//    printf("radioButtonChanged: %s\n",wtag.toLocal8Bit().constData());
    if (ui->radioButton_clean->isChecked()) {
        mode = 1;
        ui->label_CM->setEnabled(false);
        ui->groupBox_statistics->setEnabled(false);
        ui->label_statistics->setEnabled(false);
        ui->groupBox_traverse->setEnabled(false);
        ui->label_sfactor->setEnabled(false);
        ui->lineEdit_sfactor->setEnabled(false);
    } else if (ui->radioButton_connect->isChecked()) {
        mode = 2;
        ui->label_CM->setEnabled(false);
        ui->groupBox_statistics->setEnabled(true);
        ui->label_statistics->setEnabled(true);
        ui->groupBox_traverse->setEnabled(false);
        ui->label_sfactor->setEnabled(true);
        ui->lineEdit_sfactor->setEnabled(true);
    } else if (ui->radioButton_cm->isChecked()) {
        mode = 3;
        ui->label_CM->setEnabled(true);
        ui->groupBox_statistics->setEnabled(true);
        ui->label_statistics->setEnabled(true);
        ui->groupBox_traverse->setEnabled(true);
        ui->label_sfactor->setEnabled(true);
        ui->lineEdit_sfactor->setEnabled(true);
    }
}

//------------------------------------------------------------------------------------------------------
void MainWindow::jump_radioButtonChanged(QAbstractButton *b)
{
    if (ui->radioButton_nojump->isChecked()) {
        jumpstr = "0";
    } else if (ui->radioButton_ejump->isChecked()) {
        jumpstr = "1";
    } else if (ui->radioButton_vjump->isChecked()) {
        jumpstr = "2";
    }
}

void MainWindow::conduit_analyser()
{
	int res;
    char cmdstr[2048];
    QString limitmodestr, limitvaluestr, qstr, resultstr;

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

    qstr = QCoreApplication::applicationDirPath() + "/exec/conduit_analyse ";
	qstr += inputFileName;
	qstr += " ";
	qstr += outputFileName;

    if (mode == 2) {
        qstr += " ";
        qstr += ui->lineEdit_sfactor->text();
        qstr += " ";
        qstr += ui->lineEdit_maxlen->text();
        qstr += " ";
    } else if (mode == 3) {
        qstr += " ";
        qstr += ui->lineEdit_sfactor->text();
        qstr += " ";
        qstr += ui->lineEdit_pow->text();
        qstr += " ";
        qstr += ui->lineEdit_ntrials->text();
        qstr += " ";
        qstr += ui->lineEdit_x0->text();
        qstr += " ";
        qstr += ui->lineEdit_y0->text();
        qstr += " ";
        qstr += ui->lineEdit_z0->text();
        qstr += " ";
        qstr += ui->lineEdit_radius->text();
        qstr += " ";
        qstr += ui->lineEdit_deadend_radius->text();
        qstr += " ";
        qstr += ui->lineEdit_speed->text();
        qstr += " ";
        qstr += ui->lineEdit_CV->text();
        qstr += " ";
        qstr += ui->lineEdit_npaths->text();
        qstr += " ";
        qstr += jumpstr;
        qstr += " ";
        qstr += ui->lineEdit_jumpprob->text();
    }
    if (mode > 1) {
        qstr += " ";
        qstr += limitmodestr;
        qstr += " ";
        qstr += limitvaluestr;
        qstr += " ";
        qstr += ui->lineEdit_ddiam->text();
        qstr += " ";
        qstr += ui->lineEdit_dlen->text();
    }

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
        resultstr = "FAILED: read error on input file";
	else if (res == 3)
        resultstr = "FAILED: write error on output file";
	else if (res == 4)
        resultstr = "FAILED: error in CreateDistributions";
    else if (res == 5)
        resultstr = "FAILED: MAXBLOCK exceeded";
    else if (res == 6)
        resultstr = "FAILED: MAX_SFIBRES exceeded: reduce radius";
    ui->labelResult->setText(resultstr);
}

