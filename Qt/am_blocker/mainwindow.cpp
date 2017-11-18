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
    connect(ui->radioButton_block_range, SIGNAL(toggled(bool)), this, SLOT(onToggled()));
    connect(ui->radioButton_block_centre, SIGNAL(toggled(bool)), this, SLOT(onToggled()));
    connect(ui->radioButton_sphere, SIGNAL(toggled(bool)), this, SLOT(onToggled()));

	ui->textEdit->setReadOnly(true);
        QString infoFile = QCoreApplication::applicationDirPath() + "/info/am_blocker_info.txt";
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

void MainWindow::makeRanges(void)
{
    int x0, x1, x2, y0, y1, y2, z0, z1, z2, xw, yw, zw, xr, yr, zr;

    x0str = ui->lineEdit_xcentre->text();
    y0str = ui->lineEdit_ycentre->text();
    z0str = ui->lineEdit_zcentre->text();
    xsizestr = ui->lineEdit_xsize->text();
    ysizestr = ui->lineEdit_ysize->text();
    zsizestr = ui->lineEdit_zsize->text();
    x0 = x0str.toInt();
    y0 = y0str.toInt();
    z0 = z0str.toInt();
    xw = xsizestr.toInt();
    yw = ysizestr.toInt();
    zw = zsizestr.toInt();
    xr = xw/2;
    yr = yw/2;
    zr = zw/2;
    x1 = x0 - xr;
    x2 = x0 + xr;
    y1 = y0 - yr;
    y2 = y0 + yr;
    z1 = z0 - zr;
    z2 = z0 + zr;
    x1str = QString::number(x1);
    y1str = QString::number(y1);
    z1str = QString::number(z1);
    x2str = QString::number(x2);
    y2str = QString::number(y2);
    z2str = QString::number(z2);
    ui->lineEdit_x1->setText(x1str);
    ui->lineEdit_y1->setText(y1str);
    ui->lineEdit_z1->setText(z1str);
    ui->lineEdit_x2->setText(x2str);
    ui->lineEdit_y2->setText(y2str);
    ui->lineEdit_z2->setText(z2str);
}

//void MainWindow::on_radioButton_centre_toggled(bool checked)
//{
//    if (checked) {
//        ui->groupBox_centre->setEnabled(true);
//        ui->groupBox_range->setEnabled(false);
//        ui->groupBox_sphere->setEnabled(false);
//    } else {
//        ui->groupBox_centre->setEnabled(false);
//        ui->groupBox_range->setEnabled(true);
//    }
//}

void MainWindow::onToggled()
{
    if (ui->radioButton_block_centre->isChecked()) {
        ui->groupBox_centre->setEnabled(true);
        ui->groupBox_range->setEnabled(false);
        ui->groupBox_sphere->setEnabled(false);
    }
    if (ui->radioButton_block_range->isChecked()) {
        ui->groupBox_centre->setEnabled(false);
        ui->groupBox_range->setEnabled(true);
        ui->groupBox_sphere->setEnabled(false);
    }
    if (ui->radioButton_sphere->isChecked()) {
        ui->groupBox_centre->setEnabled(false);
        ui->groupBox_range->setEnabled(false);
        ui->groupBox_sphere->setEnabled(true);
    }
}

void MainWindow::blocker()
{
	int res;
	QString qstr, resultstr;
    char cmdstr[2048];

    if (ui->radioButton_block_range->isChecked()) {
        x1str = ui->lineEdit_x1->text();
        x2str = ui->lineEdit_x2->text();
        y1str = ui->lineEdit_y1->text();
        y2str = ui->lineEdit_y2->text();
        z1str = ui->lineEdit_z1->text();
        z2str = ui->lineEdit_z2->text();
    } else if (ui->radioButton_block_centre->isChecked()) {
        makeRanges();
    } else {
        x0str = ui->lineEdit_x0->text();
        y0str = ui->lineEdit_y0->text();
        z0str = ui->lineEdit_z0->text();
        Rstr = ui->lineEdit_radius->text();
    }
    qstr = QCoreApplication::applicationDirPath() + "/exec/am_block ";
    qstr += infileName;
    qstr += " ";
    qstr += outfileName;
    qstr += " ";
    if (ui->radioButton_sphere->isChecked()) {
        qstr += x0str;
        qstr += " ";
        qstr += y0str;
        qstr += " ";
        qstr += z0str;
        qstr += " ";
        qstr += Rstr;
        qstr += " ";
    } else {
        qstr += x1str;
        qstr += " ";
        qstr += x2str;
        qstr += " ";
        qstr += y1str;
        qstr += " ";
        qstr += y2str;
        qstr += " ";
        qstr += z1str;
        qstr += " ";
        qstr += z2str;
        qstr += " ";
    }
    if (ui->checkBoxConnect->isChecked())           // connect_flag
        qstr += "1 ";
    else
        qstr += "0 ";
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
	strcpy(cmdstr, qstr.toLocal8Bit().constData());

	res = system(cmdstr);
	if (res == 0)
		resultstr = "SUCCESS";
	else if (res == 1)
		resultstr = "FAILED: wrong number of arguments";
	else if (res == 2)
        resultstr = "FAILED: Read error on input AM file";
    else if (res == 3)
        resultstr = "FAILED: CreateSelectNet error";
    else if (res == 4)
        resultstr = "FAILED: Write error on output AM file";
    else if (res == 5)
        resultstr = "FAILED: Write error on output CMGUI files";
    else
        resultstr = "WTF?";

    ui->labelResult->setText(resultstr);
}


