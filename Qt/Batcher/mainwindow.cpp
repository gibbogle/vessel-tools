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
    // Set fixed width to sum of column widths + 26
    ui->tableWidget->setFixedWidth(286);
    ui->tableWidget->setColumnWidth(0,100);
    ui->tableWidget->setColumnWidth(1,40);
    ui->tableWidget->setColumnWidth(2,40);
    ui->tableWidget->setColumnWidth(3,40);
    ui->tableWidget->setColumnWidth(4,40);
    // Set fixed height to sum of row heights + 26
    ui->tableWidget->setFixedHeight(226);
    for (int i=0; i<10; i++) {
        ui->tableWidget->setRowHeight(i,20);
    }
}

void MainWindow:: populateTable() {
    ui->tableWidget->clear();
    QFile sitefile(ui->labelSiteFile->text());
    if (!sitefile.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    nsites = 0;
    while (!sitefile.atEnd()) {
         QString aline = sitefile.readLine();
         processLine(aline);
    }
}

void MainWindow::processLine(QString aline) {
    static int row = 0;
    QStringList ss = aline.split(QRegExp("\\s+"),QString::SkipEmptyParts);

    if(ui->tableWidget->rowCount() < row + 1)
        ui->tableWidget->setRowCount(row + 1);
    if(ui->tableWidget->columnCount() < ss.size())
        ui->tableWidget->setColumnCount( ss.size() );

    for( int column = 0; column < ss.size(); column++)
    {
        QTableWidgetItem *newItem = new QTableWidgetItem( ss.at(column) );
        ui->tableWidget->setItem(row, column, newItem);
    }

    row++;
    if (ss.at(0).compare("") != 0) nsites = row;
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

void MainWindow::siteFileSelecter()
{
	ui->labelResult->setText("");
    siteFileName = QFileDialog::getOpenFileName(this,
        tr("Site file"), ".", tr("Site Files (*.txt)"));
    ui->labelSiteFile->setText(siteFileName);
    populateTable();
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

//void MainWindow::savepaths_checkbox(){
//    if (ui->checkBox_savepaths->isChecked()) {
//        ui->lineEdit_npaths->setEnabled(true);
//    } else {
//        ui->lineEdit_npaths->setDisabled(true);
//    }
//}

void MainWindow::on_radioButton_len_limit_toggled(bool checked)
{
    ui->lineEdit_len_limit->setEnabled(checked);
}

void MainWindow::on_radioButton_len_diam_limit_toggled(bool checked)
{
    ui->lineEdit_len_diam_limit->setEnabled(checked);
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

    for (int isite = 1; isite<=nsites; isite++) {
    QString name_str = ui->tableWidget->item(isite-1,0)->text();
    QString x0_str = ui->tableWidget->item(isite-1,1)->text();
    QString y0_str = ui->tableWidget->item(isite-1,2)->text();
    QString z0_str = ui->tableWidget->item(isite-1,3)->text();
    QString R_str = ui->tableWidget->item(isite-1,4)->text();
    QString outfile_str = "cube";
//    outfile_str += ".out";
    QString numstr = QString("%1").arg(isite, 2, 10, QChar('0'));
    outfile_str += numstr;
    outfile_str += ".out";
    qstr = QCoreApplication::applicationDirPath() + "/exec/conduit_analyse ";
	qstr += inputFileName;
	qstr += " ";
    qstr += outfile_str;
    qstr += " ";
    qstr += ui->lineEdit_sfactor->text();
    qstr += " ";
    qstr += ui->lineEdit_pow->text();
    qstr += " ";
    qstr += ui->lineEdit_ntrials->text();
    qstr += " ";

    qstr += x0_str;
    qstr += " ";
    qstr += y0_str;
    qstr += " ";
    qstr += z0_str;
    qstr += " ";
    qstr += R_str;

    qstr += " ";
    qstr += ui->lineEdit_speed->text();
    qstr += " ";
    qstr += "0";    //ui->lineEdit_npaths->text();
    qstr += " ";
    qstr += limitmodestr;
    qstr += " ";
    qstr += limitvaluestr;
    qstr += " ";
    qstr += ui->lineEdit_ddiam->text();
    qstr += " ";
    qstr += ui->lineEdit_dlen->text();

/*
    qstr += outputFileName;
    qstr += " ";
    qstr += ui->lineEdit_sfactor->text();
    qstr += " ";
    if (ui->checkBox_join->isChecked()) {
        qstr += ui->lineEdit_maxlen->text();
    } else {
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
        qstr += ui->lineEdit_speed->text();
        qstr += " ";
        qstr += ui->lineEdit_npaths->text();
    }
    qstr += " ";
    qstr += limitmodestr;
    qstr += " ";
    qstr += limitvaluestr;
    qstr += " ";
    qstr += ui->lineEdit_ddiam->text();
    qstr += " ";
    qstr += ui->lineEdit_dlen->text();
*/
	if (qstr.size()>(int)sizeof(cmdstr)-1) {
		printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
		resultstr = "FAILED: cmdstr not big enough for the command";
		ui->labelResult->setText(resultstr);
		return;
	}
    ui->labelResult->setText(qstr);

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
    ui->labelResult->setText(resultstr);
    }
}

