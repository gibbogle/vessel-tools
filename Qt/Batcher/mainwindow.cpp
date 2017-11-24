#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    headers << "Site name" << "x0" << "y0" << "z0" << "R";
    ui->tableWidget->setHorizontalHeaderLabels(headers);
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
    ui->textEdit->setReadOnly(true);
    QString infoFile = QCoreApplication::applicationDirPath() + "/info/batcher_info.txt";
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

void MainWindow:: populateTable() {
    ui->tableWidget->clear();
    QFile sitefile(ui->labelSiteFile->text());
    if (!sitefile.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    while (!sitefile.atEnd()) {
         QString aline = sitefile.readLine();
         processLine(aline);
    }
    ui->tableWidget->setHorizontalHeaderLabels(headers);
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
    outputBaseFileName = QFileDialog::getSaveFileName(this,
        tr("Output base name"), ".", tr("Output base name (*.)"));
    ui->labelOutputFile->setText(outputBaseFileName);
}

void MainWindow::on_radioButton_len_limit_toggled(bool checked)
{
    ui->lineEdit_len_limit->setEnabled(checked);
}

void MainWindow::on_radioButton_len_diam_limit_toggled(bool checked)
{
    ui->lineEdit_len_diam_limit->setEnabled(checked);
}

void MainWindow::on_checkBoxImage_toggled(bool checked)
{
    ui->lineEditThreshold->setEnabled(checked);
}

int MainWindow::blocker(QString numstr)
{
    QString qstr, resultstr;
    char cmdstr[2048];
    int res;

    qstr = QCoreApplication::applicationDirPath() + "/exec/am_block ";
    qstr += inputFileName;
    qstr += " ";
    qstr += "subregion_" + numstr + ".am";

    qstr += " ";
    qstr += x0_str;
    qstr += " ";
    qstr += y0_str;
    qstr += " ";
    qstr += z0_str;
    qstr += " ";
    qstr += R_str;
    qstr += " 0 0";
    if (qstr.size()>(int)sizeof(cmdstr)-1) {
        printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
        resultstr = "FAILED: cmdstr not big enough for the command";
        ui->labelResult->setText(resultstr);
        return -1;
    }

    ui->labelCommand->setText(qstr);
    strcpy(cmdstr, qstr.toLocal8Bit().constData());
    res = system(cmdstr);

    if (res == 0) {
        resultstr = "SUCCESS";
    } else {
        resultstr = "am_block failed to crop for site: " + numstr;
        ui->labelResult->setText(resultstr);
    }
    return res;
}

int MainWindow::distancer(QString numstr)
{
    int res;
    QString qstr, resultstr;
    char cmdstr[2048];
    tempfileName = "imagedata.bin";

    qstr = QCoreApplication::applicationDirPath() + "/exec/vdistancef ";
    qstr += inputFileName;
    qstr += " ";
    qstr += outputBaseFileName + "_dist_" + numstr + ".out";
    qstr += " ";
    qstr += ui->lineEditGrid_dx->text();
    qstr += " ";
    qstr += "4";    //ui->lineEdit_ncpu->text();
    qstr += " ";
//    qstr += "0";    // no need for random grid points
//    qstr += " ";

    if (ui->checkBoxImage->isChecked()) {
        qstr += ui->lineEditThreshold->text();
        qstr += " ";
        qstr += tempfileName;
    } else {
        qstr += "0 dummy";
    }
    qstr += " ";
    qstr += ui->lineEditVx->text();
    qstr += " ";
    qstr += ui->lineEditVy->text();
    qstr += " ";
    qstr += ui->lineEditVz->text();

//    if (ui->checkBoxSubregion->isChecked()) {
//        if (ui->radioButtonCube->isChecked()) {
//            qstr += " C ";
//        } else {
            qstr += " S ";
//        }
        qstr += x0_str; //ui->lineEditX0->text();
        qstr += " ";
        qstr += y0_str; //ui->lineEditY0->text();
        qstr += " ";
        qstr += z0_str; //ui->lineEditZ0->text();
        qstr += " ";
        qstr += R_str;  //ui->lineEditRadius->text();
//    }
    // NOT always use close file
//    if (ui->checkBox_use_close->isChecked()) {
//        qstr += " ";
//        qstr += closefileName;
//    }

    ui->labelCommand->setText(qstr);
    if (qstr.size()>(int)sizeof(cmdstr)-1) {
        printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
        resultstr = "FAILED: cmdstr not big enough for the command";
        ui->labelResult->setText(resultstr);
        return -1;
    }
    strcpy(cmdstr, qstr.toLocal8Bit().constData());

    res = system(cmdstr);
    if (res == 0)
        resultstr = "SUCCESS";
    else if (res == 1)
        resultstr = "FAILED: Command line error (1)";
    else if (res == 2)
        resultstr = "FAILED: Command line argument count error (2)";
    else if (res == 3)
        resultstr = "FAILED: Command line argument error (3)";
    else if (res == 4)
        resultstr = "FAILED: Allocate error: increase grid spacing";
    else if (res == 5)
        resultstr = "FAILED: indx out of range";
    else if (res == 6)
        resultstr = "FAILED: Open failure on close file";
    else if (res == 7)
        resultstr = "FAILED: Supplied close file has incorrect format";
    else if (res == 8)
        resultstr = "FAILED: Access error in closedata: kbytes > nmbytes";
    else if (res == 9)
        resultstr = "FAILED: Exceeded dimension of array imagedata";
    else
        resultstr = "WTF?";

    ui->labelResult->setText(resultstr);
    if (res != 0) return res;
    if (ui->checkBoxImage->isChecked()) {
        res = tiffer(numstr);
    }
//    if (res == 0) {
//        tiff_ready = true;
//    }
//    if (ui->checkBoxImage->isChecked() && is_tiff && tiff_ready){
//        ui->pushButtonMakeTiff->setEnabled(true);
//    }
    return res;
}

int MainWindow::tiffer(QString numstr) {
    int res;
    QString qstr, resultstr;
    char cmdstr[2048];

    qstr = QCoreApplication::applicationDirPath() + "/exec/maketiff ";
    qstr += tempfileName;
    qstr += " ";
    qstr += outputBaseFileName + "_" + numstr + ".tif";   // tifffileName;
    qstr += " 1";
    if (qstr.size()>(int)sizeof(cmdstr)-1) {
        printf("Failed to convert qstr->cmdstr since qstr didn't fit\n");
        resultstr = "FAILED: cmdstr not big enough for the command";
        ui->labelResultImage->setText(resultstr);
        return -1;
    }
    strcpy(cmdstr, qstr.toLocal8Bit().constData());

    res = system(cmdstr);
    if (res == 0)
        resultstr = "Image SUCCESS";
    else if (res == 1)
        resultstr = "Image FAILED: wrong number of arguments";
    else if (res == 2)
        resultstr = "Image FAILED: reading image data file";
    else if (res == 3)
        resultstr = "Image FAILED: writing tiff file";
    ui->labelResultImage->setText(resultstr);
    return res;
}

void MainWindow::batch_analyser()
{
	int res;
    char cmdstr[2048];
    QString limitmodestr, limitvaluestr, qstr, resultstr;
    QString outfile_str;

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

    // Two sweeps through the subregions.
    // First sweep:
    // To compute Cm at the specified start location, using Cm sphere start radius, with the full network
    // Second sweep:
    // Crop out a spherical subregion of specified radius, to compute network statistics for the subregion network
    // (Note: do we want to use the connected network?)

    for (int irun=0; irun<2; irun++) {
        int isite = 0;
        for (;;) {
            isite++;
            QTableWidgetItem* item = ui->tableWidget->item(isite-1,0);
            if (!item || item->text().isEmpty()) break;
            x0_str = ui->tableWidget->item(isite-1,1)->text();
            y0_str = ui->tableWidget->item(isite-1,2)->text();
            z0_str = ui->tableWidget->item(isite-1,3)->text();
            R_str = ui->tableWidget->item(isite-1,4)->text();
            QString numstr = QString("%1").arg(isite, 2, 10, QChar('0'));
            if (irun == 0) {
                outfile_str = outputBaseFileName + "_";
                if (ui->checkBoxDistance->isChecked()) {
                    int res = distancer(numstr);
                    if (res != 0) return;
                }
            } else {
                outfile_str = outputBaseFileName + "_sub_";

                // Crop the full network to give a spherical subregion network
                subregionFileName = "subregion_" + numstr + ".am";
                int res = blocker(numstr);
                if (res != 0) return;
            }
            outfile_str += numstr;
            outfile_str += ".out";
            qstr = QCoreApplication::applicationDirPath() + "/exec/conduit_analyse ";
            if (irun == 0)
                qstr += inputFileName;
            else
                qstr += subregionFileName;
            qstr += " ";
            qstr += outfile_str;
            qstr += " ";
            qstr += ui->lineEdit_sfactor->text();
            qstr += " ";
            qstr += ui->lineEdit_pow->text();
            qstr += " ";
            if (irun == 0)
                qstr += ui->lineEdit_ntrials->text();
            else
                qstr += "100";
            qstr += " ";
            qstr += x0_str;
            qstr += " ";
            qstr += y0_str;
            qstr += " ";
            qstr += z0_str;
            qstr += " ";
            if (irun == 0)
                qstr += ui->lineEdit_Cm_radius->text();
            else
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

            ui->labelCommand->setText(qstr);
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
                resultstr = "FAILED: read error on input file";
            else if (res == 3)
                resultstr = "FAILED: write error on output file";
            else if (res == 4)
                resultstr = "FAILED: error in CreateDistributions";
            else if (res == 5)
                resultstr = "FAILED: MAXBLOCK exceeded";
            ui->labelResult->setText(resultstr);
            if (res != 0) return;
        }
    }
}

