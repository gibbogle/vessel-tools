#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
	void inputFileSelecter();
	void outputFileSelecter();
	void checkBox_x();
	void checkBox_y();
	void checkBox_z();
	void cropper();

private:
    Ui::MainWindow *ui;
    void makeRanges();
    QString x0str, x1str, x2str, y0str, y1str, y2str, z0str, z1str, z2str;
    QString xsizestr, ysizestr, zsizestr;

public:
	QString inputFileName;
	QString outputFileName;

private slots:
    void on_radioButton_centre_toggled(bool checked);
    void on_radioButton_block_toggled(bool checked);
};

#endif // MAINWINDOW_H
