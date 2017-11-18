#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#ifdef HAVE_QT5
#include <QtWidgets/QMainWindow>
#else
#include <QtGui/QMainWindow>
#endif

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
	void inFileSelecter();
	void outFileSelecter();
    void selecter();

private:
    Ui::MainWindow *ui;
    void makeRanges();
    QString x0str, x1str, x2str, y0str, y1str, y2str, z0str, z1str, z2str;
    QString xsizestr, ysizestr, zsizestr;

public:
	QString infileName;
	QString outfileName;

private slots:
    void on_radioButton_centre_toggled(bool checked);

};

#endif // MAINWINDOW_H
