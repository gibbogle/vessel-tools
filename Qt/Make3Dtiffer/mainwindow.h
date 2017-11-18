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
    void directorySelecter();
	void outputFileSelecter();
    void make3Dtiffer();

private:
    Ui::MainWindow *ui;

public:
    QString directory;
	QString inputFileName;
	QString outputFileName;

};

#endif // MAINWINDOW_H
