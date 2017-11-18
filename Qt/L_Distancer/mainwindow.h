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
    void cboxBdryFileToggled();
    void inAmiraFileSelecter();
    void lymphaticFileSelecter();
    void bdryFileSelecter();
    void outFileSelecter();
    void L_distancer();
//    void sphereOption();

private:
    Ui::MainWindow *ui;
    void inBdryFileSelecter();
    void outBdryFileSelecter();

public:
    QString inAmiraFileName;
    QString lymphaticFileName;
    QString bdryFileName;
    QString inBdryFileName;
    QString outBdryFileName;
    QString outFileName;


};

#endif // MAINWINDOW_H
