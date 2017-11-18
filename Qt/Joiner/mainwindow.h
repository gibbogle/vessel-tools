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
    void inputLEFileSelecter();
    void inputGTFileSelecter();
    void outputFileSelecter();
    void joiner();

private:
    Ui::MainWindow *ui;

public:
    QString inputLEFileName;
    QString inputGTFileName;
    QString outputFileName;

private slots:

};

#endif // MAINWINDOW_H
