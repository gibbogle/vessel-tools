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
	void inputFileSelecter();
//	void outputFileSelecter();
    void peeler();

private:
    Ui::MainWindow *ui;

public:
	QString inputFileName;
//	QString outputFileName;
    QString prefixpath;
};

#endif // MAINWINDOW_H
