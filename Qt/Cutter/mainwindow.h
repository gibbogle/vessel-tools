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
    void cutter();

private:
    Ui::MainWindow *ui;

public:
	QString inputFileName;

private slots:

};

#endif // MAINWINDOW_H
