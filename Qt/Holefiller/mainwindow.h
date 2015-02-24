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
    void holefiller();

private:
    Ui::MainWindow *ui;

public:
	QString inputFileName;
	QString outputFileName;

};

#endif // MAINWINDOW_H
