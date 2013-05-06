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
    void inverter();

private:
    Ui::MainWindow *ui;

public:
	QString inputFileName;
    QString outputFileName;

private slots:

};

#endif // MAINWINDOW_H
