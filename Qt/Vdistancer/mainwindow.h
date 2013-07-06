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
	void inFileSelecter();
    void outFileSelecter();
    void distancer();
    void sphereOption();
    void randomOption();

private:
    Ui::MainWindow *ui;

public:
	QString infileName;
    QString outfileName;


};

#endif // MAINWINDOW_H
