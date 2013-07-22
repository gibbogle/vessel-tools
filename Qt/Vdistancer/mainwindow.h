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
    void closeFileSelecter();
    void outFileSelecter();
    void tiffFileSelecter();
    void distancer();
    void sphereOption();
//    void randomOption();
    void imageOption();

private:
    Ui::MainWindow *ui;

public:
	QString infileName;
    QString closefileName;
    QString outfileName;
    QString tifffileName;


};

#endif // MAINWINDOW_H
