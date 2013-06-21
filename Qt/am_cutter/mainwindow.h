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
    void am_cutter();

private:
    Ui::MainWindow *ui;

public:
	QString infileName;

};

#endif // MAINWINDOW_H
