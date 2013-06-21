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
	void outFileSelecter();
    void am_joiner();

private:
    Ui::MainWindow *ui;

public:
	QString outfileName;


};

#endif // MAINWINDOW_H
