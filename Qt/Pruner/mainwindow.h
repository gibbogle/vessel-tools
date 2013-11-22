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
    void on_checkBox_ratio_toggled(bool checked);
	void inputFileSelecter();
	void outputFileSelecter();
	void pruner();

private:
    Ui::MainWindow *ui;

public:
	QString inputFileName;
	QString outputFileName;
//	QString outputBaseName;
};

#endif // MAINWINDOW_H
