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
	void objectFileSelecter();
	void skelFileSelecter();
	void outputFileSelecter();
	void diamCheckBox();
    void junctionRadioButton();
	void topology();

private:
    Ui::MainWindow *ui;

public:
	QString objectFileName;
	QString skelFileName;
	QString outputFileName;
	QString outputBaseName;
};

#endif // MAINWINDOW_H
