#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#ifdef HAVE_QT5
#include <QtWidgets/QMainWindow>
#else
#include <QtGui/QMainWindow>
#endif

#include <QAbstractButton>

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
    void savepaths_checkbox();
    void on_radioButton_len_limit_toggled(bool checked);
    void on_radioButton_len_diam_limit_toggled(bool checked);
    void conduit_analyser();
    void mode_radioButtonChanged(QAbstractButton *b);
    void jump_radioButtonChanged(QAbstractButton *b);

private:
    Ui::MainWindow *ui;

public:
	QString inputFileName;
	QString outputFileName;
    int mode;
    QString jumpstr;
};

#endif // MAINWINDOW_H
