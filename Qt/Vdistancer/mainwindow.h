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
    void tiffer();
    void imageOption();
    void imageChange();

private:
    Ui::MainWindow *ui;

public:
	QString infileName;
    QString closefileName;
    QString outfileName;
    QString tifffileName;
    bool is_tiff;
    bool tiff_ready;
    QString tempfileName;
private slots:
    void on_checkBox_use_close_toggled(bool checked);
    void on_checkBoxSubregion_toggled(bool checked);
};

#endif // MAINWINDOW_H
