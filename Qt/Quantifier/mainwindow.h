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
    void amFileSelecter();
    void closeFileSelecter();
    void resultFileSelecter();
    void computeVolume();
    void computeArea();
    void computeVessels();
    int checkSlice(int, int);
    void voxelChanged();
    void checkReady();
    void doSetup();

private:
    Ui::MainWindow *ui;

public:
    float voxelsize[3];
    int nvoxels[3];
    QString amFileName;
    QString closeFileName;
    QString resultFileName;
    bool is_amfile;
    bool is_closefile;
    bool is_resultfile;
    bool is_ready;
};

#endif // MAINWINDOW_H
