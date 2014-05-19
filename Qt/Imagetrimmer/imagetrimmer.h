#ifndef IMAGETRIMMER_H
#define IMAGETRIMMER_H

#include <QMainWindow>
#include <QFileDialog>

#include "imagewidget.h"

namespace Ui {
    class ImageTrimmer;
}

class QAction;
class QMenu;
class QMenuBar;
class QScrollArea;
class QScrollBar;
class QVBoxLayout;

class ImageTrimmer : public QMainWindow {
    Q_OBJECT

public:

    ImageTrimmer(QWidget* parent = 0);
    ~ImageTrimmer();

protected:

private slots:

    void about();

    void open1();
    void nextFrame();
    void prevFrame();
    void updateImages();
    void saveLine();
    void saveLineMatrix();
    void loadLineMatrix();

public slots:
    void setLineSaveStatus(bool);

private:

    Ui::ImageTrimmer *ui;

    void opener(int);
    void createActions();
    void createMenus();
    void createStatusBar();

    ImageWidget *imageWidget1;

    QAction *openAct1;
    QAction *exitAct;
    QAction *aboutAct;
    QAction *aboutQtAct;

    QMenu *fileMenu;
    QMenu *filterMenu;
    QMenu *helpMenu;
    
    QLabel *statusLabel;

    QString lineFileName;
	bool useITK;
    int iframe;
};

#endif
