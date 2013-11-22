#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QMainWindow>

#include "imagewidget.h"

namespace Ui {
    class ImageViewer;
}

class QAction;
class QMenu;
class QMenuBar;
class QScrollArea;
class QScrollBar;
class QVBoxLayout;

class ImageViewer : public QMainWindow {
    Q_OBJECT

public:

    ImageViewer(QWidget* parent = 0);
    ~ImageViewer();

protected:
    //    void closeEvent(QCloseEvent *event);


private slots:

    /**
     * Load and display an image from file
     */
    void open1();
    void open2();

    /**
     * Save as, to the disk, displayed image 
     */
    void saveAs();

    /**
     * Apply median filter to image
     */
//    void medianFilter();

    /**
     * apply an Anisotropic diffusion filter to image
     */
//    void gradientAnisotropicDiffusionFilter();

    /**
     * 
     */
    void about();

    void nextFrame();
    void prevFrame();
    void updateImages();
    void subtracter();
    void subtractImage();
    void imageSelecter();
    void saveImageFile();
    bool imageEquivalence();
    void on_pushButtonAddImage_clicked();

private:

    Ui::ImageViewer *ui;

    void opener(int);
    void createActions();
    void createMenus();
    void createStatusBar();

    ImageWidget *imageWidget1;
    ImageWidget *imageWidget2;

    QAction *openAct1;
    QAction *openAct2;
    QAction *saveAsAct;
    QAction *exitAct;
    QAction *aboutAct;
    QAction *aboutQtAct;

    QAction *medianFilterAct;
    QAction *GADFilterAct;

    QMenu *fileMenu;
    QMenu *filterMenu;
    QMenu *helpMenu;
    
    QLabel *statusLabel;

	bool useITK;
    int iframe;
};

#endif
