#include <QApplication>

#include "3Dviewer.h"

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
//    ImageViewer *imageViewer = new ImageViewer;
//    imageViewer->show();
    ImageViewer imageViewer;
    imageViewer.show();
    return app.exec();
}
