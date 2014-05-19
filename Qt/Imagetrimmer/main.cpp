#include <QApplication>

#include "imagetrimmer.h"

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    ImageTrimmer imageTrimmer;
    imageTrimmer.show();
    return app.exec();
}
