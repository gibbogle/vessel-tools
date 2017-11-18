#-------------------------------------------------
#
# Project created by QtCreator 2011-04-18T11:11:25
#
#-------------------------------------------------

greaterThan(QT_MAJOR_VERSION, 4) {
    QT += widgets
    DEFINES += HAVE_QT5
}
QT       += core gui

TARGET = am_smoother
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp

HEADERS  += mainwindow.h

FORMS    += mainwindow.ui

DESTDIR = ../../bin
