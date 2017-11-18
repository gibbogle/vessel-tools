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

TARGET    = quantifier
TEMPLATE  = app


SOURCES  += main.cpp mainwindow.cpp quantify.cpp imageviewer.cpp

HEADERS  += mainwindow.h quantify.h network.h imageviewer.h

FORMS    += mainwindow.ui

DEFINES  += _CRT_SECURE_NO_DEPRECATE

DESTDIR   = ../../bin
