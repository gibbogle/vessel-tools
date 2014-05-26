#-------------------------------------------------
#
# Project created by QtCreator 2011-04-18T11:11:25
#
#-------------------------------------------------

QT       += core gui

TARGET    = tiffquantifier
TEMPLATE  = app

INCLUDEPATH += C:/Program Files/ITK/include/ITK-4.0

SOURCES  += main.cpp mainwindow.cpp tiffquantify.cpp imageviewer.cpp

HEADERS  += mainwindow.h tiffquantify.h imageviewer.h

FORMS    += mainwindow.ui

DEFINES  += _CRT_SECURE_NO_DEPRECATE

DESTDIR   = ../../bin
