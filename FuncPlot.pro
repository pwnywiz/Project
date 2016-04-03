#-------------------------------------------------
#
# Project created by QtCreator 2016-01-20T16:15:01
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = FuncPlot
TEMPLATE = app
LIBS += -lQtGnuplot -llapacke

CONFIG += debug


SOURCES += main.cpp\
        visualization.cpp \
    cpp/compan.cpp \
    cpp/generalized.cpp \
    cpp/initfuncs.cpp \
    cpp/interpolation.cpp \
    cpp/polyfuncs.cpp \
    cpp/polys.cpp \
    cpp/sy.cpp \
    cpp/sylvester.cpp

HEADERS  += visualization.h \
    cpp/compan.h \
    cpp/generalized.h \
    cpp/initfuncs.h \
    cpp/interpolation.h \
    cpp/polyfuncs.h \
    cpp/polys.h \
    cpp/sy.h \
    cpp/sylvester.h

FORMS    += \
    visualization.ui

OTHER_FILES +=
