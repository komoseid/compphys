TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
INCLUDEPATH += /usr/local/include
LIBS+=  -L/usr/local/lib -lblas -llapack -larmadillo
SOURCES += main.cpp \
    jacobi.cpp

HEADERS += \
    jacobi.h
