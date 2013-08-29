TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -std=c++11

ROOT_DIR = $$PWD
SRC_DIR = $$PWD/src
INCLUDEPATH += $$ROOT_DIR

SOURCES += main.cpp \
    src/mainapplication.cpp \
    src/lib.cpp \
    src/atom.cpp \
    src/box.cpp \
    src/state.cpp \
    src/integrator.cpp \
    src/generator.cpp \
    src/processor.cpp

HEADERS += \
    src/mainapplication.h \
    src/lib.h \
    src/inlines.h \
    src/atom.h \
    src/box.h \
    src/state.h \
    src/linked_list.h \
    src/integrator.h \
    src/generator.h \
    src/processor.h

LIBS += -larmadillo

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_LFLAGS -= -O1
    QMAKE_LFLAGS += -O3
    QMAKE_LFLAGS_RELEASE -= -O1
    QMAKE_LFLAGS_RELEASE += -O3
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O3
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}


QMAKE_LFLAGS -= -O1
QMAKE_CXXFLAGS -= -O2
QMAKE_LFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc
QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
