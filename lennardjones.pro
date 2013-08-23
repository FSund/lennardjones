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
    src/generator.cpp

HEADERS += \
    src/mainapplication.h \
    src/lib.h \
    src/inlines.h \
    src/atom.h \
    src/box.h \
    src/state.h \
    src/linked_list.h \
    src/integrator.h \
    src/generator.h

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
