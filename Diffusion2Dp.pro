TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
#CONFIG -= qt

LIBS += -lquadmath

SOURCES += main.cpp \
    node2d.cpp \
    imc_dfm.cpp \
    diffusion2dp.cpp \
    functor2d.cpp \
    boundary2d.cpp \
    data.cpp

HEADERS += \
    node2d.h \
    imc_dfm.h \
    diffusion2dp.h \
    functor2d.h \
    boundary2d.h \
    data.h

FORMS +=

