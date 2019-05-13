 CONFIG      += debug_and_release
#  CONFIG      += release

 QT          += opengl

 HEADERS     = glwidget.h \
               window.h \
               Mouse.h \
               Grid.h \
               Common.h \
               pcgsolver/pcg_solver.h \

 SOURCES     = glwidget.cpp \
               glwidget_levelset.cpp \
               main.cpp \
               window.cpp \
               Mouse.cpp \
               Grid.cpp \
               Common.cpp

QMAKE_CXXFLAGS += -std=c++0x

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3 -funroll-loops -ffast-math -fomit-frame-pointer
QMAKE_LFLAGS_RELEASE += -O3

QMAKE_CXXFLAGS_DEBUG += -pg -g -O2

LIBS += -lblas

 # install
 #target.path = $$[QT_INSTALL_EXAMPLES]/opengl/2dpainting
 #sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS flip2d.pro
 #sources.path = $$[QT_INSTALL_EXAMPLES]/opengl/2dpainting
 #INSTALLS += target sources

 symbian: include($$QT_SOURCE_TREE/examples/symbianpkgrules.pri)
