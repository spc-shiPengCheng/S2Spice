QT -= gui

TEMPLATE = lib
DEFINES += VECTORFITTING_LIBRARY

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    vectorfitting.cpp

HEADERS += \
    VectorFitting_global.h \
    vectorfitting.h

TRANSLATIONS += \
    VectorFitting_zh_CN.ts

INCLUDEPATH += ../eigen-5.0.0/eigen-5.0.0

# Default rules for deployment.
unix {
    target.path = /usr/lib
}
!isEmpty(target.path): INSTALLS += target
