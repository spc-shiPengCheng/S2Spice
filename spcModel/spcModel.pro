QT       += core gui charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

QT += printsupport

CONFIG += c++17

CONFIG += openmp       # 启用 OpenMP

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += $$PWD/Library/eigen-5.0.0/eigen-5.0.0


SOURCES += \
    $$PWD/Library/VectorFitting/check_sparameter.cpp \
    $$PWD/Library/VectorFitting/spice_netlist_gengerate.cpp \
    $$PWD/Library/VectorFitting/vectorfitting.cpp \
    GUI/mainwindow.cpp \
    GUI/st_generate_netlist_dialog.cpp \
    GUI/st_setting_spice_solver_dialog.cpp \
    GUI/st_show_sparameter_singular_dialog.cpp \
    GUI/st_show_syz_parameter.cpp \
    GUI/st_vectorfitting_parameter_setting.cpp \
    GUI\qcustomplot.cpp\
    Tool/st_convert_syz.cpp \
    Tool/st_file_io.cpp \
    Tool/st_logger.cpp \
    main.cpp

HEADERS += \
    $$PWD/Library/VectorFitting/check_sparameter.h \
    $$PWD/Library/VectorFitting/spice_netlist_gengerate.h \
    $$PWD/Library/VectorFitting/vectorfitting.h \
    GUI/mainwindow.h \
    GUI/st_generate_netlist_dialog.h \
    GUI/st_setting_spice_solver_dialog.h \
    GUI/st_show_sparameter_singular_dialog.h \
    GUI/st_show_syz_parameter.h \
    GUI/st_vectorfitting_parameter_setting.h \
    GUI\qcustomplot.h\
    Tool/st_convert_syz.h \
    Tool/st_file_io.h \
    Tool/st_logger.h \
    Tool/st_struct_file.h

FORMS += \
    GUI/mainwindow.ui \
    GUI/st_generate_netlist_dialog.ui \
    GUI/st_setting_spice_solver_dialog.ui \
    GUI/st_show_sparameter_singular_dialog.ui \
    GUI/st_vectorfitting_parameter_setting.ui \
    GUI/stshowsyzparameter.ui

TRANSLATIONS += \
    S2Spice_zh_CN.ts
CONFIG += lrelease
CONFIG += embed_translations

win32:QMAKE_CXXFLAGS += /bigobj


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target


DISTFILES += \
    GUI/style.qss

RESOURCES += \
    GUI/resource.qrc

RC_ICONS = app.ico
