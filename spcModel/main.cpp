#include "./GUI/mainwindow.h"

#include <QApplication>
#include <QLocale>
#include <QTranslator>
#include <QFile>
#include <QStyleFactory>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    // 设置全局窗口图标
    a.setWindowIcon(QIcon(":/icon/app.png"));

    // // 在 main 函数中设置应用的基础风格
    // QApplication::setStyle(QStyleFactory::create("Fusion"));

    // 加载样式表文件
    QFile styleFile(":/style.qss"); // 资源路径
    if (styleFile.open(QFile::ReadOnly)) {
        QString styleSheet = QLatin1String(styleFile.readAll());
        a.setStyleSheet(styleSheet); // 全局样式表生效！
        styleFile.close();
    } else {
        // 文件打开失败处理
        qDebug() << "Failed to load style sheet.";
    }

    QTranslator translator;
    const QStringList uiLanguages = QLocale::system().uiLanguages();
    for (const QString &locale : uiLanguages) {
        const QString baseName = "S2Spice_" + QLocale(locale).name();
        if (translator.load(":/i18n/" + baseName)) {
            a.installTranslator(&translator);
            break;
        }
    }
    MainWindow w;
    w.show();
    return a.exec();
}
