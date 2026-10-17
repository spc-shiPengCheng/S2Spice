#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "../Tool/st_struct_file.h"
#include <QFutureWatcher>
#include <atomic>

class QProcess;
class QFile;

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    //Init dialog.
    void initDialog();

    //Init signal.
    void initSignal();

    //Open VectorFittingUI dialog.
    void openVectorFittingDialog();

    void initPortTable(int portNum);

    void initSparameterInformation();

    void runSpiceFileSimulation(QString spiceFilePath, std::vector<std::vector<std::vector<std::complex<double>>>>& sparametrRISpice,int setPortIndex);

private slots:
    //Import S parameter.
    void open();

    //Export S parameter
    void onExportOriginSnpFile();
    void onExportVFSnpFile();

    void showSparameter(int row,int col);

    void vectorfitting();

    void passivityTestData();
    void passivityTestModel();

    void causalityTestData();

    void causality();

    void passivityEnforceAction();

    void logInformation(const QString& message);

    void switchCurrentSnp();

    void onDeleteData();

    void onDeleteModel();

    void onSetSpiceNetlist();

    void onRunSpice(int portIndex);

    void onGenerateSpiceNetlist();

    void onDeleteTempSpiceFile();

    void showModelInformation();

    void on_stop();

    void onTaskFinished();


private:
    Ui::MainWindow *ui;

    QString sparameterFilePath;

    struct snp sparameterData;

    std::vector<snp> sparameterDataSet;
    QProcess* process;
    QString Path = "";

    QFutureWatcher<void> m_watcher;
    std::atomic<bool>    m_cancelled{false};
    QString m_netlistPath;
    QString m_spicePath;
};
#endif // MAINWINDOW_H
