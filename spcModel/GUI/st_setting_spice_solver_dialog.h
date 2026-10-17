#ifndef ST_SETTING_SPICE_SOLVER_DIALOG_H
#define ST_SETTING_SPICE_SOLVER_DIALOG_H

#include <QDialog>

namespace Ui {
class StSettingSpiceSolverDialog;
}

class StSettingSpiceSolverDialog : public QDialog
{
    Q_OBJECT

public:
    explicit StSettingSpiceSolverDialog(QWidget *parent = nullptr);
    ~StSettingSpiceSolverDialog();

    void setNgspiceAndNetlistPath(QString ngspicePath,QString netlistPath);

    QString getNgspicePath();
    QString getNetlistPath();

private slots:
    void on_loadSpiceSolver();
    void on_loadNetlist();
    void on_ok();

private:
    Ui::StSettingSpiceSolverDialog *ui;

    QString m_netlistPath;
    QString m_spicePath;
};

#endif // ST_SETTING_SPICE_SOLVER_DIALOG_H
