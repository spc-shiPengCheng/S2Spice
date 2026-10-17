#ifndef ST_GENERATE_NETLIST_DIALOG_H
#define ST_GENERATE_NETLIST_DIALOG_H

#include <QDialog>
#include "../Tool/st_struct_file.h"

namespace Ui {
class StGenerateNetlistDialog;
}

class StGenerateNetlistDialog : public QDialog
{
    Q_OBJECT

public:
    explicit StGenerateNetlistDialog(QWidget *parent = nullptr);
    ~StGenerateNetlistDialog();

    void getNetlistType(NetlistType netlistType);

private slots:
        void on_ok();

private:
    Ui::StGenerateNetlistDialog *ui;
};

#endif // ST_GENERATE_NETLIST_DIALOG_H
