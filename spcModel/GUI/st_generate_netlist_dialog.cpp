#include "st_generate_netlist_dialog.h"
#include "ui_st_generate_netlist_dialog.h"
#include <QButtonGroup>

StGenerateNetlistDialog::StGenerateNetlistDialog(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::StGenerateNetlistDialog)
{
    ui->setupUi(this);

    setWindowTitle("Generate Netlist");

    QButtonGroup *group = new QButtonGroup(this);
    group->addButton(ui->radioButton_controlledSource);
    group->addButton(ui->radioButton_RLC);
    group->setExclusive(true); // 默认就是 true

    ui->radioButton_controlledSource->setChecked(true);

    connect(ui->pushButton_ok, &QPushButton::clicked,this,&StGenerateNetlistDialog::on_ok);
}

StGenerateNetlistDialog::~StGenerateNetlistDialog()
{
    delete ui;
}

void StGenerateNetlistDialog::getNetlistType(NetlistType netlistType)
{
    if(ui->radioButton_controlledSource->isChecked())
    {
        netlistType = NetlistType::ControlledSource;
    }
    else if(ui->radioButton_RLC->isChecked())
    {
        netlistType = NetlistType::RLC;
    }
    else
    {
        netlistType = NetlistType::ControlledSource;
    }
}

void StGenerateNetlistDialog::on_ok()
{
    accept();
}
