#include "st_setting_spice_solver_dialog.h"
#include "ui_st_setting_spice_solver_dialog.h"

#include <QFileDialog>

StSettingSpiceSolverDialog::StSettingSpiceSolverDialog(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::StSettingSpiceSolverDialog)
{
    ui->setupUi(this);

    ui->lineEdit_netlist->setReadOnly(true);
    ui->lineEdit_spiceSolver->setReadOnly(true);

    connect(ui->pushButton_spiceSolver, &QPushButton::clicked,this,&StSettingSpiceSolverDialog::on_loadSpiceSolver);

    connect(ui->pushButton_netlist, &QPushButton::clicked,this,&StSettingSpiceSolverDialog::on_loadNetlist);

    connect(ui->pushButton_ok, &QPushButton::clicked,this,&StSettingSpiceSolverDialog::on_ok);
}

StSettingSpiceSolverDialog::~StSettingSpiceSolverDialog()
{
    delete ui;
}

void StSettingSpiceSolverDialog::setNgspiceAndNetlistPath(QString ngspicePath,QString netlistPath)
{
    m_spicePath = ngspicePath;
    m_netlistPath = netlistPath;
    ui->lineEdit_netlist->setText(m_netlistPath);
    ui->lineEdit_spiceSolver->setText(m_spicePath);
}

QString StSettingSpiceSolverDialog::getNgspicePath()
{
    return m_spicePath;
}

QString StSettingSpiceSolverDialog::getNetlistPath()
{
    return m_netlistPath;
}

void StSettingSpiceSolverDialog::on_loadSpiceSolver()
{
    // 1. 打开文件选择对话框
    QString filePath = QFileDialog::getOpenFileName(
        this,                   // 父窗口
        "导入文件",              // 对话框标题
        "",                     // 默认打开目录（空字符串表示系统默认）
        "执行文件 (*.exe);"  // 文件过滤器
        );

    // 2. 用户点击取消时返回空字符串
    if (filePath.isEmpty()) {
        return;
    }
    ui->lineEdit_spiceSolver->setText(filePath);

    m_spicePath = filePath;
}

void StSettingSpiceSolverDialog::on_loadNetlist()
{
    // 1. 打开文件选择对话框
    QString filePath = QFileDialog::getOpenFileName(
        this,                   // 父窗口
        "导入文件",              // 对话框标题
        "",                     // 默认打开目录（空字符串表示系统默认）
        "所有文件 (*);;网表文件 (*.cir);"  // 文件过滤器
        );

    // 2. 用户点击取消时返回空字符串
    if (filePath.isEmpty()) {
        return;
    }

    ui->lineEdit_netlist->setText(filePath);
    m_netlistPath = filePath;
}

void StSettingSpiceSolverDialog::on_ok()
{
    accept();
}
