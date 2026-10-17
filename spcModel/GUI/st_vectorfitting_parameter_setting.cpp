#include "./GUI/st_vectorfitting_parameter_setting.h"
#include "ui_st_vectorfitting_parameter_setting.h"

StVectorFittingParameterSetting::StVectorFittingParameterSetting(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::StVectorFittingParameterSetting)
{
    ui->setupUi(this);

    initSignal();
}

StVectorFittingParameterSetting::~StVectorFittingParameterSetting()
{
    delete ui;
}

bool StVectorFittingParameterSetting::getVectorFittingParameter(VectorFittingParameter &VFparameter)
{
    bool isDouble = false;
    double rms = ui->lineEdit_rms->text().trimmed().toDouble(&isDouble);
    if(!isDouble)
    {
        qDebug() << QString("Vector Fitting rms is not double!");
        VFparameter.rms = 0.001;
        return false;
    }
    else if(rms <= 0)
    {
        qDebug() << QString("Vector Fitting rms must greater than 0!");
        VFparameter.rms = 0.001;
        return false;
    }
    else
    {
        VFparameter.rms = ui->lineEdit_rms->text().trimmed().toDouble();
    }

    bool isInt = false;
    int minOrder = ui->lineEditMinOrder->text().trimmed().toInt(&isInt);
    if(!isInt)
    {
        qDebug() << QString("Vector Fitting minimum order is not int!");
        VFparameter.minOrder = 2;
        return false;
    }
    else if(minOrder <= 0)
    {
        qDebug() << QString("Vector Fitting minimum order must greater than 0!");
        VFparameter.minOrder = 2;
        return false;
    }
    else
    {
        VFparameter.minOrder = ui->lineEditMinOrder->text().trimmed().toInt();
    }

    int maxOrder = ui->lineEditMaxOrder->text().trimmed().toInt(&isInt);
    if(!isInt)
    {
        qDebug() << QString("Vector Fitting maximum order is not int!");
        VFparameter.maxOrder = 50;
        return false;
    }
    else if(maxOrder <=0)
    {
        qDebug() << QString("Vector Fitting maximum order must greater than 0!");
        VFparameter.maxOrder = 50;
        return false;
    }
    else if(maxOrder <= minOrder)
    {
        qDebug() << QString("Vector Fitting maximum order must greater than minimum order!");
        VFparameter.maxOrder = 50;
        return false;
    }
    else
    {
        VFparameter.maxOrder = ui->lineEditMaxOrder->text().trimmed().toInt();
    }

    return true;
}

void StVectorFittingParameterSetting::setVectorFittingParameter(VectorFittingParameter &VFparameter)
{

}

void StVectorFittingParameterSetting::initSignal()
{
    connect(ui->pushButton_run,&QPushButton::clicked,this,&StVectorFittingParameterSetting::oK);
}

void StVectorFittingParameterSetting::oK()
{
    accept();
}
