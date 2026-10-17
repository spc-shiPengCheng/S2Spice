#ifndef ST_VECTORFITTING_PARAMETER_SETTING_H
#define ST_VECTORFITTING_PARAMETER_SETTING_H

#include <QDialog>

#include "../Tool/st_struct_file.h"

namespace Ui {
class StVectorFittingParameterSetting;
}

class StVectorFittingParameterSetting : public QDialog
{
    Q_OBJECT

public:
    explicit StVectorFittingParameterSetting(QWidget *parent = nullptr);
    ~StVectorFittingParameterSetting();

    bool getVectorFittingParameter(struct VectorFittingParameter& VFparameter);

    void setVectorFittingParameter(struct VectorFittingParameter& VFparameter);

private:
    void initSignal();

private slots:
    void oK();

private:
    Ui::StVectorFittingParameterSetting *ui;
};

#endif // ST_VECTORFITTING_PARAMETER_SETTING_H
