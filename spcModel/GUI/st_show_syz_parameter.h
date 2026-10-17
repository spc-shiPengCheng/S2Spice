#ifndef ST_SHOW_SYZ_PARAMETER_H
#define ST_SHOW_SYZ_PARAMETER_H

#include <QChart>
#include <QDialog>
#include "../Tool/st_struct_file.h"

namespace Ui {
class StShowSYZParameter;
}

class StShowSYZParameter : public QDialog
{
    Q_OBJECT

public:
    explicit StShowSYZParameter(QWidget *parent = nullptr);
    explicit StShowSYZParameter(int row,
                                int col,
                                DisplayParameterType type,
                                bool isLinearXAxis,
                                bool isLinearYAxis,
                                struct snp& sparameterData,
                                DisplayDataType displayDataType,
                                QWidget *parent = nullptr);
    ~StShowSYZParameter();

    void initDialog();

private:
    void paintCurve(int row,
                    int col,
                    DisplayParameterType type,
                    bool isLinearXAxis,
                    bool isLinearYAxis,
                    struct snp& sparameterData,
                    DisplayDataType displayDataType);

private:
    Ui::StShowSYZParameter *ui;

    QChart* m_chartTop;

    QChart* m_chartBottom;
};

#endif // ST_SHOW_SYZ_PARAMETER_H
