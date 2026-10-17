#pragma once

#include <QDialog>

class QCustomPlot;

namespace Ui {
class StShowSParameterSingularDialog;
}

class StShowSParameterSingularDialog : public QDialog
{
    Q_OBJECT

public:
    explicit StShowSParameterSingularDialog(QWidget *parent = nullptr);
    ~StShowSParameterSingularDialog();

    void setGraphData(const std::vector<double>& freq, const std::vector<double>& singular);

private:
    void initDialog();

    void initStyleSheet();

    void initSignals();

private:
    Ui::StShowSParameterSingularDialog *ui;
    QCustomPlot* m_customPlot = nullptr;
};

