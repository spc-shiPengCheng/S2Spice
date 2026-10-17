#include "st_show_sparameter_singular_dialog.h"
#include "ui_st_show_sparameter_singular_dialog.h"
#include "qcustomplot.h"


StShowSParameterSingularDialog::StShowSParameterSingularDialog(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::StShowSParameterSingularDialog)
{
    ui->setupUi(this);

    initDialog();

    initSignals();

    initStyleSheet();
}

StShowSParameterSingularDialog::~StShowSParameterSingularDialog()
{
    delete ui;
}

void StShowSParameterSingularDialog::setGraphData(const std::vector<double>& freq, const std::vector<double>& singular)
{
    QCPTextElement *title = new QCPTextElement(m_customPlot, "Singular values", QFont("Arial", 16, QFont::Bold));
    m_customPlot->plotLayout()->insertRow(0);
    m_customPlot->plotLayout()->addElement(0, 0, title);

    // -------- 1. 添加第一条曲线 ----------
    m_customPlot->addGraph();  // 创建一个新的曲线对象，由 QCustomPlot 管理

    // -------- 2. 设置曲线的样式 ----------
    // 设置曲线颜色和粗细，这里使用蓝色，线宽为 2
    m_customPlot->graph(0)->setPen(QPen(Qt::blue, 2));
    // （可选）设置曲线下方填充颜色，这里用半透明的蓝色
    // m_customPlot->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20)));

    // -------- 3. 将数据设置到曲线上 ----------
    m_customPlot->graph(0)->setData(QVector<double>(freq.begin(), freq.end()), QVector<double>(singular.begin(), singular.end()));

    m_customPlot->addGraph();
    QVector<double>refLineFreq;
    refLineFreq.push_back(freq[0]);
    refLineFreq.push_back(freq[freq.size()-1]);
    QVector<double>refLineSingular;
    refLineSingular.push_back(1);
    refLineSingular.push_back(1);
    m_customPlot->graph(1)->setPen(QPen(Qt::green, 2));
    m_customPlot->graph(1)->setData(refLineFreq, refLineSingular);

    // -------- 4. 设置坐标轴标签 ----------
    m_customPlot->xAxis->setLabel("Frequency");

    // -------- 5. 自动调整坐标轴范围以显示所有数据 ----------
    // m_customPlot->rescaleAxes();
    m_customPlot->xAxis->rescale();

    // -------- 6. 重绘刷新界面 ----------
    m_customPlot->replot();

}

void StShowSParameterSingularDialog::initDialog()
{
    setWindowTitle("Singular");
    resize(578, 522);
    m_customPlot = ui->customPlot;
}

void StShowSParameterSingularDialog::initStyleSheet()
{

}

void StShowSParameterSingularDialog::initSignals()
{

}
