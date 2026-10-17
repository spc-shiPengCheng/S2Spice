#include "./GUI/st_show_syz_parameter.h"
#include "ui_stshowsyzparameter.h"

#include <QLineSeries>
#include <QLogValueAxis>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>

#include "../Tool/st_struct_file.h"

StShowSYZParameter::StShowSYZParameter(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::StShowSYZParameter)
{
    ui->setupUi(this);

    // m_chart = new QChart();
    // ui->graphicsView->setChart(m_chart);
}

StShowSYZParameter::StShowSYZParameter(int row,
                                       int col,
                                       DisplayParameterType type,
                                       bool isLinearXAxis,
                                       bool isLinearYAxis,
                                       struct snp& sparameterData,
                                       DisplayDataType displayDataType,
                                       QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::StShowSYZParameter)
{

    ui->setupUi(this);

    initDialog();

    paintCurve(row, col, type, isLinearXAxis, isLinearYAxis, sparameterData, displayDataType);

}

StShowSYZParameter::~StShowSYZParameter()
{
    delete m_chartTop;
    if(m_chartTop != nullptr)
    {
        m_chartTop = nullptr;
    }

    delete m_chartBottom;
    if(m_chartBottom != nullptr)
    {
        m_chartBottom = nullptr;
    }

    delete ui;
}

void StShowSYZParameter::initDialog()
{
    // resize(578, 522);
}

void StShowSYZParameter::paintCurve(int row,
                                    int col,
                                    DisplayParameterType type,
                                    bool isLinearXAxis,
                                    bool isLinearYAxis,
                                    struct snp& sparameterData,
                                    DisplayDataType displayDataType)
{

    m_chartTop = new QChart();
    m_chartBottom = new QChart();

    std::vector<std::vector<std::vector<std::complex<double>>>> parameterDataRI;
    std::vector<std::vector<std::vector<std::complex<double>>>> parameterDataRIVF;
    std::vector<std::vector<std::vector<std::complex<double>>>> parameterDataRINetlist;
    std::vector<double> currentDataAbs;
    std::vector<std::complex<double>> currentData;
    QString title = "";
    if(type == DisplayParameterType::S)
    {
        parameterDataRI = sparameterData.sparametrRI;
        parameterDataRIVF = sparameterData.sparametrRIVF;
        parameterDataRINetlist = sparameterData.sparametrRINetlist;

        for(int i=0;i<(int)sparameterData.sparametrRI.size();i++)
        {
            currentDataAbs.push_back(std::abs(sparameterData.sparametrRI[i][row][col]));
        }

        title = "S";
    }
    else if(type == DisplayParameterType::Y)
    {
        parameterDataRI = sparameterData.yparametrRI;
        parameterDataRIVF = sparameterData.yparametrRIVF;
        parameterDataRINetlist = sparameterData.yparametrRINetlist;

        for(int i=0;i<(int)sparameterData.yparametrRI.size();i++)
        {
            currentDataAbs.push_back(std::abs(sparameterData.yparametrRI[i][row][col]));
        }

        title = "Y";
    }
    else if(type == DisplayParameterType::Z)
    {
        parameterDataRI = sparameterData.zparametrRI;
        parameterDataRIVF = sparameterData.zparametrRIVF;
        parameterDataRINetlist = sparameterData.sparametrRINetlist;

        for(int i=0;i<(int)sparameterData.zparametrRI.size();i++)
        {
            currentDataAbs.push_back(std::abs(sparameterData.zparametrRI[i][row][col]));
        }

        title = "Z";
    }
    else
    {
        parameterDataRI = sparameterData.sparametrRI;
        parameterDataRIVF = sparameterData.sparametrRIVF;
        parameterDataRINetlist = sparameterData.sparametrRINetlist;

        for(int i=0;i<(int)sparameterData.sparametrRI.size();i++)
        {
            currentDataAbs.push_back(std::abs(sparameterData.sparametrRI[i][row][col]));
        }

        title = "S";
    }
    setWindowTitle(title + QString("(%1,%2)").arg(row+1).arg(col+1));



    for(int i=0;i<(int)parameterDataRI.size();i++)
    {
        currentData.push_back(parameterDataRI[i][row][col]);
    }

    QLineSeries* seriesTop = new QLineSeries();
    QLineSeries* seriesBottom = new QLineSeries();

    std::vector<double> frequency = sparameterData.frequency;
    for(int i=0;i<frequency.size();i++)
    {
        frequency[i] = frequency[i] / sparameterData.frequencyUnit;
    }
    if(displayDataType == DisplayDataType::REAL_IMAG)
    {
        for(int i=0;i<(int)currentData.size();i++)
        {
            *seriesTop << QPointF(frequency[i],currentData[i].real());
            *seriesBottom << QPointF(frequency[i],currentData[i].imag());
        }

        seriesTop->setName(title + "data" + QString("(%1,%2)").arg(row+1).arg(col+1) + "-Real");
        seriesBottom->setName(title + "data" +QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Imag");
    }
    else if(displayDataType == DisplayDataType::MAG_ANGLE)
    {
        for(int i=0;i<(int)currentData.size();i++)
        {
            *seriesTop << QPointF(frequency[i],std::abs(currentData[i]));
            *seriesBottom << QPointF(frequency[i],std::atan2(currentData[i].imag(), currentData[i].real()) * 180/M_PI);
        }

        seriesTop->setName(title + "data" + QString("(%1,%2)").arg(row+1).arg(col+1) + "-Mag");
        seriesBottom->setName(title + "data" +QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Angle");
    }
    else if(displayDataType == DisplayDataType::DB_ANGLE)
    {
        for(int i=0;i<(int)currentData.size();i++)
        {
            *seriesTop << QPointF(frequency[i],20*log10(std::abs(currentData[i])));
            *seriesBottom << QPointF(frequency[i],std::atan2(currentData[i].imag(), currentData[i].real()) * 180/M_PI);
        }

        seriesTop->setName(title + "data" + QString("(%1,%2)").arg(row+1).arg(col+1) + "-dB");
        seriesBottom->setName(title + "data" +QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Angle");
    }
    else if(displayDataType == DisplayDataType::SMITH)
    {
        return;
    }
    m_chartTop->addSeries(seriesTop);
    m_chartBottom->addSeries(seriesBottom);


    QLineSeries* seriesTopVf = nullptr;
    QLineSeries* seriesBottomVf = nullptr;
    if(!parameterDataRIVF.empty())
    {
        seriesTopVf = new QLineSeries();
         seriesBottomVf = new QLineSeries();

        std::vector<std::complex<double>> currentData;
        for(int i=0;i<(int)parameterDataRIVF.size();i++)
        {
            currentData.push_back(parameterDataRIVF[i][row][col]);
        }


        if(displayDataType == DisplayDataType::REAL_IMAG)
        {
            for(int i=0;i<(int)currentData.size();i++)
            {
                *seriesTopVf << QPointF(frequency[i],currentData[i].real());
                *seriesBottomVf<< QPointF(frequency[i],currentData[i].imag());

                seriesTopVf->setName(title + "model" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Real");
                seriesBottomVf->setName(title + "model" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Imag");
            }
        }
        else if(displayDataType == DisplayDataType::MAG_ANGLE)
        {
            for(int i=0;i<(int)currentData.size();i++)
            {
                *seriesTopVf<< QPointF(frequency[i],std::abs(currentData[i]));
                *seriesBottomVf<< QPointF(frequency[i],std::atan2(currentData[i].imag(), currentData[i].real()) * 180/M_PI);

                seriesTopVf->setName(title + "model" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Mag");
                seriesBottomVf->setName(title + "model" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Angle");
            }
        }
        else if(displayDataType == DisplayDataType::DB_ANGLE)
        {
            for(int i=0;i<(int)currentData.size();i++)
            {
                *seriesTopVf << QPointF(frequency[i],20*log10(std::abs(currentData[i])));
                *seriesBottomVf<< QPointF(frequency[i],std::atan2(currentData[i].imag(), currentData[i].real()) * 180/M_PI);

                seriesTopVf->setName(title + "model" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-dB");
                seriesBottomVf->setName(title + "model" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Angle");
            }
        }
        else if(displayDataType == DisplayDataType::SMITH)
        {

        }


        m_chartTop->addSeries(seriesTopVf);
        m_chartBottom->addSeries(seriesBottomVf);
    }




    QLineSeries* seriesTopNetlist = nullptr;
    QLineSeries* seriesBottomNetlist = nullptr;
    if(!parameterDataRINetlist.empty())
    {
        seriesTopNetlist = new QLineSeries();
        seriesBottomNetlist = new QLineSeries();

        std::vector<std::complex<double>> currentData;
        for(int i=0;i<(int)parameterDataRINetlist.size();i++)
        {
            currentData.push_back(parameterDataRINetlist[i][row][col]);
        }


        if(displayDataType == DisplayDataType::REAL_IMAG)
        {
            for(int i=0;i<(int)currentData.size();i++)
            {
                *seriesTopNetlist << QPointF(frequency[i],currentData[i].real());
                *seriesBottomNetlist<< QPointF(frequency[i],currentData[i].imag());

                seriesTopNetlist->setName(title + "netlist" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Real");
                seriesBottomNetlist->setName(title + "netlist" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Imag");
            }
        }
        else if(displayDataType == DisplayDataType::MAG_ANGLE)
        {
            for(int i=0;i<(int)currentData.size();i++)
            {
                *seriesTopNetlist<< QPointF(frequency[i],std::abs(currentData[i]));
                *seriesBottomNetlist<< QPointF(frequency[i],std::atan2(currentData[i].imag(), currentData[i].real()) * 180/M_PI);

                seriesTopNetlist->setName(title + "netlist" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Mag");
                seriesBottomNetlist->setName(title + "netlist" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Angle");
            }
        }
        else if(displayDataType == DisplayDataType::DB_ANGLE)
        {
            for(int i=0;i<(int)currentData.size();i++)
            {
                *seriesTopNetlist << QPointF(frequency[i],20*log10(std::abs(currentData[i])));
                *seriesBottomNetlist<< QPointF(frequency[i],std::atan2(currentData[i].imag(), currentData[i].real()) * 180/M_PI);

                seriesTopNetlist->setName(title + "netlist" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-dB");
                seriesBottomNetlist->setName(title + "netlist" + QString("(%1,%2)").arg(row+1).arg(col+1)+ "-Angle");
            }
        }
        else if(displayDataType == DisplayDataType::SMITH)
        {

        }


        m_chartTop->addSeries(seriesTopNetlist);
        m_chartBottom->addSeries(seriesBottomNetlist);
    }




    // 假设我们有一个 QChartView 和一个 QChart
    QLogValueAxis *axisXTop = new QLogValueAxis();

    // 假设我们有一个 QChartView 和一个 QChart
    QLogValueAxis *axisXBottom = new QLogValueAxis();

    // 假设我们有一个 QChartView 和一个 QChart
    QValueAxis *axisXTopLinear = new QValueAxis();

    // 假设我们有一个 QChartView 和一个 QChart
    QValueAxis *axisXBottomLinear = new QValueAxis();

    if(!isLinearXAxis)
    {

        // 使用 'e' 或 'E' 格式，指定小数点后的位数
        // 例如: "%.1e" 会显示 1.0e+03, 1.0e+04 等
        axisXTop->setLabelFormat("%.1e");

        // 设置对数轴的范围和基数
        axisXTop->setBase(10.0);
        axisXTop->setMin((frequency[0]));
        axisXTop->setMax((frequency[frequency.size()-1]));
        axisXTop->setTitleText(QString("Frequency[%1]").arg(sparameterData.frequencyUnitLabel));


        m_chartTop->addAxis(axisXTop, Qt::AlignBottom);
        seriesTop->attachAxis(axisXTop);
        if(seriesTopVf != nullptr)
        {
             seriesTopVf->attachAxis(axisXTop);
        }




        // 使用 'e' 或 'E' 格式，指定小数点后的位数
        // 例如: "%.1e" 会显示 1.0e+03, 1.0e+04 等
        axisXBottom->setLabelFormat("%.1e");

        // 设置对数轴的范围和基数
        axisXBottom->setBase(10.0);
        axisXBottom->setMin((frequency[0]));
        axisXBottom->setMax((frequency[frequency.size()-1]));
        axisXBottom->setTitleText(QString("Frequency[%1]").arg(sparameterData.frequencyUnitLabel));


        m_chartBottom->addAxis(axisXBottom, Qt::AlignBottom);
        seriesBottom->attachAxis(axisXBottom);
        if(seriesBottomVf != nullptr)
        {
            seriesBottomVf->attachAxis(axisXTop);
        }
    }
    else
    {
        axisXTopLinear->setMin((frequency[0]));
        axisXTopLinear->setMax((frequency[frequency.size()-1]));
        axisXTopLinear->setTitleText(QString("Frequency[%1]").arg(sparameterData.frequencyUnitLabel));


        m_chartTop->addAxis(axisXTopLinear, Qt::AlignBottom);
        seriesTop->attachAxis(axisXTopLinear);
        if(seriesTopVf != nullptr)
        {
            seriesTopVf->attachAxis(axisXTop);
        }




        axisXBottomLinear->setMin((frequency[0]));
        axisXBottomLinear->setMax((frequency[frequency.size()-1]));
        axisXBottomLinear->setTitleText(QString("Frequency[%1]").arg(sparameterData.frequencyUnitLabel));


        m_chartBottom->addAxis(axisXBottomLinear, Qt::AlignBottom);
        seriesBottom->attachAxis(axisXBottomLinear);
        if(seriesBottomVf != nullptr)
        {
            seriesBottomVf->attachAxis(axisXTop);
        }
    }


    // 假设我们有一个 QChartView 和一个 QChart
    QLogValueAxis *axisYTop = new QLogValueAxis();

    // 假设我们有一个 QChartView 和一个 QChart
    QLogValueAxis *axisYBottom = new QLogValueAxis();


    // 假设我们有一个 QChartView 和一个 QChart
    QValueAxis *axisYTopLinear = new QValueAxis();

    // 假设我们有一个 QChartView 和一个 QChart
    QValueAxis *axisYBottomLinear = new QValueAxis();
    if(!isLinearYAxis)
    {

        // 使用 'e' 或 'E' 格式，指定小数点后的位数
        // 例如: "%.1e" 会显示 1.0e+03, 1.0e+04 等
        axisYTop->setLabelFormat("%e");


        // 设置对数轴的范围和基数
        axisYTop->setBase(10.0);
        // auto maxIt = std::max_element(currentDataAbs.begin(),currentDataAbs.end());
        // auto mminIt = std::min_element(currentDataAbs.begin(),currentDataAbs.end());
        // axisYTop->setMin((*mminIt));
        // axisYTop->setMax((*maxIt));

        // 将轴添加到图表...

        // m_chart->addAxis(axisY, Qt::AlignLeft);
        m_chartTop->addAxis(axisYTop, Qt::AlignLeft);
        seriesTop->attachAxis(axisYTop);
        if(seriesTopVf != nullptr)
        {
            seriesTopVf->attachAxis(axisXTop);
        }



        // 使用 'e' 或 'E' 格式，指定小数点后的位数
        // 例如: "%.1e" 会显示 1.0e+03, 1.0e+04 等
        axisYBottom->setLabelFormat("%e");


        // 设置对数轴的范围和基数
        axisYBottom->setBase(10.0);
        // auto maxIt1 = std::max_element(currentDataAbs.begin(),currentDataAbs.end());
        // auto mminIt1 = std::min_element(currentDataAbs.begin(),currentDataAbs.end());
        // axisYBottom->setMin((*mminIt1));
        // axisYBottom->setMax((*maxIt1));

        // 将轴添加到图表...

        m_chartBottom->addAxis(axisYBottom, Qt::AlignLeft);
        seriesBottom->attachAxis(axisYBottom);
        if(seriesBottomVf!= nullptr)
        {
            seriesBottomVf->attachAxis(axisXTop);
        }
    }
    else
    {
        // auto maxIt = std::max_element(currentDataAbs.begin(),currentDataAbs.end());
        // auto mminIt = std::min_element(currentDataAbs.begin(),currentDataAbs.end());
        // axisYTopLinear->setMin((*mminIt));
        // axisYTopLinear->setMax((*maxIt));

        // 将轴添加到图表...

        // m_chart->addAxis(axisY, Qt::AlignLeft);
        m_chartTop->addAxis(axisYTopLinear, Qt::AlignLeft);
        seriesTop->attachAxis(axisYTopLinear);
        if(seriesTopVf!= nullptr)
        {
            seriesTopVf->attachAxis(axisXTop);
        }



        // auto maxIt1 = std::max_element(currentDataAbs.begin(),currentDataAbs.end());
        // auto mminIt1 = std::min_element(currentDataAbs.begin(),currentDataAbs.end());
        // axisYBottomLinear->setMin((*mminIt1));
        // axisYBottomLinear->setMax((*maxIt1));

        // 将轴添加到图表...

        m_chartBottom->addAxis(axisYBottomLinear, Qt::AlignLeft);

        seriesBottom->attachAxis(axisYBottomLinear);
        if(seriesBottomVf!= nullptr)
        {
            seriesBottomVf->attachAxis(axisXTop);
        }
    }



    ui->graphicsView_top->setChart(m_chartTop);
    ui->graphicsView_bottom->setChart(m_chartBottom);

    ui->graphicsView_top->setRubberBand(QChartView::RectangleRubberBand);
    ui->graphicsView_bottom->setRubberBand(QChartView::RectangleRubberBand);
}
