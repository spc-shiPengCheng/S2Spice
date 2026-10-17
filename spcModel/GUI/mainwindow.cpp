#include "./GUI/mainwindow.h"
#include "ui_mainwindow.h"

#include "../Library/VectorFitting/vectorfitting.h"
#include "../Library/VectorFitting/check_sparameter.h"

#include "../Tool/st_file_io.h"
#include "../Tool/st_convert_syz.h"
#include "./st_show_syz_parameter.h"
#include "./st_show_sparameter_singular_dialog.h"
#include "./st_vectorfitting_parameter_setting.h"
#include "./st_setting_spice_solver_dialog.h"
#include "./st_generate_netlist_dialog.h"
#include "../../Tool/st_logger.h"
#include "../Library/VectorFitting/spice_netlist_gengerate.h"

#include <QChartView>
#include <QFileDialog>
#include <QLineSeries>
#include <QLogValueAxis>
#include <QValueAxis>
#include <QMessageBox>
#include <QCategoryAxis>
#include <QProcess>
#include <QTemporaryFile>
#include <QRegularExpression>
#include <Windows.h>
#include <QThread>
#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>
using namespace std;
// #include <process.h>

static QStringList readColumn(const QString &filePath, int columnIndex,
                              const QString &separator = " ")
{
    QStringList columnData;
    QFile file(filePath);

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qWarning() << "无法打开文件:" << filePath;
        qDebug() << "Failed to open file:" << file.errorString();
        return columnData;
    }

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty()) continue;  // 跳过空行

        QStringList parts = line.split(separator, Qt::SkipEmptyParts);
        if (columnIndex < parts.size()) {
            columnData.append(parts.at(columnIndex));
        } else {
            // 该行列数不足，可以选择跳过或添加空字符串
            columnData.append(QString());
        }
    }
    file.close();
    return columnData;
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    initDialog();

    initSignal();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::initDialog()
{
    this->setWindowTitle(tr("spcModel"));

    // resize(1000, 820);
    setMinimumSize(1098,815);


    ui->comboBox_DataType->addItem("Real/Imag", DisplayDataType::REAL_IMAG);
    ui->comboBox_DataType->addItem("Mag/Angle", DisplayDataType::MAG_ANGLE);
    ui->comboBox_DataType->addItem("dB/Angle", DisplayDataType::DB_ANGLE);
    ui->comboBox_DataType->addItem("Smith", DisplayDataType::SMITH);


    ui->label_fittingError->setText("");

    process = new QProcess(this);
}

void MainWindow::initSignal()
{
    connect(ui->actionTouchstoneImport, SIGNAL(triggered()), this, SLOT (open()));

    connect(ui->actionTouchstoneExportData, SIGNAL(triggered()), this, SLOT(onExportOriginSnpFile()));

    connect(ui->actionTouchstoneExportModel, SIGNAL(triggered()), this, SLOT(onExportVFSnpFile()));

    connect(ui->actionCheckPassivityData, SIGNAL(triggered()), this, SLOT(passivityTestData()));

    connect(ui->actionCheckCausalityData, SIGNAL(triggered()), this, SLOT(causalityTestData()));

    connect(ui->actionVF, SIGNAL(triggered()), this, SLOT (vectorfitting()));

    connect(ui->actionCheckPassivityModel, SIGNAL(triggered()), this, SLOT (passivityTestModel()));

    connect(ui->actionPassivity_Enforce, SIGNAL(triggered()), this, SLOT (passivityEnforceAction()));

    connect(ui->tableWidget_port, &QTableWidget::cellClicked, this, &MainWindow::showSparameter);

    connect(StLogger::instance(), &StLogger::logMessage, this, &MainWindow::logInformation);

    connect(ui->comboBox_data, &QComboBox::currentIndexChanged, this, &MainWindow::switchCurrentSnp);

    connect(ui->comboBox_model, &QComboBox::currentIndexChanged, this, &MainWindow::showModelInformation);

    connect(ui->actionDelete, SIGNAL(triggered()), this, SLOT (onDeleteData()));

    connect(ui->actionDeleteModel, SIGNAL(triggered()), this, SLOT (onDeleteModel()));

    connect(ui->actionRun_Spice, SIGNAL(triggered()), this, SLOT (onSetSpiceNetlist()));

    connect(ui->actionspice, SIGNAL(triggered()), this, SLOT(onGenerateSpiceNetlist()));

    connect(ui->actionStop, SIGNAL(triggered()), this, SLOT(on_stop()));
}

void MainWindow::initPortTable(int portNum)
{
    ui->tableWidget_port->setRowCount(portNum);
    ui->tableWidget_port->setColumnCount(portNum);

    for(int row=0;row<portNum;row++)
    {
        for(int col = 0;col<portNum;col++)
        {
            ui->tableWidget_port->setItem(row,col,new QTableWidgetItem(QString("S(%1,%2)").arg(row+1).arg(col+1)));
        }
    }
}

void MainWindow::initSparameterInformation()
{
    ui->plainTextEdit_dataInfo->clear();

    ui->plainTextEdit_dataInfo->appendPlainText(QString("Port Num: %1").arg(sparameterData.portNum));
    ui->plainTextEdit_dataInfo->appendPlainText(QString("Reference Impedance: %1").arg(sparameterData.referenceZ));
    ui->plainTextEdit_dataInfo->appendPlainText(QString("Max Frequency: %1 Hz").arg(*std::max_element(sparameterData.frequency.begin(), sparameterData.frequency.end())));
    ui->plainTextEdit_dataInfo->appendPlainText(QString("Min Frequency: %1 Hz").arg(*std::min_element(sparameterData.frequency.begin(), sparameterData.frequency.end())));
    ui->plainTextEdit_dataInfo->appendPlainText(QString("Frequency Points Num: %1").arg(sparameterData.frequency.size()));
    ui->plainTextEdit_dataInfo->appendPlainText(QString("Reciprocal: %1").arg(isReciprocal(sparameterData.sparametrRI) ? "True":"False"));
    ui->plainTextEdit_dataInfo->appendPlainText(QString("Passive: %1").arg(isPassive(sparameterData.sparametrRI) ? "True":"False"));
    ui->plainTextEdit_dataInfo->appendPlainText(QString("Causal: %1").arg(isCausal(sparameterData.sparametrRI) ? "True":"False"));
    ui->plainTextEdit_dataInfo->appendPlainText(QString("Size: %1 byte").arg(sparameterData.fileSize));
    ui->plainTextEdit_dataInfo->appendPlainText(QString("Snp: s%1p file").arg(sparameterData.portNum));
}

void MainWindow::runSpiceFileSimulation(QString spiceFilePath,
                                        std::vector<std::vector<std::vector<std::complex<double>>>>& sparametrRISpice,
                                        int setPortIndex)
{

    connect(process, &QProcess::finished, this, &MainWindow::onDeleteTempSpiceFile);

    sparametrRISpice.clear();


    QFile input(spiceFilePath);
    if (!input.open(QIODevice::ReadOnly | QIODevice::Text))
        return;

    QTextStream inGetParam(&input);
    QStringList portList  ={};
    QString netlistName = QFileInfo(spiceFilePath).completeBaseName();
    while (!inGetParam.atEnd()) {
        QString line = inGetParam.readLine().simplified();

        if (line.contains(".SUBCKT " + netlistName)) {
            QStringList mainSubckt= line.split(" ");
            portList = mainSubckt.mid(2);
        }
    }
    int portNum = portList.size();
    sparametrRISpice.resize(201,std::vector<std::vector<std::complex<double>>>(portNum, std::vector<std::complex<double>>(portNum,0)));

    int currentPortIndex = setPortIndex;
    QString portName = portList[setPortIndex];
    // for(QString portName : portList)
    {
        QString portIndex = portName.mid(1);
        // 创建临时文件
        Path = QFileInfo(spiceFilePath).absolutePath() + "/tempSpice.cir";
        QFile tempFile(QFileInfo(spiceFilePath).absolutePath() + "/tempSpice.cir");

        if (!tempFile.open(QIODevice::ReadWrite | QIODevice::Text))
            return;

        // 将文件指针移回文件开头
        input.seek(0);
        // 重要：QTextStream 内部缓冲区需要重置，重新关联文件
        // in.seek(0);  // 或者直接重新构造 QTextStream: QTextStream in2(&file);
        QTextStream in(&input);
        QTextStream out(&tempFile);

        while (!in.atEnd()) {
            QString line = in.readLine();
            if (line.contains(".SUBCKT " + QFileInfo(spiceFilePath).completeBaseName())) {
                out << "Vin" << portIndex <<" " << portName << " 0 DC 0 AC 1" << "\n";
                for(int portIndexCol = 0;portIndexCol< portList.size();portIndexCol++)
                {
                    if(portIndex != QString::number(portIndexCol+1))
                    {
                        out << "Vin" << portIndexCol+1 <<" " << "p" << portIndexCol+1 << " 0 DC 0 AC 0" << "\n";
                    }
                }
            }
            else if(line.contains(".ENDS " + QFileInfo(spiceFilePath).completeBaseName()))
            {
                out << ".control" << "\n";
                out << QString(tr("ac lin 201 %1 %2").arg(sparameterData.frequency[0]).arg(sparameterData.frequency[sparameterData.frequency.size()-1])) << "\n";
                // out << QString(tr("ac lin 201 1e7 8e9")) << "\n";
                // out << "ac lin 201 1e7 6e8" << "\n";
                out << "let Z0 = 50"  << "\n";
                // out << "let Vport = v(" << portName << ")" << "\n";
                // out << "let Iin = i(F"<< portIndex << << ")" << "\n";
                // out << "let a1 = (Vport + Z0*Iin) / (2*sqrt(Z0))" << "\n";
                // out << "let b1 = (Vport - Z0*Iin) / (2*sqrt(Z0))" << "\n";
                // out << "let S11 = b1 / a1   " << "\n";
                // out << "wrdata s11_result.txt real(S11) imag(S11)" << "\n";


                // out << "let Vport2 = v(" << portName << ")" << "\n";
                // out << "let Iin2 = i(V"<< 2 << ")" << "\n";
                // out << "let a2 = (Vport2 + Z0*Iin2) / (2*sqrt(Z0))" << "\n";
                // out << "let b2 = (Vport2 - Z0*Iin2) / (2*sqrt(Z0))" << "\n";
                // out << "let S21 = b2 / a2   " << "\n";
                // out << "wrdata s21_result.txt real(S21) imag(S21)" << "\n";

                for(int portIndexCol = 0;portIndexCol< portList.size();portIndexCol++)
                {
                    if(QString::number(portIndexCol+1) == portIndex)
                    {
                        out << QString(tr("let Vport%1 = v(%2)").arg(portIndexCol+1).arg(portName)) << "\n";
                        out << QString(tr("let Iin%1 = i(V%2)").arg(portIndexCol+1).arg(portIndexCol+1)) << "\n";
                        out << QString(tr("let a%1 = (Vport%1 + Z0*Iin%1) / (2*sqrt(Z0))").arg(portIndexCol+1)) << "\n";
                        out << QString(tr("let b%1 = (Vport%1 - Z0*Iin%1) / (2*sqrt(Z0))").arg(portIndexCol+1)) << "\n";
                        out << QString(tr("let S%1%2 = b%1 / a%1   ").arg(portIndexCol+1).arg(portIndex)) << "\n";
                        out << QString(tr("wrdata s%1%2_result.txt real(S%1%2) imag(S%1%2)").arg(portIndexCol+1).arg(portIndex)) << "\n";
                    }
                    else
                    {
                        out << QString(tr("let Vport%1 = v(%2)").arg(portIndexCol+1).arg(portName)) << "\n";
                        out << QString(tr("let Iin%1 = i(V%2)").arg(portIndexCol+1).arg(portIndexCol+1)) << "\n";
                        out << QString(tr("let a%1 = (Vport%1) / (2*sqrt(Z0))").arg(portIndexCol+1)) << "\n";
                        out << QString(tr("let b%1 = (-Z0*Iin%1) / (2*sqrt(Z0))").arg(portIndexCol+1)) << "\n";
                        out << QString(tr("let S%1%2 = b%1 / a%1   ").arg(portIndexCol+1).arg(portIndex)) << "\n";
                        out << QString(tr("wrdata s%1%2_result.txt real(S%1%2) imag(S%1%2)").arg(portIndexCol+1).arg(portIndex)) << "\n";
                    }

                }


                out << ".endc" << "\n";
                out << ".end" << "\n";
            }
            else
            {
                out << line << "\n";
            }
        }


        tempFile.close();



        QStringList commond;

        commond << "-b"  << tempFile.fileName();
        qDebug() << tempFile.fileName();

        // process->setCreateProcessArgumentsModifier(
        //     [](QProcess::CreateProcessArguments *args)
        //     {
        //         // 拿到 Windows 启动信息结构体
        //         LPSTARTUPINFO si = args->startupInfo;
        //         // 关键1：标记使用 wShowWindow 控制窗口显示
        //         si->dwFlags |= STARTF_USESHOWWINDOW;
        //         // 关键2：设置窗口为隐藏
        //         si->wShowWindow = SW_HIDE;
        //     }
        //     );

        // // ===================== 2. 进程通道配置（必须设置） =====================
        // // 独立管道，不共享控制台、不转发窗口
        // process->setProcessChannelMode(QProcess::SeparateChannels);





        QString program = m_spicePath;
        process->start(program, commond);


        if(!process->waitForStarted(3000)) // 3秒超时判断
        {
            qDebug() << "ngspice 启动失败：" << process->errorString();
            process->deleteLater();
        }



        process->waitForFinished(-1);


        QString appDir = QCoreApplication::applicationDirPath();


        QFileInfo fi(appDir);
        // ===== 关键步骤：刷新目录缓存 =====
        QDir dir = fi.absoluteDir();
        dir.refresh();   // 强制重读目录

        // QThread::sleep(5);

        for(int portIndexCol = 0;portIndexCol< portList.size();portIndexCol++)
        {
            QString appDirNew = appDir + QString(tr("/s%1%2_result.txt").arg(portIndexCol+1).arg(portIndex));

            QStringList res1 = readColumn(appDirNew, 0, " ");
            QStringList res2 = readColumn(appDirNew, 1, " ");
            QStringList res4 = readColumn(appDirNew, 3, " ");
            for(int i=0;i<res2.size();i++)
            {
                sparametrRISpice[i][portIndexCol][currentPortIndex] = std::complex<double>(res2[i].toDouble() , res4[i].toDouble());
            }
            QFile::remove(appDirNew);
        }
        currentPortIndex++;
    }

    input.close();
}

void MainWindow::vectorfitting()
{
    if(sparameterData.sparametrRI.empty())
    {
        return;
    }

    StVectorFittingParameterSetting* StVectorFittingParameterSettingPtr = new StVectorFittingParameterSetting(this);

    int result = StVectorFittingParameterSettingPtr->exec();

    if(result == QDialog::Accepted)
    {

        m_cancelled.store(false);   // 复位
        QFuture<void> f = QtConcurrent::run([=]() {

            // 耗时计算，不能操作 UI
            struct VectorFittingParameter vfparameter;
            bool flag = StVectorFittingParameterSettingPtr->getVectorFittingParameter(vfparameter);
            if(!flag)
            {
                QMessageBox msgBox;
                msgBox.setWindowTitle("Info");
                msgBox.setText(QString("Vector fitting parameter illegal."));

                msgBox.setIcon(QMessageBox::Information);

                // 手动添加自定义按钮，并分配角色 (YesRole 或 NoRole)
                QPushButton *customYesButton = msgBox.addButton("Ok", QMessageBox::YesRole);

                // 设置默认按钮
                msgBox.setDefaultButton(customYesButton);

                // 显示对话框并获取用户操作结果
                msgBox.exec();
                return;
            }


            struct VectorFittingInputData vfInputData;

            vfInputData.portNum = sparameterData.portNum;
            vfInputData.frequency = sparameterData.frequency;
            vfInputData.sparametrRI = sparameterData.sparametrRI;
            // vfInputData.realPolesNum = ceil(vfparameter.minOrder/2);
            // vfInputData.complexPolesNum = ceil(vfparameter.minOrder/2/2);
            vfInputData.realPolesNum = 2;
            vfInputData.complexPolesNum = 2;
            double realSum = 0;
            double imagSum = 0;
            for(int i=0;i<sparameterData.portNum;i++)
            {
                for(int j=0;j<sparameterData.portNum;j++)
                {
                    for(int k=0;k<sparameterData.frequency.size();k++)
                    {
                        realSum += std::abs(sparameterData.sparametrRI[k][i][j].real());
                        imagSum += std::abs(sparameterData.sparametrRI[k][i][j].imag());
                    }
                }
            }
            vfInputData.complexPolesRealCoefficient = imagSum/realSum;
            vfInputData.isComputeComplexPoles = false;
            vfInputData.dataPrecision = 1e-2;
            vfInputData.isProportional = vfparameter.isProportional;
            vfInputData.isConstant = vfparameter.isConstant;
            vfInputData.isNormalFrequency = vfparameter.isFrequencyNormalized;
            vfInputData.iterationNum = vfparameter.iterations;
            vfInputData.fittingError = vfparameter.rms;

            vfInputData.isWeight = vfparameter.isWeigth;
            vfInputData.isDirectComputePoles;
            std::vector<double> singular = getSingular(sparameterData.sparametrRI);
            double maxSingular = *max_element(singular.begin(),singular.end());
            vfInputData.isEnforcePolesStable = maxSingular > 1.1 ? false:true;
            vfInputData.isRelax = vfparameter.isRelaxed;
            std::vector<std::complex<double>> poles;
            std::vector<std::vector<std::complex<double>>> residuse;
            std::vector<double> constant;
            std::vector<double> proportional;
            vfInputData.poles1 = poles;
            vfInputData.residuse = residuse;
            vfInputData.proportional = proportional;
            vfInputData.constant = constant;



            int iterationNum = 1;
            std::vector<int> iteration;
            std::vector<int> order;
            std::vector<double> vfError;


            QElapsedTimer timer;
            timer.start();


            double error = VFImprove(vfInputData);

            iteration.push_back(iterationNum);
            order.push_back(vfInputData.realPolesNum + vfInputData.complexPolesNum*2);
            vfError.push_back(error);



            int maxComplexPolesNum = 0;
            std::vector<std::complex<double>> maxComplexPolesReference;
            computeComplexPolesNum(vfInputData.frequency,
                                   vfInputData.sparametrRI,
                                   1e-2,
                                   1e-2,
                                   maxComplexPolesNum,
                                   maxComplexPolesReference);


            while((error > vfInputData.fittingError
                    && (vfInputData.realPolesNum + vfInputData.complexPolesNum*2) < vfparameter.maxOrder)
                   || std::isnan(error))
            {
                if (m_cancelled.load()) {
                    qDebug() << "任务被取消.";
                    return;   // 从这里安全退出，析构会正常执行
                }

                if(maxComplexPolesNum > vfInputData.complexPolesNum)
                {
                    vfInputData.complexPolesNum = vfInputData.complexPolesNum + 1;
                }
                else
                {
                    vfInputData.realPolesNum++;
                }

                iterationNum++;

                error = VFImprove(vfInputData);

                StLogger::instance()->log(QString("The rms of iteration: %1").arg(error));

                iteration.push_back(iterationNum);
                order.push_back(vfInputData.realPolesNum + vfInputData.complexPolesNum*2);
                vfError.push_back(error);
            }


            qint64 elapsed = timer.elapsed();   // 返回毫秒，自 start 以来
            qDebug() << "拟合耗时:" << elapsed << "ms";


            QMetaObject::invokeMethod(this, [=]{
                // 1. 创建数据序列（一条曲线）
                QLineSeries *series = new QLineSeries();
                series->setName("Order&Error");
                for(int i=0;i<iteration.size();i++)
                {
                    *series << QPointF(iteration[i], vfError[i]);
                }
                series->setPointsVisible(true);   // 可选：显示数据点

                // 1. 创建数据序列（一条曲线）
                QLineSeries *seriesReference = new QLineSeries();
                seriesReference->setName("Reference error");

                *seriesReference << QPointF(iteration[0]-1, vfInputData.fittingError);
                *seriesReference << QPointF(iteration[iteration.size()-1] + 2, vfInputData.fittingError);


                // 2. 创建图表并添加序列
                QChart *chart = new QChart();
                chart->addSeries(series);
                chart->addSeries(seriesReference);
                chart->setTitle("Iteration process");

                // 3. 创建并配置底部 X 轴（主 X 轴）
                QCategoryAxis *axisX = new QCategoryAxis;
                axisX->setTitleText("Iterations");
                axisX->setRange(iteration[0]-1, iteration[iteration.size()-1] + 2);
                // 手动定义每个刻度的位置和显示的标签文本
                axisX->append("", iteration[0]-1);
                for(int i=0;i<order.size();i++)
                {
                    axisX->append(QString::number((iteration[i])), iteration[i]);
                }
                axisX->append(QString::number(iteration[iteration.size()-1]+ 1), iteration[iteration.size()-1]+ 1);
                axisX->append(QString::number(iteration[iteration.size()-1]+2), iteration[iteration.size()-1]+2);

                axisX->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
                chart->addAxis(axisX, Qt::AlignBottom);
                series->attachAxis(axisX);      // 曲线附加到这个 X 轴
                seriesReference->attachAxis(axisX);

                // 4. 创建并配置左侧 Y 轴（主 Y 轴）
                QLogValueAxis *axisYLeft = new QLogValueAxis;
                axisYLeft->setTitleText("Error");
                axisYLeft->setBase(10.0);     // ✅ 设置底数为10
                // axisYLeft->setRange(vfInputData.fittingError/10, 1);
                axisYLeft->setMin(vfInputData.fittingError/10); // 例如，从 0.1 开始
                axisYLeft->setMax(1);
                // axisYLeft->setTickCount(5);
                axisYLeft->setLabelFormat("%.0e");
                chart->addAxis(axisYLeft, Qt::AlignLeft);
                series->attachAxis(axisYLeft);  // 曲线附加到这个 Y 轴
                seriesReference->attachAxis(axisYLeft);  // 曲线附加到这个 Y 轴



                QCategoryAxis *axisYRight = new QCategoryAxis();
                axisYRight->setTitleText("Order");
                axisYRight->setRange(iteration[0]-1, iteration[iteration.size()-1]+2);

                // 手动定义每个刻度的位置和显示的标签文本
                axisYRight->append("", iteration[0]-1);
                for(int i=0;i<order.size();i++)
                {
                    axisYRight->append(QString::number((order[i])), iteration[i]);
                }
                axisYRight->append("", iteration[iteration.size()-1]+ 1);
                axisYRight->append("", iteration[iteration.size()-1]+2);

                axisYRight->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);

                // 将轴添加到图表中
                chart->addAxis(axisYRight, Qt::AlignTop);



                // 6. 显示图表
                QChartView *chartView = new QChartView(chart);
                chartView->setRenderHint(QPainter::Antialiasing);
                chartView->resize(800, 500);
                chartView->show();





                sparameterData.fittingError = QString::number(error);
                ui->label_fittingError->setText(sparameterData.fittingError);
                std::vector<std::vector<std::vector<std::complex<double>>>> VFData;
                getSparameter(sparameterData.frequency,vfInputData.poles1,vfInputData.residuse,vfInputData.constant,vfInputData.proportional,VFData);
                sparameterData.sparametrRIVF = VFData;


                std::vector<std::vector<std::vector<std::complex<double>>>> zparameter;
                s2z(sparameterData.sparametrRI,zparameter,sparameterData.referenceZ);
                sparameterData.zparametrRI = zparameter;
                s2z(sparameterData.sparametrRIVF,zparameter,sparameterData.referenceZ);
                sparameterData.zparametrRIVF = zparameter;

                std::vector<std::vector<std::vector<std::complex<double>>>> yparameter;
                s2y(sparameterData.sparametrRI,yparameter,sparameterData.referenceZ);
                sparameterData.yparametrRI = yparameter;
                s2y(sparameterData.yparametrRIVF,yparameter,sparameterData.referenceZ);
                sparameterData.yparametrRIVF = yparameter;

                sparameterData.poles = vfInputData.poles1;
                sparameterData.residuse = vfInputData.residuse;
                sparameterData.constant = vfInputData.constant;
                sparameterData.proportional = vfInputData.proportional;

                ui->comboBox_model->addItem(sparameterData.fileName,sparameterData.fileName);
                ui->comboBox_model->setCurrentText(sparameterData.fileName);

                qDebug()<<vfInputData.poles1.size();

                qDebug()<<"pole";
                for(const auto& value : sparameterData.poles)
                {
                    qDebug()<<value.real()<<value.imag();
                }
                qDebug()<<"residue";
                for(int i=0;i<sparameterData.residuse.size();i++)
                {
                    for(int j=0;j<sparameterData.residuse[i].size();j++)
                    {
                        qDebug()<<sparameterData.residuse[i][j].real()<<sparameterData.residuse[i][j].imag();
                    }
                }

                qDebug()<<sparameterData.constant;
                qDebug()<<sparameterData.proportional;   // 回主线程更新
            }, Qt::QueuedConnection);
        });


        m_watcher.setFuture(f);

        // 任务结束（无论完成还是取消）都会触发这个信号
        connect(&m_watcher, &QFutureWatcher<void>::finished,
                this, &MainWindow::onTaskFinished, Qt::UniqueConnection);

    }
}

void MainWindow::onTaskFinished() {
    m_cancelled.store(false);       // 复位，方便下次启动
}

void MainWindow::passivityTestData()
{
    QMessageBox msgBox;
    msgBox.setWindowTitle("Passivity Test Data");

    if(sparameterData.sparametrRI.empty())
    {
        msgBox.setText(QString("S parameter is empty."));
    }
    else
    {
        bool isPassiveOrigin = isPassive(sparameterData.sparametrRI);

        ui->plainTextEdit_log->appendPlainText(QString("Current sparameter passive: %1").arg(isPassiveOrigin ? "true" : "false"));

        if(isPassiveOrigin)
        {
            msgBox.setText(QString("Current sparameter is passive."));
        }
        else
        {
            msgBox.setText(QString("Current sparameter is not passive."));
        }


        std::vector<double>singular = getSingular(sparameterData.sparametrRI);

        StShowSParameterSingularDialog* dlg = new StShowSParameterSingularDialog(this);
        dlg->setGraphData(sparameterData.frequency, singular);
        // dlg.setWindowModality(Qt::ApplicationModal);
        dlg->show();
    }


    msgBox.setIcon(QMessageBox::Information);

    // 手动添加自定义按钮，并分配角色 (YesRole 或 NoRole)
    QPushButton *customYesButton = msgBox.addButton("Ok", QMessageBox::YesRole);

    // 设置默认按钮
    msgBox.setDefaultButton(customYesButton);

    // 显示对话框并获取用户操作结果
    msgBox.exec();
}

void MainWindow::passivityTestModel()
{
    QMessageBox msgBox;
    msgBox.setWindowTitle("Passivity Test Model");

    if(sparameterData.sparametrRIVF.empty())
    {
        msgBox.setText(QString("S parameter model is empty."));
    }
    else
    {
        bool isPassiveVF = isPassive(sparameterData.sparametrRIVF);

        ui->plainTextEdit_log->appendPlainText(QString("Current Model(sparameter) passive: %1").arg(isPassiveVF ? "true" : "false"));

        if(isPassiveVF)
        {
            msgBox.setText(QString("Current sparameter mdoel is passive."));
        }
        else
        {
            msgBox.setText(QString("Current sparameter model is not passive."));
        }


        std::vector<double>singular = getSingular(sparameterData.sparametrRIVF);

        StShowSParameterSingularDialog* dlg = new StShowSParameterSingularDialog(this);
        dlg->setGraphData(sparameterData.frequency, singular);
        // dlg.setWindowModality(Qt::ApplicationModal);
        dlg->show();
    }


    msgBox.setIcon(QMessageBox::Information);

    // 手动添加自定义按钮，并分配角色 (YesRole 或 NoRole)
    QPushButton *customYesButton = msgBox.addButton("Ok", QMessageBox::YesRole);

    // 设置默认按钮
    msgBox.setDefaultButton(customYesButton);

    // 显示对话框并获取用户操作结果
    msgBox.exec();
}

void MainWindow::causalityTestData()
{
    QMessageBox msgBox;
    msgBox.setWindowTitle("Causality Test Data");

    if(sparameterData.sparametrRI.empty())
    {
        msgBox.setText(QString("S parameter is empty."));
    }
    else
    {
        bool isCausalityOrigin = isCausal(sparameterData.sparametrRI);

        ui->plainTextEdit_log->appendPlainText(QString("Current sparameter causal: %1").arg(isCausalityOrigin ? "true" : "false"));

        if(isCausalityOrigin)
        {
            msgBox.setText(QString("Current sparameter is causal."));
        }
        else
        {
            msgBox.setText(QString("Current sparameter is not causal."));
        }
    }


    msgBox.setIcon(QMessageBox::Information);

    // 手动添加自定义按钮，并分配角色 (YesRole 或 NoRole)
    QPushButton *customYesButton = msgBox.addButton("Ok", QMessageBox::YesRole);

    // 设置默认按钮
    msgBox.setDefaultButton(customYesButton);

    // 显示对话框并获取用户操作结果
    msgBox.exec();
}

void MainWindow::causality()
{
    bool isCausality = isCausal(sparameterData.sparametrRI);

    ui->plainTextEdit_log->appendPlainText(QString("Current sparameter causality: %1").arg(isCausality ? "true" : "false"));
}

void MainWindow::passivityEnforceAction()
{
    std::vector<std::vector<std::complex<double>>> residueRes;
    std::vector<double> constantRes;
    passivityEnforceModel(sparameterData.portNum,sparameterData.frequency,
                          sparameterData.poles,
                          sparameterData.residuse,
                          sparameterData.proportional,
                          sparameterData.constant,sparameterData.proportional.empty()?0:1,sparameterData.constant.empty()?0:1,100,
                          residueRes,constantRes);

    sparameterData.residuse.clear();
    sparameterData.residuse = residueRes;
    sparameterData.constant.clear();
    sparameterData.constant = constantRes;
    std::vector<std::vector<std::vector<std::complex<double>>>> VFData;
    getSparameter(sparameterData.frequency,sparameterData.poles,sparameterData.residuse,sparameterData.constant,sparameterData.proportional,VFData);
    sparameterData.sparametrRIVF.clear();
    sparameterData.sparametrRIVF = VFData;

    std::vector<std::vector<std::vector<std::complex<double>>>> zparameter;
    s2z(sparameterData.sparametrRIVF,zparameter,sparameterData.referenceZ);
    sparameterData.zparametrRIVF = zparameter;

    std::vector<std::vector<std::vector<std::complex<double>>>> yparameter;
    s2y(sparameterData.yparametrRIVF,yparameter,sparameterData.referenceZ);
    sparameterData.yparametrRIVF = yparameter;

    passivityTestModel();

}

void MainWindow::logInformation(const QString &message)
{
    ui->plainTextEdit_log->appendPlainText(message);
}

void MainWindow::switchCurrentSnp()
{
    for(snp var : sparameterDataSet)
    {
        if(ui->comboBox_data->currentData().toString() == var.fileName)
        {
            sparameterData = var;
            initPortTable(sparameterData.portNum);
            initSparameterInformation();
            ui->label_fittingError->setText(sparameterData.fittingError);
            break;
        }
    }
}

void MainWindow::onDeleteData()
{
    for(int index=0;index< sparameterDataSet.size();index++)
    {
        if(sparameterData.fileName == sparameterDataSet[index].fileName)
        {
            sparameterDataSet.erase(sparameterDataSet.begin() + index);

            int dataIndex = ui->comboBox_data->findData(sparameterData.fileName);
            if(dataIndex != -1)
            {
                ui->comboBox_data->removeItem(dataIndex);
            }

            int modelIndex = ui->comboBox_model->findData(sparameterData.fileName);
            if(modelIndex != -1)
            {
                ui->comboBox_model->removeItem(modelIndex);
            }
        }
    }
    sparameterData.reset();

    if(!sparameterDataSet.empty())
    {
        sparameterData = sparameterDataSet[0];

        initPortTable(sparameterData.portNum);
        initSparameterInformation();
    }
    else
    {
        ui->plainTextEdit_log->clear();
        ui->plainTextEdit_dataInfo->clear();
        ui->plainTextEdit_modelInfo->clear();
        ui->comboBox_data->clear();
        ui->tableWidget_port->clear();
        ui->tableWidget_port->setRowCount(0);
        ui->tableWidget_port->setColumnCount(0);
    }

}

void MainWindow::onDeleteModel()
{
    for(int index=0;index< sparameterDataSet.size();index++)
    {
        if(sparameterData.fileName == sparameterDataSet[index].fileName)
        {
            sparameterDataSet[index].sparametrRIVF.clear();
            sparameterDataSet[index].zparametrRIVF.clear();

            int modelIndex = ui->comboBox_model->findData(sparameterData.fileName);
            if(modelIndex != -1)
            {
                ui->comboBox_model->removeItem(modelIndex);
            }
        }
    }
    sparameterData.sparametrRIVF.clear();
    sparameterData.yparametrRIVF.clear();
    sparameterData.zparametrRIVF.clear();

    ui->plainTextEdit_modelInfo->clear();
}

void MainWindow::onSetSpiceNetlist()
{
    StSettingSpiceSolverDialog* StSettingSpiceSolverDialogPtr = new StSettingSpiceSolverDialog(this);
    StSettingSpiceSolverDialogPtr->setNgspiceAndNetlistPath(m_spicePath,m_netlistPath);
    int result = StSettingSpiceSolverDialogPtr->exec();
    if(result == QDialog::Accepted)
    {
        m_spicePath = StSettingSpiceSolverDialogPtr->getNgspicePath();
        m_netlistPath = StSettingSpiceSolverDialogPtr->getNetlistPath();
    }
}

void MainWindow::onRunSpice(int portIndex)
{
    if(m_netlistPath.isEmpty())
    {
        return;
    }
    std::vector<std::vector<std::vector<std::complex<double>>>> sparameterRISpice = sparameterData.sparametrRINetlist;
    runSpiceFileSimulation(m_netlistPath, sparameterRISpice,portIndex);

    sparameterData.sparametrRINetlist = sparameterRISpice;
    qDebug() << sparameterRISpice.size();
}

void MainWindow::onGenerateSpiceNetlist()
{
    if(sparameterData.poles.empty())
    {
        qDebug() << "S model is empty!";
        return;
    }

    StGenerateNetlistDialog* stGenerateNetlistDialogPtr = new StGenerateNetlistDialog(this);

    int result = stGenerateNetlistDialogPtr->exec();

    if(result == QDialog::Accepted)
    {
        NetlistType netlistType = NetlistType::ControlledSource;
        stGenerateNetlistDialogPtr->getNetlistType(netlistType);

        // 选择单个图片文件
        QString filePath = QFileDialog::getSaveFileName(
            this,                     // 父窗口
            tr("选择文件"),           // 标题
            QDir::homePath(),        // 初始目录（用户主目录）
            tr("spice file(*.cir);;所有文件(*)")
            );

        QFileInfo fileInfo(filePath);
        QString fileNameWithoutExtension = fileInfo.baseName();  // 结果为 "file"

        if(filePath.isEmpty())
        {
            return;
        }

        if(netlistType == NetlistType::ControlledSource)
        {
            write_spice_subcircuit_s_controlled_source(sparameterData.portNum,
                                                       sparameterData.referenceZ,
                                                       sparameterData.poles,
                                                       sparameterData.residuse,
                                                       sparameterData.proportional,
                                                       sparameterData.constant,
                                                       filePath.toStdString(),
                                                       fileNameWithoutExtension.toStdString());
        }
        else if(netlistType == NetlistType::RLC)
        {
            write_spice_subcircuit_s_RLC(sparameterData.portNum,
                                         sparameterData.referenceZ,
                                         sparameterData.poles,
                                         sparameterData.residuse,
                                         sparameterData.proportional,
                                         sparameterData.constant,
                                         filePath.toStdString(),
                                         fileNameWithoutExtension.toStdString());
        }
        else
        {
            assert(true);
        }
    }
}

void MainWindow::onDeleteTempSpiceFile()
{
    QFile a(Path);
    // a.open(QIODevice::ReadOnly);
    a.remove();
}

void MainWindow::showModelInformation()
{
    if(sparameterData.poles.empty())
    {
        return;
    }
    ui->plainTextEdit_modelInfo->clear();

    ui->plainTextEdit_modelInfo->appendPlainText(QString("Port Num: %1").arg(sparameterData.portNum));
    ui->plainTextEdit_modelInfo->appendPlainText(QString("Reference Impedance: %1").arg(sparameterData.referenceZ));
    int modelOrder = 0;
    for(int i=0;i<sparameterData.poles.size();i++)
    {
        if(std::abs(sparameterData.poles[i].imag()) < equalZero)
        {
            modelOrder++;
        }
        else
        {
            modelOrder += 2;
        }
    }
    ui->plainTextEdit_modelInfo->appendPlainText(QString("Model Order: %1").arg(modelOrder));
    ui->plainTextEdit_modelInfo->appendPlainText(QString("Passive: %1").arg(isPassive(sparameterData.sparametrRIVF) ? "True":"False"));
    ui->plainTextEdit_modelInfo->appendPlainText(QString("Source: %1").arg(sparameterData.fileName));
}

void MainWindow::on_stop()
{
    m_cancelled.store(true);        // 通知工作线程退出
}

void MainWindow::open()
{
    sparameterData.filePath = "";
    sparameterData.frequency.clear();
    sparameterData.sparametrRI.clear();
    sparameterData.sparametrRIVF.clear();

    // 1. 弹出对话框，获取用户选择的文件路径
    sparameterFilePath = QFileDialog::getOpenFileName(
        this,                     // 父窗口
        "选择要打开的文件",       // 对话框标题
        "/home",                  // 初始目录（例如：“C:\\”、“/home”）
        "文本文件 (*.s*p);;所有文件 (*)" // 文件过滤器
        );


    // 检查用户是否取消了选择
    if (sparameterFilePath.isEmpty()) {
        return; // 用户点击了“取消”
    }

    readSparameter(sparameterFilePath,sparameterData);
    sparameterDataSet.push_back(sparameterData);

    initPortTable(sparameterData.portNum);

    initSparameterInformation();

    ui->comboBox_data->addItem(sparameterData.fileName,sparameterData.fileName);
    ui->comboBox_data->setCurrentText(sparameterData.fileName);
}

void MainWindow::onExportOriginSnpFile()
{
    if(sparameterData.sparametrRI.empty())
    {
        qDebug() << "S parameter is empty!";
        return;
    }

    // 选择单个图片文件
    QString filePath = QFileDialog::getSaveFileName(
        this,                     // 父窗口
        tr("选择图片"),           // 标题
        QDir::homePath(),        // 初始目录（用户主目录）
        tr("S parameter file(*.s*p);;所有文件(*)")
        );
    if (!filePath.isEmpty()) {
        writeSparameter(filePath, sparameterData, DATAORMODEL::DATA);
    }
}

void MainWindow::onExportVFSnpFile()
{
    if(sparameterData.sparametrRIVF.empty())
    {
        qDebug() << "Model S parameter is empty!";
        return;
    }

    // 选择单个图片文件
    QString filePath = QFileDialog::getSaveFileName(
        this,                     // 父窗口
        tr("选择图片"),           // 标题
        QDir::homePath(),        // 初始目录（用户主目录）
        tr("S parameter file(*.s*p);;所有文件(*)")
        );
    if (!filePath.isEmpty()) {
        writeSparameter(filePath, sparameterData, DATAORMODEL::MODEL);
    }
}

void MainWindow::showSparameter(int row, int col)
{
    onRunSpice(col);
    if(!ui->tableWidget_port->item(row,col))
    {
        return;
    }

    bool isLinearXAxis = false;
    if(ui->radioButton_xAxisLinear->isChecked())
    {
        isLinearXAxis = true;
    }
    else
    {
        isLinearXAxis = false;
    }

    bool isLinearYAxis = false;
    if(ui->radioButton_yAxisLinear->isChecked())
    {
        isLinearYAxis = true;
    }
    else
    {
        isLinearYAxis = false;
    }

    DisplayDataType dsiplayDataType = (DisplayDataType)ui->comboBox_DataType->currentData().toInt();

    DisplayParameterType displayParameterType = DisplayParameterType::S;
    if(ui->radioButton_Sparameter->isChecked())
    {
        displayParameterType = DisplayParameterType::S;
    }
    else if(ui->radioButton_Yparameter->isChecked())
    {
        displayParameterType = DisplayParameterType::Y;
    }
    else if(ui->radioButton_Zparameter->isChecked())
    {
        displayParameterType = DisplayParameterType::Z;
    }
    else
    {
        qDebug() << QString("No diaplay parameter type is checked!");
        return;
        // assert(false);
    }

    StShowSYZParameter* StShowSYZParameterPtr = new StShowSYZParameter(row,col,displayParameterType,
                                                                       isLinearXAxis,isLinearYAxis,
                                                                       sparameterData,dsiplayDataType,this);

    StShowSYZParameterPtr->show();
}
