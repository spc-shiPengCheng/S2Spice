#include "st_file_io.h"

#include <QFile>
#include <QTextStream>
#include <QDebug>

#include "./st_struct_file.h"
#include "Tool/st_convert_syz.h"

void readSparameter(QString filePath,struct snp& data)
{
    // 2. 创建QFile对象并打开文件
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        // 打开失败，弹出错误信息
        return;
    }

    // 3. 使用QTextStream读取文本内容（适用于文本文件）
    QTextStream in(&file);
    // 设置编码，推荐UTF-8
    // in.setCodec("UTF-8");

    // QString fileContent = in.readAll(); // 一次性读取全部内容

    int portNum = filePath.at(filePath.size()-1-1).digitValue();
    QString sparameterFormat = "";
    std::vector<double> frequency;
    std::vector<std::vector<std::complex<double>>> sparameter = std::vector<std::vector<std::complex<double>>>(portNum,std::vector<std::complex<double>>(portNum,0));
    int rowCount  = 0;

    if(portNum == 1 || portNum == 2)
    {
        while (!in.atEnd()) { // 循环直到文件末尾
            QString lineOrigin = in.readLine(); // 核心：读取一行
            QString line = lineOrigin.simplified();

            rowCount ++;
            // 可以根据条件进行特定处理
            if (!line.startsWith("!")) {
                if(line.isEmpty())
                {
                    continue;
                }

                if(line.startsWith("#"))
                {
                    QStringList fileFormatData = line.simplified().split(" ");
                    data.referenceZ = fileFormatData[5].toDouble();
                    sparameterFormat = fileFormatData[3];

                    QString FrequencyUnit = fileFormatData[1];
                    if(!FrequencyUnit.compare("Hz",Qt::CaseInsensitive))
                    {
                        data.frequencyUnit = 1;
                        data.frequencyUnitLabel = "Hz";
                    }
                    else if(!FrequencyUnit.compare("KHz",Qt::CaseInsensitive))
                    {
                        data.frequencyUnit = 1e3;
                        data.frequencyUnitLabel = "KHz";
                    }
                    else if(!FrequencyUnit.compare("MHz",Qt::CaseInsensitive))
                    {
                        data.frequencyUnit = 1e6;
                        data.frequencyUnitLabel = "MHz";
                    }
                    else if(!FrequencyUnit.compare("GHz",Qt::CaseInsensitive))
                    {
                        data.frequencyUnit = 1e9;
                        data.frequencyUnitLabel = "GHz";
                    }
                    else
                    {
                        data.frequencyUnit = 1;
                        data.frequencyUnitLabel = "Hz";
                    }
                    // qDebug()<<fileFormatData;
                    continue;
                }

                if(sparameterFormat == "RI")
                {
                    QStringList sparameterData = line.split(" ");
                    frequency.push_back(sparameterData[0].toDouble()*data.frequencyUnit);

                    sparameterData.removeAt(0);
                    for(int i =0; i<sparameterData.size();)
                    {
                        sparameter[i/2/portNum][i/2%portNum] = std::complex<double>(sparameterData[i].toDouble(),sparameterData[i+1].toDouble());
                        i+=2;
                    }
                }
                else if(sparameterFormat == "DB")
                {
                    QStringList sparameterData = line.split(" ");
                    if(sparameterData.size() != 9)
                    {
                        qDebug() << rowCount;
                    }
                    frequency.push_back(sparameterData[0].toDouble()*data.frequencyUnit);

                    sparameterData.removeAt(0);
                    for(int i =0; i<sparameterData.size();)
                    {
                        double mag = std::pow(10,sparameterData[i].toDouble()/20);
                        double angle = sparameterData[i+1].toDouble()*(M_PI/180);
                        sparameter[i/2/portNum][i/2%portNum] = std::complex<double>(mag*std::cos(angle),mag*std::sin(angle));
                        i+=2;
                    }
                }
                else if(sparameterFormat == "MA")
                {
                    QStringList sparameterData = line.split(" ");
                    sparameterData.removeAll("");
                    if(sparameterData.size() != 9)
                    {
                        qDebug() << rowCount;
                    }
                    frequency.push_back(sparameterData[0].toDouble()*data.frequencyUnit);

                    sparameterData.removeAt(0);
                    for(int i =0; i<sparameterData.size();)
                    {
                        double mag = sparameterData[i].toDouble();
                        double angle = sparameterData[i+1].toDouble()*(M_PI/180);
                        sparameter[i/2/portNum][i/2%portNum] = std::complex<double>(mag*std::cos(angle),mag*std::sin(angle));
                        i+=2;
                    }
                }
                else
                {
                    assert(false);
                    return;
                }

                data.sparametrRI.push_back(sparameter);
            }
        }

    }
    else
    {
        while (!in.atEnd()) { // 循环直到文件末尾
            QString lineOrigin = in.readLine(); // 核心：读取一行
            lineOrigin.replace('\t', " ");
            QString line = lineOrigin.simplified();


            if (!line.startsWith("!")) {
                if(line.isEmpty())
                {
                    continue;
                }

                if(line.startsWith("#"))
                {
                    QStringList fileFormatData = line.split(" ");
                    data.referenceZ = fileFormatData[5].toDouble();
                    sparameterFormat = fileFormatData[3];

                    QString FrequencyUnit = fileFormatData[1];
                    if(!FrequencyUnit.compare("Hz",Qt::CaseInsensitive))
                    {
                        data.frequencyUnit = 1;
                        data.frequencyUnitLabel = "Hz";
                    }
                    else if(!FrequencyUnit.compare("KHz",Qt::CaseInsensitive))
                    {
                        data.frequencyUnit = 1e3;
                        data.frequencyUnitLabel = "KHz";
                    }
                    else if(!FrequencyUnit.compare("MHz",Qt::CaseInsensitive))
                    {
                        data.frequencyUnit = 1e6;
                        data.frequencyUnitLabel = "MHz";
                    }
                    else if(!FrequencyUnit.compare("GHz",Qt::CaseInsensitive))
                    {
                        data.frequencyUnit = 1e9;
                        data.frequencyUnitLabel = "GHz";
                    }
                    else
                    {
                        data.frequencyUnit = 1;
                        data.frequencyUnitLabel = "Hz";
                    }
                    // qDebug()<<fileFormatData;
                    continue;
                }


                if(sparameterFormat == "RI")
                {
                    QStringList sparameterData = line.split(" ");
                    frequency.push_back(sparameterData[0].toDouble()*data.frequencyUnit);

                    sparameterData.removeAt(0);

                    int readRowNum = (qCeil(portNum*2/std::min((double)sparameterData.size(),8.0)))*portNum;
                    for(int row = 1;row<readRowNum;row++)
                    {
                        QString lineOrigin = in.readLine(); // 核心：读取一行
                        QString line = lineOrigin.simplified();
                        line.replace('\t', " ");
                        if (!line.startsWith("!")) {
                            if(line.isEmpty())
                            {
                                continue;
                            }
                        }
                        QStringList sparameterDataNew = line.split(" ");
                        sparameterData.append(sparameterDataNew);
                    }

                    for(int i =0; i<sparameterData.size();)
                    {
                        sparameter[i/2/portNum][i/2%portNum] = std::complex<double>(sparameterData[i].toDouble(),sparameterData[i+1].toDouble());
                        i+=2;
                    }
                }
                else if(sparameterFormat == "DB")
                {
                    QStringList sparameterData = line.split(" ");
                    if(sparameterData.size() != 9)
                    {
                        qDebug() << rowCount;
                    }
                    frequency.push_back(sparameterData[0].toDouble()*data.frequencyUnit);

                    sparameterData.removeAt(0);

                    int readRowNum = (qCeil(portNum*2/std::min((double)sparameterData.size(),8.0)))*portNum;
                    for(int row = 1;row<readRowNum;row++)
                    {
                        QString lineOrigin = in.readLine(); // 核心：读取一行
                        QString line = lineOrigin.simplified();
                        line.replace('\t', " ");
                        if (!line.startsWith("!")) {
                            if(line.isEmpty())
                            {
                                continue;
                            }
                        }
                        QStringList sparameterDataNew = line.split(" ");
                        sparameterData.append(sparameterDataNew);

                        for(int i =0; i<sparameterData.size();)
                        {
                            double mag = std::pow(10,sparameterData[i].toDouble()/20);
                            double angle = sparameterData[i+1].toDouble()*(M_PI/180);
                            sparameter[i/2/portNum][i/2%portNum] = std::complex<double>(mag*std::cos(angle),mag*std::sin(angle));
                            i+=2;
                        }
                    }
                }
                else if(sparameterFormat == "MA")
                {
                    QStringList sparameterData = line.split(" ");
                    sparameterData.removeAll("");
                    if(sparameterData.size() != 9)
                    {
                        qDebug() << rowCount;
                    }
                    frequency.push_back(sparameterData[0].toDouble()*data.frequencyUnit);

                    sparameterData.removeAt(0);

                    int readRowNum = (qCeil(portNum*2/std::min((double)sparameterData.size(),8.0)))*portNum;
                    for(int row = 1;row<readRowNum;row++)
                    {
                        QString lineOrigin = in.readLine(); // 核心：读取一行
                        QString line = lineOrigin.simplified();
                        line.replace('\t', " ");
                        if (!line.startsWith("!")) {
                            if(line.isEmpty())
                            {
                                continue;
                            }
                        }
                        QStringList sparameterDataNew = line.split(" ");
                        sparameterData.append(sparameterDataNew);

                        for(int i =0; i<sparameterData.size();)
                        {
                            double mag = sparameterData[i].toDouble();
                            double angle = sparameterData[i+1].toDouble()*(M_PI/180);
                            sparameter[i/2/portNum][i/2%portNum] = std::complex<double>(mag*std::cos(angle),mag*std::sin(angle));
                            i+=2;
                        }
                    }
                }
                else
                {
                    assert(false);
                    return;
                }

                data.sparametrRI.push_back(sparameter);
            }
        }
    }

    data.frequency = frequency;
    data.filePath = filePath;
    QStringList temp = filePath.split("/");
    data.fileName = temp[temp.size()-1].split(".")[0];
    data.portNum = portNum;
    data.fileSize = file.size();

    file.close(); // 操作完成后关闭文件

    std::vector<std::vector<std::vector<std::complex<double>>>> zparameter;
    s2z(data.sparametrRI,zparameter,data.referenceZ);
    data.zparametrRI = zparameter;
    std::vector<std::vector<std::vector<std::complex<double>>>> yparameter;
    s2y(data.sparametrRI,yparameter,data.referenceZ);
    data.yparametrRI = yparameter;
}

void writeSparameter(QString filePath, snp &data, DATAORMODEL dataOrModel)
{
    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qWarning() << "无法打开文件：" << file.errorString();
        return;
    }

    QTextStream out(&file);
    QString lineData = "";
    if(data.portNum == 1 || data.portNum == 2)
    {
        out << "#" << " " << data.frequencyUnitLabel << " " << "S" << " "<< "RI"<< " " <<"R"<< " " << data.referenceZ << "\n";// 自动处理编码，默认 UTF-8
        size_t frequencyCount = data.frequency.size();

        if(dataOrModel == DATAORMODEL::DATA)
        {
            for(int i=0;i<frequencyCount;i++)
            {
                lineData = "";
                lineData += QString::number(data.frequency[i]/data.frequencyUnit) + " ";
                for(int j=0;j<data.portNum;j++)
                {
                    for(int k=0;k<data.portNum;k++)
                    {
                        lineData += QString::number(data.sparametrRI[i][j][k].real()) + " " + QString::number(data.sparametrRI[i][j][k].imag()) + " ";
                    }
                }
                lineData.chop(1);
                out << lineData << "\n";
            }
        }
        else if(dataOrModel == DATAORMODEL::MODEL)
        {
            for(int i=0;i<frequencyCount;i++)
            {
                lineData = "";
                lineData += QString::number(data.frequency[i]/data.frequencyUnit) + " ";
                for(int j=0;j<data.portNum;j++)
                {
                    for(int k=0;k<data.portNum;k++)
                    {
                        lineData += QString::number(data.sparametrRIVF[i][j][k].real()) + " " + QString::number(data.sparametrRIVF[i][j][k].imag()) + " ";
                    }
                }
                lineData.chop(1);
                out << lineData << "\n";
            }
        }

        out.flush();             // 确保数据写入磁盘
    }
    else
    {
        out << "#" << " " << data.frequencyUnitLabel << " " << "S" << " "<< "RI"<< " " <<"R"<< " " << data.referenceZ << "\n";// 自动处理编码，默认 UTF-8
        size_t frequencyCount = data.frequency.size();
        int lineCounter = 0;
        int freqLineCounter =0;

        if(dataOrModel == DATAORMODEL::DATA)
        {
            for(int i=0;i<frequencyCount;i++)
            {
                lineData = "";
                lineCounter = 0;
                freqLineCounter = 0;
                lineData += QString::number(data.frequency[i]/data.frequencyUnit) + " ";
                for(int j=0;j<data.portNum;j++)
                {
                    for(int k=0;k<data.portNum;k++)
                    {
                        lineCounter++;
                        freqLineCounter++;
                        lineData += QString::number(data.sparametrRI[i][j][k].real()) + " " + QString::number(data.sparametrRI[i][j][k].imag()) + " ";
                        if(lineCounter == 4 || freqLineCounter == data.portNum)
                        {
                            if(freqLineCounter == data.portNum)
                            {
                                freqLineCounter = 0;
                            }
                            lineCounter = 0;
                            lineData.chop(1);
                            out << lineData << "\n";
                            lineData = "";
                        }
                    }
                }
                lineData.chop(1);
                out << lineData << "\n";
            }
        }
        else if(dataOrModel == DATAORMODEL::MODEL)
        {
            for(int i=0;i<frequencyCount;i++)
            {
                lineData = "";
                lineCounter = 0;
                freqLineCounter = 0;
                lineData += QString::number(data.frequency[i]/data.frequencyUnit) + " ";
                for(int j=0;j<data.portNum;j++)
                {
                    for(int k=0;k<data.portNum;k++)
                    {
                        lineCounter++;
                        freqLineCounter++;
                        lineData += QString::number(data.sparametrRIVF[i][j][k].real()) + " " + QString::number(data.sparametrRIVF[i][j][k].imag()) + " ";
                        if(lineCounter == 4 || freqLineCounter == data.portNum)
                        {
                            if(freqLineCounter == data.portNum)
                            {
                                freqLineCounter = 0;
                            }
                            lineCounter = 0;
                            lineData.chop(1);
                            out << lineData << "\n";
                            lineData = "";
                        }
                    }
                }
                lineData.chop(1);
                out << lineData << "\n";
            }
        }

        out.flush();             // 确保数据写入磁盘
    }


    file.close();

    if (out.status() != QTextStream::Ok) {
        qWarning() << "写入过程中出错";
        return;
    }
}
