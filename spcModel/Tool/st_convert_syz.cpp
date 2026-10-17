#include "st_convert_syz.h"
#include <Eigen/Dense>
#include <QDebug>
#include <QString>



void s2z(std::vector<std::vector<std::vector<std::complex<double> > > > &sparameterData,
         std::vector<std::vector<std::vector<std::complex<double> > > > &zparameterData,
         std::complex<double> referenceZ)

{

    int portNum = sparameterData[0].size();
    Eigen::MatrixXcd G(portNum,portNum);
    G.setZero();
    for(int i=0;i<portNum;i++)
    {
        G(i,i) = std::complex<double>(1/std::sqrt(std::abs(referenceZ.real())),0);
    }

    double UValue = std::sqrt(referenceZ.real())/std::abs(referenceZ);

    Eigen::MatrixXcd U(portNum,portNum);
    U.setZero();
    for(int i=0;i<portNum;i++)
    {
        U(i,i) = std::complex<double>(referenceZ.real(),0);
    }

    zparameterData.clear();
    zparameterData = sparameterData;
    Eigen::MatrixXcd eye = Eigen::MatrixXcd::Identity(portNum,portNum);
    for(int freq=0;freq<sparameterData.size();freq++)
    {
        Eigen::MatrixXcd SMatrix(portNum,portNum);
        for(int i=0;i<sparameterData[0].size();i++)
        {
            for(int j=0;j<sparameterData[0][0].size();j++)
            {
                SMatrix(i,j) = sparameterData[freq][i][j];
            }
        }


        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(eye - SMatrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Eigen::VectorXd& singular_values = svd.singularValues();
        double threshold = 1e-8; // 根据你的数值精度需求调整
        Eigen::MatrixXcd pinv_A = svd.matrixV() *
                                 (singular_values.array().abs() > threshold)
                                     .select(singular_values.array().inverse(), 0)
                                     .matrix()
                                     .asDiagonal() *
                                 svd.matrixU().adjoint();
        Eigen::MatrixXcd ZMatrix = referenceZ.real() * (eye + SMatrix) * pinv_A;
        for(int i=0;i<portNum;i++)
        {
            for(int j=0;j<portNum;j++)
            {
                zparameterData[freq][i][j] = ZMatrix(i,j);
            }
        }
    }
}




void s2y(const std::vector<std::vector<std::vector<std::complex<double> > > > &sParameterDataRI, std::vector<std::vector<std::vector<std::complex<double> > > > &yParameterDataRI, double referenceZ)
{
    if(sParameterDataRI.empty())
    {
        qDebug() << QString("S parameter is empty.");
        return ;
    }

    size_t portNum = sParameterDataRI[0].size();
    size_t frequencyNum = sParameterDataRI.size();
    yParameterDataRI.clear();

    Eigen::MatrixXcd I = Eigen::MatrixXcd::Identity(portNum,portNum);
    for(int freq=0;freq<frequencyNum;freq++)
    {
        Eigen::MatrixXcd SMatrix(portNum,portNum);
        for(int i=0;i<portNum;i++)
        {
            for(int j=0;j<portNum;j++)
            {
                SMatrix(i,j) = sParameterDataRI[freq][i][j];
            }
        }

        Eigen::MatrixXcd YMatrix = 1/referenceZ * (I + SMatrix).inverse() * (I - SMatrix);

        std::vector<std::vector<std::complex<double>>> yParameter(portNum,
                                                                  std::vector<std::complex<double>>(portNum,0));
        for(int i=0;i<portNum;i++)
        {
            for(int j=0;j<portNum;j++)
            {
                yParameter[i][j] = YMatrix(i,j);
            }
        }
        yParameterDataRI.push_back(yParameter);
    }
}

