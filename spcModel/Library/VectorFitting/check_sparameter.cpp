#include "check_sparameter.h"

#include <QDebug>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <QDebug>
#include <cmath>

bool isReciprocal(const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterDataRI, double threshold)
{
    if(sparameterDataRI.empty())
    {
        qDebug()<< "S parameter is empty.";
        return false;
    }

    size_t portNum = sparameterDataRI[0][0].size();
    size_t frequencyNum = sparameterDataRI.size();


    for(int k=0;k<frequencyNum;k++)
    {
        for(int i=0;i<portNum;i++)
        {
            for(int j=i;j<portNum;j++)
            {
                if(abs(sparameterDataRI[k][i][j] - sparameterDataRI[k][j][i]) > threshold)
                {
                    return false;
                }
            }
        }
    }

    return true;
}

bool isPassive(const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterDataRI, double threshold)
{
    if(sparameterDataRI.empty())
    {
        qDebug() << QString("S parameter is empty.");
        return false;
    }

    bool isPassive = true;
    size_t portNum = sparameterDataRI[0].size();
    size_t frequencyNum = sparameterDataRI.size();

    for(int freq=0;freq<frequencyNum;freq++)
    {
        Eigen::MatrixXcd SMatrix(portNum,portNum);
        for(int i=0;i<portNum;i++)
        {
            for(int j=0;j<portNum;j++)
            {
                SMatrix(i,j) = sparameterDataRI[freq][i][j];
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(SMatrix);
        Eigen::MatrixXd singularValues = svd.singularValues();
        double maxSingular = singularValues.maxCoeff();

        if(maxSingular > 1+threshold)
        {
            isPassive = false;
        }
    }

    return isPassive;
}

bool isCausal(const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterDataRI, double threshold)
{
    if(sparameterDataRI.empty())
    {
        qDebug() << QString("S parameter is empty.");
        return false;
    }


    double cqmi = calculate_cqmi_matrix(sparameterDataRI);

    if(cqmi > 80)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::vector<double> getSingular(const std::vector<std::vector<std::vector<std::complex<double> > > > &sparameterDataRI)
{
    if(sparameterDataRI.empty())
    {
        qDebug() << QString("S parameter is empty.");
        return std::vector<double>();
    }

    size_t portNum = sparameterDataRI[0].size();
    size_t frequencyNum = sparameterDataRI.size();
    std::vector<double> singularData(frequencyNum);

    for(int freq=0;freq<frequencyNum;freq++)
    {
        Eigen::MatrixXcd SMatrix(portNum,portNum);
        for(int i=0;i<portNum;i++)
        {
            for(int j=0;j<portNum;j++)
            {
                SMatrix(i,j) = sparameterDataRI[freq][i][j];
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(SMatrix);
        Eigen::MatrixXd singularValues = svd.singularValues();
        double maxSingular = singularValues.maxCoeff();
        singularData[freq] = maxSingular;
    }

    return singularData;
}


// 计算单个 S 参数序列的 CQMi
// s: 按频率递增排列的复数序列 S(f_k)
double calculate_cqmi_1port(const std::vector<std::complex<double>>& s) {
    const size_t n = s.size();
    if (n < 3) return 100.0;  // 数据点不足

    // 1. 相邻差分 ΔS_k = S(f_{k+1}) - S(f_k)
    std::vector<std::complex<double>> delta(n - 1);
    for (size_t k = 0; k + 1 < n; ++k) {
        delta[k] = s[k + 1] - s[k];
    }

    // 2. 计算相邻差分向量的叉积虚部 Cross_k = Im( conj(ΔS_k) * ΔS_{k+1} )
    double R_total = 0.0;      // 总旋转量
    double R_negative = 0.0;   // 负向（顺时针）旋转量

    for (size_t k = 0; k + 1 < delta.size(); ++k) {
        std::complex<double> cross = std::conj(delta[k]) * delta[k + 1];
        double cross_im = cross.imag();
        double abs_cross = std::abs(cross_im);
        R_total += abs_cross;
        if (cross_im < 0.0) {
            R_negative += abs_cross;
        }
    }

    if (R_total == 0.0) return 100.0;
    return (R_negative / R_total) * 100.0;
}

// 计算整个多端口 S 参数矩阵的 CQMi（取所有元素中的最小值）
// s_matrix[k][i][j]：第 k 个频率点，端口 i 到端口 j 的 S 参数
// 频率索引 k 从 0 到 nfreq-1，按频率递增排列
double calculate_cqmi_matrix(const std::vector<std::vector<std::vector<std::complex<double>>>>& s_matrix) {
    if (s_matrix.empty()) return 100.0;

    const size_t nfreq = s_matrix.size();
    const size_t nports = s_matrix[0].size();  // 假设每个频率点都是 nports x nports

    double min_cqmi = 100.0;

    for (size_t i = 0; i < nports; ++i) {
        for (size_t j = 0; j < nports; ++j) {
            // 提取 S_ij 随频率变化的序列
            std::vector<std::complex<double>> s_ij(nfreq);
            for (size_t k = 0; k < nfreq; ++k) {
                s_ij[k] = s_matrix[k][i][j];
            }

            double cqmi = calculate_cqmi_1port(s_ij);
            if (cqmi < min_cqmi) {
                min_cqmi = cqmi;
            }
        }
    }

    return min_cqmi;
}
