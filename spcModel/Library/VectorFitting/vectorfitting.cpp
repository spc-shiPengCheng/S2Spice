#include "vectorfitting.h"

#include <Eigen/Dense>
#include <QString>
#include <unsupported/Eigen/CXX11/Tensor>
#include <QDebug>
#include <omp.h>
// #define EIGEN_STACK_ALLOCATION_LIMIT 0 // 0表示无限制

#include "../../Tool/st_struct_file.h"

#include "../../Tool/st_logger.h"




double VF(int realPolesNum,int complexPolesNum,bool isProportional,bool isConstant,
        int portNum, std::vector<double> frequency,
        std::vector<std::vector<std::vector<std::complex<double> > > > sparameterData,
        int iterations, double rms,bool isRelaxed,
        std::vector<std::complex<double> > &poles1,
        std::vector<std::vector<std::complex<double> > > &residuse,
        std::vector<double> &constant,
        std::vector<double>& proportional)
{
    if(frequency.empty() ||
       sparameterData.empty())
    {
        qDebug() << QString("Input VF data is empty!");
        return 0;
    }
    int frequencySize = frequency.size();
    std::vector<std::vector<std::complex<double>>> freqResponse(frequencySize,
                                                                std::vector<std::complex<double>>(portNum*portNum,0));

    double freqAverage = 0;
    for(int freq = 0;freq<frequencySize;freq++)
    {
        freqAverage += frequency[freq];
        for(int row = 0;row<portNum;row++)
        {
            for(int col=0;col<portNum;col++)
            {
                freqResponse[freq][row*portNum+col] = sparameterData[freq][row][col];
            }
        }
    }
    freqAverage /= frequencySize;


    double fMin = frequency[0];
    double fMax = frequency[frequencySize-1];

    double interval = (fMax - fMin)/(realPolesNum-1);
    std::vector<std::complex<double>> init_poles(realPolesNum + complexPolesNum,0);
    for(int i=0;i<realPolesNum;i++)
    {
        init_poles[i] = std::complex<double>(-2*M_PI*(fMin+i*interval),0);
    }

    if(complexPolesNum == 1)
    {
        init_poles[realPolesNum] = std::complex<double>(-0.01*2*M_PI*(fMin+fMax)/2,2*M_PI*(fMin+fMax)/2);
    }
    else
    {
        interval = (fMax - fMin)/(complexPolesNum-1);
        for(int i=0;i<complexPolesNum;i++)
        {
            init_poles[i + realPolesNum] = std::complex<double>(-0.01*2*M_PI*(fMin+i*interval),2*M_PI*(fMin+i*interval));
        }
    }


    for(int i=0;i<init_poles.size();i++)
    {
        init_poles[i] = init_poles[i]/freqAverage;
        qDebug()<<QString("Pole real: %1").arg(init_poles[i].real());
        qDebug()<<QString("Pole imag: %1").arg(init_poles[i].imag());
    }


    std::vector<std::complex<double>> sFreq(frequencySize,0);
    for(int i=0;i<frequencySize;i++)
    {
        sFreq[i] = std::complex<double>(0,2*M_PI*frequency[i])/freqAverage;
        // sFreq[i] = std::complex<double>(0,2*M_PI*frequency[i]);
    }


    double finalRms = 0;
    while(iterations>0)
    {
        int colOrigin = 0;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                colOrigin += 1;
            }
            else
            {
                colOrigin += 2;
            }
        }

        int colSigma = 0;
        //relaxed
        if(isRelaxed)
        {
            colSigma = colOrigin + 1;
        }
        else
        {
            // colSigma = colOrigin;
            colSigma = colOrigin + 1;
        }


        int indexConstant = colOrigin;
        if(isConstant)
        {
            colOrigin++;
        }

        int indexProportional = colOrigin;
        if(isProportional)
        {
            colOrigin++;
        }

        std::vector<int> indexRealPolse;
        std::vector<int> indexComplexPoles;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                indexRealPolse.push_back(i);
            }
            else
            {
                indexComplexPoles.push_back(i);
            }
        }

        int realPolesNum = indexRealPolse.size();
        int complexPolesNum = indexComplexPoles.size();

        std::vector<int> idx_A_real(realPolesNum);
        for(int i=0;i<realPolesNum;i++)
        {
            idx_A_real[i] = i;
        }
        std::vector<int> idx_A_complex_re(complexPolesNum);
        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_re[i] = realPolesNum + 2*i;
        }
        std::vector<int> idx_A_complex_im(complexPolesNum);
        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_im[i] = idx_A_complex_re[i] + 1;
        }

        // 假设 s, poles 是 Eigen 列向量，index_real_poles 是 std::vector<int> 或 Eigen::VectorXi
        Eigen::VectorXcd s(sFreq.size()); // 可能是复数，例如在频域分析中
        for(int i=0;i<frequencySize;i++)
        {
            s[i]= sFreq[i];
        }

        Eigen::VectorXcd poles(init_poles.size());
        for(int i=0;i<init_poles.size();i++)
        {
            poles[i]= init_poles[i];
        }

        Eigen::VectorXi index_real_poles(indexRealPolse.size()); // 或 std::vector<int>
        for(int i=0;i<realPolesNum;i++)
        {
            index_real_poles[i]= indexRealPolse[i];
        }

        // 1. 从 poles 中选出指定的极点，并转换为行向量
        Eigen::RowVectorXcd selected_poles = poles(index_real_poles).transpose();
        // poles(index_real_poles) 是一个列向量，.transpose() 将其转为行向量

        // 2. 获取维度
        int m = s.size();
        int k = selected_poles.size(); // 也等于 index_real_poles.size()

        // 3. 通过复制扩展创建两个 (m, k) 的矩阵
        // s 作为列向量，将其复制扩展为 m x k 矩阵，每一列都是 s
        Eigen::MatrixXcd s_expanded = s.replicate(1, k);
        // selected_poles 作为行向量，将其复制扩展为 m x k 矩阵，每一行都是 selected_poles
        Eigen::MatrixXcd p_expanded = selected_poles.replicate(m, 1);

        // 4. 执行逐元素减法，得到差值矩阵
        Eigen::MatrixXcd diff_matrix = s_expanded - p_expanded;
        // 此时 diff_matrix(i, j) = s(i) - selected_poles(j)

        // 5. 对差值矩阵的每个元素求倒数（使用 .array() 进行逐元素操作）
        Eigen::MatrixXcd coeff_real = diff_matrix.array().inverse();





        Eigen::VectorXi index_complex_poles(complexPolesNum); // 复数极点的索引，k x 1
        for(int i=0;i<complexPolesNum;i++)
        {
            index_complex_poles[i]= indexComplexPoles[i];
        }

        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cR = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj = selected_poles_cR.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cR.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cR = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cR = selected_poles_cR.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded = selected_poles_conj.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数
        Eigen::MatrixXcd term1_matrix = (s_expanded_cR - p_expanded_cR).array().inverse();
        Eigen::MatrixXcd term2_matrix = (s_expanded_cR - p_conj_expanded).array().inverse();

        // 7. 两个项相加得到最终结果
        Eigen::MatrixXcd coeff_complex_re = term1_matrix + term2_matrix;






        // 定义虚数单位 1j（Python中为1j，C++中为complex(0,1)）
        std::complex<double> I(0.0, 1.0);

        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cI = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj_I = selected_poles_cI.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cI.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cI = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cI = selected_poles_cI.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded_cI = selected_poles_conj_I.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数，然后乘以虚数单位I
        Eigen::MatrixXcd term1_matrix_cI = I * (s_expanded_cI - p_expanded_cI).array().inverse();
        Eigen::MatrixXcd term2_matrix_cI = I * (s_expanded_cI - p_conj_expanded_cI).array().inverse();

        // 7. 两个项相减得到最终结果（term1 - term2）
        Eigen::MatrixXcd coeff_complex_im = term1_matrix_cI - term2_matrix_cI;




        Eigen::Tensor<std::complex<double>, 3> A(portNum*portNum, frequencySize, colOrigin+colSigma);
        A.setZero();


        // coeff_real 是二维矩阵（需要转换为Tensor）
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> coeff_real_tensor(
            coeff_real.data(), coeff_real.rows(), coeff_real.cols());

        for (int index = 0;index<idx_A_real.size();index++) {
            int idx = idx_A_real[index];
            if (idx >= 0 && idx < A.dimension(2)) {
                // 对每个第三维索引进行赋值
                for (int i = 0; i < A.dimension(0); i++) {
                    for (int j = 0; j < A.dimension(1); j++) {
                        A(i, j, idx) = coeff_real(j, index);
                    }
                }
            }
        }


        // coeff_real 是二维矩阵（需要转换为Tensor）
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> coeff_complexReal_tensor(
            coeff_complex_re.data(), coeff_complex_re.rows(), coeff_complex_re.cols());

        for (int index = 0;index<idx_A_complex_re.size();index++) {
            int idx = idx_A_complex_re[index];
            if (idx >= 0 && idx < A.dimension(2)) {
                // 对每个第三维索引进行赋值
                for (int i = 0; i < A.dimension(0); i++) {
                    for (int j = 0; j < A.dimension(1); j++) {
                        A(i, j, idx) = coeff_complexReal_tensor(j, index);
                    }
                }
            }
        }



        // coeff_real 是二维矩阵（需要转换为Tensor）
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> coeff_complexImag_tensor(
            coeff_complex_im.data(), coeff_complex_im.rows(), coeff_complex_im.cols());


        for (int index = 0;index<idx_A_complex_im.size();index++) {
            int idx = idx_A_complex_im[index];
            if (idx >= 0 && idx < A.dimension(2)) {
                // 对每个第三维索引进行赋值
                for (int i = 0; i < A.dimension(0); i++) {
                    for (int j = 0; j < A.dimension(1); j++) {
                        A(i, j, idx) = coeff_complexImag_tensor(j, index);
                    }
                }
            }
        }


        if(isConstant)
        {
            for(int i=0;i<A.dimension(0);i++)
            {
                for(int j=0;j<A.dimension(1);j++)
                {
                    A(i,j,indexConstant) = std::complex<double>(1,0);
                }
            }
        }

        if(isProportional)
        {
            for(int i=0;i<A.dimension(0);i++)
            {
                for(int j=0;j<A.dimension(1);j++)
                {
                    A(i,j,indexProportional) = sFreq[j];
                }
            }
        }


        for(int i=0;i<A.dimension(0);i++)
        {
            for(int j=0;j<A.dimension(1);j++)
            {
                for(int k=0;k<(int)idx_A_real.size();k++)
                {
                    A(i,j,idx_A_real[k]+colOrigin) = std::complex<double>(-1,0)*freqResponse[j][i]*coeff_real(j,k);
                }

                for(int k=0;k<(int)idx_A_complex_re.size();k++)
                {
                    A(i,j,idx_A_complex_re[k]+colOrigin) = std::complex<double>(-1,0)*freqResponse[j][i]*coeff_complex_re(j,k);
                }

                for(int k=0;k<idx_A_complex_im.size();k++)
                {
                    A(i,j,idx_A_complex_im[k]+colOrigin) = std::complex<double>(-1,0)*freqResponse[j][i]*coeff_complex_im(j,k);
                }

                // if(isRelaxed)
                {
                    A(i,j,A.dimension(2)-1) = std::complex<double>(-1,0)*freqResponse[j][i];
                }
            }
        }

        // A.reshape(portNum*portNum*frequencySize,colOrigin+colSigma);

        Eigen::Map<Eigen::MatrixXcd> AReCombine(A.data(), portNum*portNum*frequencySize,colOrigin+colSigma);

        // Eigen::MatrixXcd AReCombine;

        // 1. 获取矩阵维度
        int rows = AReCombine.rows();
        int cols = AReCombine.cols();

        int AReCombineRealImagRows = 0;
        if(isRelaxed)
        {
            AReCombineRealImagRows = 2 * rows+1;
        }
        else
        {
            AReCombineRealImagRows = 2 * rows;
        }

        // 2. 创建结果矩阵（行数为原矩阵2倍，列数相同）
        Eigen::MatrixXd AReCombineRealImag(AReCombineRealImagRows, cols);

        // 3. 将实部放在上半部分
        // AReCombineRealImag.topRows(rows) = AReCombine.real();
        AReCombineRealImag.block(0,0,rows,cols) = AReCombine.real();

        // 4. 将虚部放在下半部分
        // AReCombineRealImag.bottomRows(rows) = AReCombine.imag();
        AReCombineRealImag.block(rows,0,rows,cols) = AReCombine.imag();

        if(isRelaxed)
        {
            Eigen::RowVectorXd row1(cols);
            row1.setConstant(frequencySize);

            AReCombineRealImag.row(2*rows) = row1;
        }



        int bRows = 0;
        if(isRelaxed)
        {
            bRows = 2*portNum*portNum*frequencySize + 1;
        }
        else
        {
            bRows = 2*portNum*portNum*frequencySize;
        }
        Eigen::VectorXd b(bRows);
        b.setZero();
        if(isRelaxed)
        {
            double sumFrequencyReal = 0;
            double freqResponseNorm2 = 0;
            for(int i=0;i<frequencySize;i++)
            {
                for(int j=0;j<portNum;j++)
                {
                    for(int k=0;k<portNum;k++)
                    {
                        sumFrequencyReal += sparameterData[i][j][k].real();

                        freqResponseNorm2 += std::pow(std::abs(sparameterData[i][j][k]),2);
                    }
                }
            }
            b(bRows - 1) = sumFrequencyReal;

            freqResponseNorm2 = std::sqrt(freqResponseNorm2);

            AReCombineRealImag.row(2*rows) = AReCombineRealImag.row(2*rows) * freqResponseNorm2;
            b(bRows - 1) *= freqResponseNorm2;
        }


        Eigen::VectorXd x;
        if(isRelaxed)
        {
            x = AReCombineRealImag.colPivHouseholderQr().solve(b);
        }
        else
        {
            // x = getNonZeroSolutionQR(AReCombineRealImag);

            b = AReCombineRealImag.col(AReCombineRealImag.cols()-1);
            Eigen::MatrixXd dealA = AReCombineRealImag.block(0,0,AReCombineRealImag.rows(),AReCombineRealImag.cols()-1);
            x = dealA.colPivHouseholderQr().solve(b);
        }


        std::vector<double> c_res;
        for(int i=colOrigin;i<x.size();i++)
        {
            c_res.push_back(x[i]);
        }

        double d_res = 0;
        if(isRelaxed)
        {
            d_res = x[x.size()-1];
            if(std::abs(d_res)<1e-8)
            {
                if(!d_res)
                {
                    d_res = 1e-8;
                }
                else
                {
                    d_res = 1e-8*std::abs(d_res)/d_res;
                }
            }
        }





        Eigen::MatrixXcd H(c_res.size(),c_res.size());
        H.setZero();

        std::vector<std::complex<double>> realPoles;
        std::vector<std::complex<double>> complexPoles;

        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                realPoles.push_back(init_poles[i]);
            }
            else
            {
                complexPoles.push_back(init_poles[i]);
            }
        }

        for(int i=0;i<idx_A_real.size();i++)
        {
            H(idx_A_real[i],idx_A_real[i]) = realPoles[i].real();
        }

        for(int i=0;i<idx_A_real.size();i++)
        {
            for(int j=0;j<H.cols();j++)
            {
                if(isRelaxed)
                {
                    H(idx_A_real[i],j) -= c_res[j]/d_res;
                }
                else
                {
                     H(idx_A_real[i],j) -= c_res[j];
                }
            }
        }

        for(int i=0;i<(int)idx_A_complex_re.size();i++)
        {
            for(int j=0;j<idx_A_complex_re.size();j++)
            {
                H(idx_A_complex_re[i],idx_A_complex_re[j]) = complexPoles[j].real();
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                H(idx_A_complex_re[i],idx_A_complex_im[j]) = complexPoles[j].imag();
            }
        }

        for(int i=0;i<idx_A_complex_im.size();i++)
        {
            for(int j=0;j<idx_A_complex_re.size();j++)
            {
                H(idx_A_complex_im[i],idx_A_complex_re[j]) = -1*complexPoles[j].imag();
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                H(idx_A_complex_im[i],idx_A_complex_im[j]) = complexPoles[j].real();
            }
        }

        for(int i=0;i<idx_A_complex_re.size();i++)
        {
            for(int j=0;j<H.cols();j++)
            {
                // H(idx_A_complex_re[i],j) -= 2*c_res[j]/d_res;
                if(isRelaxed)
                {
                     H(idx_A_complex_re[i],j) -= c_res[j]/d_res;
                }
                else
                {
                    H(idx_A_complex_re[i],j) -= c_res[j];
                }
            }
        }

        Eigen::VectorXcd polesNew = H.eigenvalues();

        for(int i=0;i<polesNew.size();i++)
        {
            if(polesNew(i).real()>0)
            {
                polesNew(i) = -polesNew(i);
            }
        }

        init_poles.clear();
        for(int i=0;i<polesNew.size();i++)
        {
            if(std::abs(polesNew[i].real())<equalZero &&
                std::abs(polesNew[i].imag())<equalZero)
            {

            }
            else
            {
                init_poles.push_back(polesNew(i));
                // qDebug()<< polesNew(i).real() << polesNew(i).imag();
            }
        }




        int n_cols = 0;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                n_cols +=1;
            }
            else
            {
                n_cols+=2;
            }
        }



        std::vector<int> realIndex;
        std::vector<int> complexIndex;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                // index_real_poles << i;
                realIndex.push_back(i);
            }
            else
            {
                // index_complex_poles << i;
                complexIndex.push_back(i);
            }
        }

        index_real_poles.resize(realIndex.size());
        index_complex_poles.resize(complexIndex.size());
        for(int i=0;i<index_real_poles.size();i++)
        {
            index_real_poles(i) = realIndex[i];
        }
        for(int i=0;i<index_complex_poles.size();i++)
        {
            index_complex_poles(i) = complexIndex[i];
        }

        if(isConstant)
        {
            indexConstant = n_cols;
            n_cols++;
        }

        if(isProportional)
        {
            indexProportional = n_cols;
            n_cols++;
        }


        int index = 0;
        idx_A_real.clear();
        idx_A_complex_re.clear();
        idx_A_complex_im.clear();
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                idx_A_real.push_back(index);
                index++;
            }
            else
            {
                idx_A_complex_re.push_back(index);
                idx_A_complex_im.push_back(index+1);
                index+=2;
            }
        }

        Eigen::MatrixXcd AResidues(frequencySize,n_cols);


        poles.resize(init_poles.size());
        for(int i=0;i<init_poles.size();i++)
        {
            poles[i]= init_poles[i];
        }

        // 1. 从 poles 中选出指定的极点，并转换为行向量
        Eigen::RowVectorXcd selected_poles_res = poles(index_real_poles).transpose();
        // poles(index_real_poles) 是一个列向量，.transpose() 将其转为行向量

        // 2. 获取维度
        m = s.size();
        k = selected_poles_res.size(); // 也等于 index_real_poles.size()

        // 3. 通过复制扩展创建两个 (m, k) 的矩阵
        // s 作为列向量，将其复制扩展为 m x k 矩阵，每一列都是 s
        Eigen::MatrixXcd s_expanded_res = s.replicate(1, k);
        // selected_poles 作为行向量，将其复制扩展为 m x k 矩阵，每一行都是 selected_poles
        Eigen::MatrixXcd p_expanded_res = selected_poles_res.replicate(m, 1);

        // 4. 执行逐元素减法，得到差值矩阵
        Eigen::MatrixXcd diff_matrix_res = s_expanded_res - p_expanded_res;
        // 此时 diff_matrix(i, j) = s(i) - selected_poles_res(j)

        // 5. 对差值矩阵的每个元素求倒数（使用 .array() 进行逐元素操作）
        Eigen::MatrixXcd coeff_real_res = diff_matrix_res.array().inverse();





        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cR_res = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj_res = selected_poles_cR_res.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cR_res.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cR_res = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cR_res = selected_poles_cR_res.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded_res = selected_poles_conj_res.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数
        Eigen::MatrixXcd term1_matrix_res = (s_expanded_cR_res - p_expanded_cR_res).array().inverse();
        Eigen::MatrixXcd term2_matrix_res = (s_expanded_cR_res - p_conj_expanded_res).array().inverse();

        // 7. 两个项相加得到最终结果
        Eigen::MatrixXcd coeff_complex_re_res = term1_matrix_res + term2_matrix_res;







        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cI_res = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj_I_res = selected_poles_cI_res.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cI_res.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cI_res = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cI_res = selected_poles_cI_res.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded_cI_res = selected_poles_conj_I_res.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数，然后乘以虚数单位I
        Eigen::MatrixXcd term1_matrix_cI_res = I * (s_expanded_cI_res - p_expanded_cI_res).array().inverse();
        Eigen::MatrixXcd term2_matrix_cI_res = I * (s_expanded_cI_res - p_conj_expanded_cI_res).array().inverse();

        // 7. 两个项相减得到最终结果（term1 - term2）
        Eigen::MatrixXcd coeff_complex_im_res = term1_matrix_cI_res - term2_matrix_cI_res;



        for(int i=0;i<AResidues.rows();i++)
        {
            for(int j=0;j<idx_A_real.size();j++)
            {
                AResidues(i,idx_A_real[j]) = coeff_real_res(i,j);
            }

            for(int j=0;j<idx_A_complex_re.size();j++)
            {
                AResidues(i,idx_A_complex_re[j]) = coeff_complex_re_res(i,j);
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                AResidues(i,idx_A_complex_im[j]) = coeff_complex_im_res(i,j);
            }

            if(isConstant)
            {
                AResidues(i,indexConstant) = std::complex<double>(1,0);
            }

            if(isProportional)
            {
                AResidues(i,indexProportional) = sFreq[i];
            }
        }



        // 1. 获取矩阵维度
        rows = AResidues.rows();
        cols = AResidues.cols();

        // 2. 创建结果矩阵（行数为原矩阵2倍，列数相同）
        Eigen::MatrixXd AResiduesCombineRealImag(2 * rows, cols);

        // 3. 将实部放在上半部分
        AResiduesCombineRealImag.topRows(rows) = AResidues.real();

        // 4. 将虚部放在下半部分
        AResiduesCombineRealImag.bottomRows(rows) = AResidues.imag();


        Eigen::MatrixXd freqResponseRI(frequencySize*2,freqResponse[0].size());
        for(int i=0;i<freqResponse.size();i++)
        {
            for(int j=0;j<freqResponse[0].size();j++)
            {
                freqResponseRI(i,j) = freqResponse[i][j].real();
                freqResponseRI(i+frequencySize,j) = freqResponse[i][j].imag();
            }
        }

        Eigen::MatrixXd xResiduse = AResiduesCombineRealImag.colPivHouseholderQr().solve(freqResponseRI);


        Eigen::MatrixXcd residuseTemp(portNum*portNum,index_real_poles.size() + 2*index_complex_poles.size());

        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<(int)(idx_A_real.size());j++)
            {
                residuseTemp(i,index_real_poles(j)) = xResiduse(idx_A_real[j],i);
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                residuseTemp(i,index_complex_poles(j)) = std::complex<double>(xResiduse(idx_A_complex_re[j],i),xResiduse(idx_A_complex_im[j],i));
            }
        }

        if(isConstant)
        {
            Eigen::VectorXd constant_coeff = xResiduse.row(indexConstant);

            constant.clear();
            constant.resize(constant_coeff.size());
            for(int i=0;i<constant_coeff.size();i++)
            {
                constant[i] = constant_coeff(i);
            }
        }

        if(isProportional)
        {
            Eigen::VectorXd proportional_coeff = xResiduse.row(indexProportional);

            proportional.clear();
            proportional.resize(proportional_coeff.size());
            for(int i=0;i<proportional_coeff.size();i++)
            {
                // proportional[i] = proportional_coeff(i)/freqAverage;
                proportional[i] = proportional_coeff(i);
            }
        }



        residuse.clear();
        residuse = std::vector<std::vector<std::complex<double>>>(residuseTemp.rows(),std::vector<std::complex<double>>(residuseTemp.cols(),0));
        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<residuseTemp.cols();j++)
            {
                residuse[i][j] = residuseTemp(i,j)*freqAverage;
                // residuse[i][j] = residuseTemp(i,j);
            }
        }

        for(int i=0;i<init_poles.size();i++)
        {
            init_poles[i] *= freqAverage;
        }



        poles1.clear();
        poles1 = init_poles;



        qDebug()<<QString("The %1th Iteration.").arg(10 - iterations + 1);


        std::vector<std::vector<std::vector<std::complex<double>>>> sparameterDataVF;
        getSparameter(frequency,init_poles,residuse,constant,proportional,sparameterDataVF);

        double currentRms = getRMS(sparameterData,sparameterDataVF);
        finalRms = currentRms;
        qDebug()<<QString("The %1th Iteration VF error(rms): %2").arg(10 - iterations + 1).arg(currentRms);

        for(int i=0;i<init_poles.size();i++)
        {
            init_poles[i] /= freqAverage;
            qDebug()<<QString("Pole real: %1   Pole imag: %2").arg(init_poles[i].real()).arg(init_poles[i].imag());
        }


        if(currentRms<rms)
        {
            break;
        }
        else
        {
            iterations--;
        }

    }

    return finalRms;
}

void getSparameter(std::vector<double> frequency,
                   std::vector<std::complex<double> > poles,
                   std::vector<std::vector<std::complex<double>>> residuse,
                   std::vector<double> constant,
                   std::vector<double> proportional,
                   std::vector<std::vector<std::vector<std::complex<double>>>> &sparameterData)
{

    int portNum = std::sqrt(residuse.size());
    sparameterData = std::vector<std::vector<std::vector<std::complex<double>>>>(frequency.size(),
                                                                                 std::vector<std::vector<std::complex<double>>>(portNum,std::vector<std::complex<double>>(portNum,0)));

    for(int freq = 0;freq<frequency.size();freq++)
    {
        for(int i=0;i<portNum;i++)
        {
            for(int j=0;j<portNum;j++)
            {
                std::complex<double> res = {0,0};
                for(int k=0;k<poles.size();k++)
                {
                    if(std::abs(poles[k].imag())<equalZero)
                    {
                        res += (residuse[i*portNum+j][k])/(std::complex<double>(0,2*M_PI*frequency[freq]) - poles[k]);
                    }
                    else
                    {
                        res += (residuse[i*portNum+j][k])/(std::complex<double>(0,2*M_PI*frequency[freq]) - poles[k]);
                        res += std::conj(residuse[i*portNum+j][k])/(std::complex<double>(0,2*M_PI*frequency[freq]) - std::conj(poles[k]));
                    }
                }

                if(!constant.empty())
                {
                    res += constant[i*portNum+j];
                }

                if(!proportional.empty())
                {
                    res += proportional[i*portNum+j] * std::complex<double>(0,2*M_PI*frequency[freq]);
                }

                sparameterData[freq][i][j] = res;
            }
        }
    }
}

double getRMS(std::vector<std::vector<std::vector<std::complex<double> > > > &originData, std::vector<std::vector<std::vector<std::complex<double> > > > &VFData)
{
    double sum = 0;
    for(int i=0;i<originData.size();i++)
    {
        for(int j=0;j<originData[0].size();j++)
        {
            for(int k=0;k<originData[0][0].size();k++)
            {
                sum += std::pow(originData[i][j][k].real() - VFData[i][j][k].real(),2);
                sum += std::pow(originData[i][j][k].imag() - VFData[i][j][k].imag(),2);
            }
        }
    }

    double rms = std::sqrt(sum/originData.size()/originData[0].size()/originData[0][0].size()/2);

    return rms;
}

Eigen::VectorXd getNonZeroSolutionQR(const Eigen::MatrixXd& A, double threshold) {
    int m = A.rows();
    int n = A.cols();

    // 进行列主元QR分解
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A.transpose());  // 注意：对A^T分解
    qr.setThreshold(threshold);

    int rank = qr.rank();
    int nullity = n - rank;  // 零空间维度

    if (nullity == 0) {
        // std::cout << "A满秩，只有零解" << std::endl;
        return Eigen::VectorXd::Zero(n);
    }

    // 获取QR分解的矩阵
    Eigen::MatrixXd Q = qr.householderQ();

    // 零空间基向量是Q的后nullity列
    Eigen::VectorXd solution = Q.col(n - 1);  // 取最后一个列向量

    // 归一化，避免数值过小
    if (solution.norm() > 0) {
        solution.normalize();
    }

    return solution;
}



void vectorfittingInvoke(int realPolesNum,int complexPolesNum,bool isProportional,bool isConstant,
                         int portNum,std::vector<double> frequency,
                         std::vector<std::vector<std::vector<std::complex<double> > > > sparameterData,
                         int iterations, double rms,bool isRelaxed,
                         std::vector<std::complex<double> > &poles,
                         std::vector<std::vector<std::complex<double> > > &residuse,
                         std::vector<double> &constant, std::vector<double> &proportional)
{
    double currentRms = VF(realPolesNum,complexPolesNum,isProportional,isConstant,portNum,frequency,
                           sparameterData,iterations, rms,isRelaxed,poles,residuse,constant,proportional);
    while(currentRms > rms && realPolesNum <=100 && complexPolesNum <=100)
    {
        realPolesNum++;
        complexPolesNum++;
        poles.clear();
        residuse.clear();
        constant.clear();
        proportional.clear();

        currentRms = VF(realPolesNum,complexPolesNum,isProportional,isConstant,portNum,frequency,
           sparameterData,iterations, rms,isRelaxed,poles,residuse,constant,proportional);
    }
}

bool isPassiveSparameter(std::vector<std::vector<std::vector<std::complex<double> > > >& sparameterDataRI)
{
    if(sparameterDataRI.empty())
    {
        qDebug() << QString("S parameter is empty.");
        return false;
    }
    int portNum = sparameterDataRI[0].size();
    int frequenctSize = sparameterDataRI.size();
    std::vector<double> allFrequencySingular(frequenctSize,0);
    for(int freq=0;freq<frequenctSize;freq++)
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
        allFrequencySingular[freq] = maxSingular;
        qDebug() << maxSingular;
    }
    auto maxIt = std::max_element(allFrequencySingular.begin(),allFrequencySingular.end());
    double maxSingular = *maxIt;

    if(maxSingular <= 1.0001)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void passivityEnforce(std::vector<std::vector<std::vector<std::complex<double> > > > &sparameterDataRI)
{
    if(sparameterDataRI.empty())
    {
        qDebug() << QString("S parameter is empty.");
        return;
    }
    int portNum = sparameterDataRI[0].size();
    int frequenctSize = sparameterDataRI.size();
    std::vector<double> allFrequencySingular(frequenctSize,0);
    for(int freq=0;freq<frequenctSize;freq++)
    {
        Eigen::MatrixXcd SMatrix(portNum,portNum);
        for(int i=0;i<portNum;i++)
        {
            for(int j=0;j<portNum;j++)
            {
                SMatrix(i,j) = sparameterDataRI[freq][i][j];
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(SMatrix,Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::MatrixXd singularValues = svd.singularValues().asDiagonal();
        double maxSingular = singularValues.maxCoeff();
        if(maxSingular >= 1.0001)
        {
            for(int i=0;i<portNum;i++)
            {
                for(int j=0;j<portNum;j++)
                {
                    if(singularValues(i,j)>= 1.0001)
                    {
                        singularValues(i,j) = 1;
                    }
                }
            }

            Eigen::MatrixXcd passivitySMatrix = svd.matrixU() * singularValues * svd.matrixV().transpose();

            for(int i=0;i<portNum;i++)
            {
                for(int j=0;j<portNum;j++)
                {
                    sparameterDataRI[freq][i][j] = passivitySMatrix(i,j);
                }
            }
        }
        allFrequencySingular[freq] = maxSingular;
        qDebug() << maxSingular;
    }
}

double VFOrigin(int realPolesNum, int complexPolesNum, bool isProportional, bool isConstant,
                int portNum,
                std::vector<double> frequency,
                std::vector<std::vector<std::vector<std::complex<double> > > > sparameterData,
                int iterations, double rms,bool isFrequencyNormalized,bool isWeight,
                std::vector<std::complex<double> > &poles1,
                std::vector<std::vector<std::complex<double> > > &residuse,
                std::vector<double> &constant,
                std::vector<double> &proportional)
{
    if(frequency.empty() ||
        sparameterData.empty())
    {
        qDebug() << QString("Input VF data is empty!");
        StLogger::instance()->log(QString("Input VF data is empty!"));
        return 0;
    }

    int frequencySize = frequency.size();
    std::vector<std::vector<std::complex<double>>> freqResponse(frequencySize,
                                                                std::vector<std::complex<double>>(portNum*portNum,0));

    for(int freq = 0;freq<frequencySize;freq++)
    {
        for(int row = 0;row<portNum;row++)
        {
            for(int col=0;col<portNum;col++)
            {
                freqResponse[freq][row*portNum+col] = sparameterData[freq][row][col];
            }
        }
    }

    std::vector<double> weightResponse(frequencySize,0);
    if(isWeight)
    {
        for(int freq = 0;freq<frequencySize;freq++)
        {
            double allPortSparameterSum = 0;
            for(int row = 0;row<freqResponse[0].size();row++)
            {
                allPortSparameterSum += std::pow(std::abs(freqResponse[freq][row]),2);
            }

            weightResponse[freq] = std::sqrt(allPortSparameterSum);
        }
    }


    double frequencyAverage = 0;
    for(int freq = 0;freq<frequencySize;freq++)
    {
        frequencyAverage += frequency[freq];
    }
    frequencyAverage /= frequencySize;

    double fMin = frequency[0];
    double fMax = frequency[frequencySize-1];


    double interval = (fMax - fMin)/(realPolesNum-1);
    std::vector<std::complex<double>> init_poles(realPolesNum + complexPolesNum,0);
    for(int i=0;i<realPolesNum;i++)
    {
        init_poles[i] = std::complex<double>(-2*M_PI*(fMin+i*interval),0);
    }


    if(complexPolesNum == 1)
    {
        init_poles[realPolesNum] = std::complex<double>(-0.01*2*M_PI*(fMin+fMax)/2,2*M_PI*(fMin+fMax)/2);
    }
    else
    {
        interval = (fMax - fMin)/(complexPolesNum-1);
        for(int i=0;i<complexPolesNum;i++)
        {
            init_poles[i + realPolesNum] = std::complex<double>(-0.01*2*M_PI*(fMin+i*interval),2*M_PI*(fMin+i*interval));
        }
    }

    if(isFrequencyNormalized)
    {
        for(int i=0;i<init_poles.size();i++)
        {
            init_poles[i] = init_poles[i] / frequencyAverage;
        }
    }


    for(int i=0;i<init_poles.size();i++)
    {
        init_poles[i] = init_poles[i];
        qDebug()<<QString("Pole real: %1").arg(init_poles[i].real());
        qDebug()<<QString("Pole imag: %1").arg(init_poles[i].imag());
    }


    std::vector<std::complex<double>> sFreq(frequencySize,0);
    for(int i=0;i<frequencySize;i++)
    {
        sFreq[i] = std::complex<double>(0,2*M_PI*frequency[i]);
    }

    if(isFrequencyNormalized)
    {
        for(int i=0;i<frequencySize;i++)
        {
            sFreq[i]  = sFreq[i]  / frequencyAverage;
        }
    }


    double finalRms = 0;
    int iterationNum = iterations;
    while(iterations>0)
    {

        StLogger::instance()->log(QString("The number of iterations: %1").arg(iterationNum - iterations +1));

        int colOrigin = 0;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                colOrigin += 1;
            }
            else
            {
                colOrigin += 2;
            }
        }

        int colSigma = colOrigin;

        int indexConstant = colOrigin;
        if(isConstant)
        {
            colOrigin++;
        }

        int indexProportional = colOrigin;
        if(isProportional)
        {
            colOrigin++;
        }

        std::vector<int> indexRealPolse;
        std::vector<int> indexComplexPoles;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                indexRealPolse.push_back(i);
            }
            else
            {
                indexComplexPoles.push_back(i);
            }
        }

        int realPolesNum = indexRealPolse.size();
        int complexPolesNum = indexComplexPoles.size();

        std::vector<int> idx_A_real(realPolesNum);
        for(int i=0;i<realPolesNum;i++)
        {
            idx_A_real[i] = i;
        }
        std::vector<int> idx_A_complex_re(complexPolesNum);
        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_re[i] = realPolesNum + 2*i;
        }
        std::vector<int> idx_A_complex_im(complexPolesNum);
        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_im[i] = idx_A_complex_re[i] + 1;
        }

        // 假设 s, poles 是 Eigen 列向量，index_real_poles 是 std::vector<int> 或 Eigen::VectorXi
        Eigen::VectorXcd s(sFreq.size()); // 可能是复数，例如在频域分析中
        for(int i=0;i<frequencySize;i++)
        {
            s[i]= sFreq[i];
        }

        Eigen::VectorXcd poles(init_poles.size());
        for(int i=0;i<init_poles.size();i++)
        {
            poles[i]= init_poles[i];
        }

        Eigen::VectorXi index_real_poles(indexRealPolse.size()); // 或 std::vector<int>
        for(int i=0;i<realPolesNum;i++)
        {
            index_real_poles[i]= indexRealPolse[i];
        }

        // 1. 从 poles 中选出指定的极点，并转换为行向量
        Eigen::RowVectorXcd selected_poles = poles(index_real_poles).transpose();
        // poles(index_real_poles) 是一个列向量，.transpose() 将其转为行向量

        // 2. 获取维度
        int m = s.size();
        int k = selected_poles.size(); // 也等于 index_real_poles.size()

        // 3. 通过复制扩展创建两个 (m, k) 的矩阵
        // s 作为列向量，将其复制扩展为 m x k 矩阵，每一列都是 s
        Eigen::MatrixXcd s_expanded = s.replicate(1, k);
        // selected_poles 作为行向量，将其复制扩展为 m x k 矩阵，每一行都是 selected_poles
        Eigen::MatrixXcd p_expanded = selected_poles.replicate(m, 1);

        // 4. 执行逐元素减法，得到差值矩阵
        Eigen::MatrixXcd diff_matrix = s_expanded - p_expanded;
        // 此时 diff_matrix(i, j) = s(i) - selected_poles(j)

        // 5. 对差值矩阵的每个元素求倒数（使用 .array() 进行逐元素操作）
        Eigen::MatrixXcd coeff_real = diff_matrix.array().inverse();





        Eigen::VectorXi index_complex_poles(complexPolesNum); // 复数极点的索引，k x 1
        for(int i=0;i<complexPolesNum;i++)
        {
            index_complex_poles[i]= indexComplexPoles[i];
        }

        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cR = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj = selected_poles_cR.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cR.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cR = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cR = selected_poles_cR.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded = selected_poles_conj.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数
        Eigen::MatrixXcd term1_matrix = (s_expanded_cR - p_expanded_cR).array().inverse();
        Eigen::MatrixXcd term2_matrix = (s_expanded_cR - p_conj_expanded).array().inverse();

        // 7. 两个项相加得到最终结果
        Eigen::MatrixXcd coeff_complex_re = term1_matrix + term2_matrix;






        // 定义虚数单位 1j（Python中为1j，C++中为complex(0,1)）
        std::complex<double> I(0.0, 1.0);

        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cI = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj_I = selected_poles_cI.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cI.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cI = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cI = selected_poles_cI.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded_cI = selected_poles_conj_I.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数，然后乘以虚数单位I
        Eigen::MatrixXcd term1_matrix_cI = I * (s_expanded_cI - p_expanded_cI).array().inverse();
        Eigen::MatrixXcd term2_matrix_cI = I * (s_expanded_cI - p_conj_expanded_cI).array().inverse();

        // 7. 两个项相减得到最终结果（term1 - term2）
        Eigen::MatrixXcd coeff_complex_im = term1_matrix_cI - term2_matrix_cI;




        Eigen::Tensor<std::complex<double>, 3> A(portNum*portNum, frequencySize, colOrigin+colSigma);
        A.setZero();


        // coeff_real 是二维矩阵（需要转换为Tensor）
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> coeff_real_tensor(
            coeff_real.data(), coeff_real.rows(), coeff_real.cols());

        for (int index = 0;index<idx_A_real.size();index++) {
            int idx = idx_A_real[index];
            if (idx >= 0 && idx < A.dimension(2)) {
                // 对每个第三维索引进行赋值
                for (int i = 0; i < A.dimension(0); i++) {
                    for (int j = 0; j < A.dimension(1); j++) {
                        A(i, j, idx) = coeff_real(j, index);
                    }
                }
            }
        }


        // coeff_real 是二维矩阵（需要转换为Tensor）
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> coeff_complexReal_tensor(
            coeff_complex_re.data(), coeff_complex_re.rows(), coeff_complex_re.cols());

        for (int index = 0;index<idx_A_complex_re.size();index++) {
            int idx = idx_A_complex_re[index];
            if (idx >= 0 && idx < A.dimension(2)) {
                // 对每个第三维索引进行赋值
                for (int i = 0; i < A.dimension(0); i++) {
                    for (int j = 0; j < A.dimension(1); j++) {
                        A(i, j, idx) = coeff_complexReal_tensor(j, index);
                    }
                }
            }
        }



        // coeff_real 是二维矩阵（需要转换为Tensor）
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> coeff_complexImag_tensor(
            coeff_complex_im.data(), coeff_complex_im.rows(), coeff_complex_im.cols());


        for (int index = 0;index<idx_A_complex_im.size();index++) {
            int idx = idx_A_complex_im[index];
            if (idx >= 0 && idx < A.dimension(2)) {
                // 对每个第三维索引进行赋值
                for (int i = 0; i < A.dimension(0); i++) {
                    for (int j = 0; j < A.dimension(1); j++) {
                        A(i, j, idx) = coeff_complexImag_tensor(j, index);
                    }
                }
            }
        }


        if(isConstant)
        {
            for(int i=0;i<A.dimension(0);i++)
            {
                for(int j=0;j<A.dimension(1);j++)
                {
                    A(i,j,indexConstant) = std::complex<double>(1,0);
                }
            }
        }

        if(isProportional)
        {
            for(int i=0;i<A.dimension(0);i++)
            {
                for(int j=0;j<A.dimension(1);j++)
                {
                    A(i,j,indexProportional) = sFreq[j];
                }
            }
        }



        int bRows = 2*portNum*portNum*frequencySize;

        Eigen::VectorXd b(bRows);
        b.setZero();

        for(int i=0;i<A.dimension(0);i++)
        {
            for(int j=0;j<A.dimension(1);j++)
            {
                for(int k=0;k<(int)idx_A_real.size();k++)
                {
                    A(i,j,idx_A_real[k]+colOrigin) = std::complex<double>(-1,0)*freqResponse[j][i]*coeff_real(j,k);
                }

                for(int k=0;k<(int)idx_A_complex_re.size();k++)
                {
                    A(i,j,idx_A_complex_re[k]+colOrigin) = std::complex<double>(-1,0)*freqResponse[j][i]*coeff_complex_re(j,k);
                }

                for(int k=0;k<idx_A_complex_im.size();k++)
                {
                    A(i,j,idx_A_complex_im[k]+colOrigin) = std::complex<double>(-1,0)*freqResponse[j][i]*coeff_complex_im(j,k);
                }


                int indexb = i*frequencySize+j;
                b(indexb) = freqResponse[j][i].real();
                b(portNum*portNum*frequencySize + indexb) = freqResponse[j][i].imag();

                if(isWeight)
                {
                    for(int k=0;k<A.dimension(2);k++)
                    {
                        A(i,j,k) *= weightResponse[j];

                        b(indexb) *= weightResponse[j];
                        b(portNum*portNum*frequencySize + indexb) *= weightResponse[j];
                    }
                }
            }
        }

        // A.reshape(portNum*portNum*frequencySize,colOrigin+colSigma);

        Eigen::Map<Eigen::MatrixXcd> AReCombine(A.data(), portNum*portNum*frequencySize,colOrigin+colSigma);

        // Eigen::MatrixXcd AReCombine;

        // 1. 获取矩阵维度
        int rows = AReCombine.rows();
        int cols = AReCombine.cols();

        int AReCombineRealImagRows = 2 * rows;

        // 2. 创建结果矩阵（行数为原矩阵2倍，列数相同）
        Eigen::MatrixXd AReCombineRealImag(AReCombineRealImagRows, cols);

        // 3. 将实部放在上半部分
        AReCombineRealImag.topRows(rows) = AReCombine.real();

        // 4. 将虚部放在下半部分
        AReCombineRealImag.bottomRows(rows) = AReCombine.imag();



        Eigen::VectorXd x;
        x = AReCombineRealImag.colPivHouseholderQr().solve(b);



        std::vector<double> c_res;
        for(int i=colOrigin;i<x.size();i++)
        {
            c_res.push_back(x[i]);
        }




        Eigen::MatrixXcd H(c_res.size(),c_res.size());
        H.setZero();

        std::vector<std::complex<double>> realPoles;
        std::vector<std::complex<double>> complexPoles;

        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                realPoles.push_back(init_poles[i]);
            }
            else
            {
                complexPoles.push_back(init_poles[i]);
            }
        }

        for(int i=0;i<idx_A_real.size();i++)
        {
            H(idx_A_real[i],idx_A_real[i]) = realPoles[i].real();
        }

        for(int i=0;i<idx_A_real.size();i++)
        {
            for(int j=0;j<H.cols();j++)
            {
                H(idx_A_real[i],j) -= c_res[j];
            }
        }

        for(int i=0;i<(int)idx_A_complex_re.size();i++)
        {
            for(int j=0;j<idx_A_complex_re.size();j++)
            {
                H(idx_A_complex_re[i],idx_A_complex_re[j]) = complexPoles[j].real();
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                H(idx_A_complex_re[i],idx_A_complex_im[j]) = complexPoles[j].imag();
            }
        }

        for(int i=0;i<idx_A_complex_im.size();i++)
        {
            for(int j=0;j<idx_A_complex_re.size();j++)
            {
                H(idx_A_complex_im[i],idx_A_complex_re[j]) = -1*complexPoles[j].imag();
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                H(idx_A_complex_im[i],idx_A_complex_im[j]) = complexPoles[j].real();
            }
        }

        for(int i=0;i<idx_A_complex_re.size();i++)
        {
            for(int j=0;j<H.cols();j++)
            {
                // H(idx_A_complex_re[i],j) -= 2*c_res[j]/d_res;
                H(idx_A_complex_re[i],j) -= c_res[j];
            }
        }

        Eigen::VectorXcd polesNew = H.eigenvalues();

        for(int i=0;i<polesNew.size();i++)
        {
            if(polesNew(i).real()>0)
            {
                polesNew(i) = -polesNew(i);
            }
        }

        init_poles.clear();
        for(int i=0;i<polesNew.size();i++)
        {
            if(std::abs(polesNew[i].real())<equalZero &&
                std::abs(polesNew[i].imag())<equalZero )
            {

            }
            else
            {
                init_poles.push_back(polesNew(i));
                // qDebug()<< polesNew(i).real() << polesNew(i).imag();
            }
        }



        int n_cols = 0;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                n_cols +=1;
            }
            else
            {
                n_cols+=2;
            }
        }



        std::vector<int> realIndex;
        std::vector<int> complexIndex;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                realIndex.push_back(i);
            }
            else
            {
                complexIndex.push_back(i);
            }
        }

        index_real_poles.resize(realIndex.size());
        index_complex_poles.resize(complexIndex.size());
        for(int i=0;i<index_real_poles.size();i++)
        {
            index_real_poles(i) = realIndex[i];
        }
        for(int i=0;i<index_complex_poles.size();i++)
        {
            index_complex_poles(i) = complexIndex[i];
        }

        if(isConstant)
        {
            indexConstant = n_cols;
            n_cols++;
        }

        if(isProportional)
        {
            indexProportional = n_cols;
            n_cols++;
        }


        int index = 0;
        idx_A_real.clear();
        idx_A_complex_re.clear();
        idx_A_complex_im.clear();
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                idx_A_real.push_back(index);
                index++;
            }
            else
            {
                idx_A_complex_re.push_back(index);
                idx_A_complex_im.push_back(index+1);
                index+=2;
            }
        }

        Eigen::MatrixXcd AResidues(frequencySize,n_cols);


        poles.resize(init_poles.size());
        for(int i=0;i<init_poles.size();i++)
        {
            poles[i]= init_poles[i];
        }

        // 1. 从 poles 中选出指定的极点，并转换为行向量
        Eigen::RowVectorXcd selected_poles_res = poles(index_real_poles).transpose();
        // poles(index_real_poles) 是一个列向量，.transpose() 将其转为行向量

        // 2. 获取维度
        m = s.size();
        k = selected_poles_res.size(); // 也等于 index_real_poles.size()

        // 3. 通过复制扩展创建两个 (m, k) 的矩阵
        // s 作为列向量，将其复制扩展为 m x k 矩阵，每一列都是 s
        Eigen::MatrixXcd s_expanded_res = s.replicate(1, k);
        // selected_poles 作为行向量，将其复制扩展为 m x k 矩阵，每一行都是 selected_poles
        Eigen::MatrixXcd p_expanded_res = selected_poles_res.replicate(m, 1);

        // 4. 执行逐元素减法，得到差值矩阵
        Eigen::MatrixXcd diff_matrix_res = s_expanded_res - p_expanded_res;
        // 此时 diff_matrix(i, j) = s(i) - selected_poles_res(j)

        // 5. 对差值矩阵的每个元素求倒数（使用 .array() 进行逐元素操作）
        Eigen::MatrixXcd coeff_real_res = diff_matrix_res.array().inverse();





        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cR_res = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj_res = selected_poles_cR_res.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cR_res.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cR_res = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cR_res = selected_poles_cR_res.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded_res = selected_poles_conj_res.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数
        Eigen::MatrixXcd term1_matrix_res = (s_expanded_cR_res - p_expanded_cR_res).array().inverse();
        Eigen::MatrixXcd term2_matrix_res = (s_expanded_cR_res - p_conj_expanded_res).array().inverse();

        // 7. 两个项相加得到最终结果
        Eigen::MatrixXcd coeff_complex_re_res = term1_matrix_res + term2_matrix_res;







        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cI_res = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj_I_res = selected_poles_cI_res.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cI_res.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cI_res = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cI_res = selected_poles_cI_res.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded_cI_res = selected_poles_conj_I_res.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数，然后乘以虚数单位I
        Eigen::MatrixXcd term1_matrix_cI_res = I * (s_expanded_cI_res - p_expanded_cI_res).array().inverse();
        Eigen::MatrixXcd term2_matrix_cI_res = I * (s_expanded_cI_res - p_conj_expanded_cI_res).array().inverse();

        // 7. 两个项相减得到最终结果（term1 - term2）
        Eigen::MatrixXcd coeff_complex_im_res = term1_matrix_cI_res - term2_matrix_cI_res;



        for(int i=0;i<AResidues.rows();i++)
        {
            for(int j=0;j<idx_A_real.size();j++)
            {
                AResidues(i,idx_A_real[j]) = coeff_real_res(i,j);
            }

            for(int j=0;j<idx_A_complex_re.size();j++)
            {
                AResidues(i,idx_A_complex_re[j]) = coeff_complex_re_res(i,j);
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                AResidues(i,idx_A_complex_im[j]) = coeff_complex_im_res(i,j);
            }

            if(isConstant)
            {
                AResidues(i,indexConstant) = std::complex<double>(1,0);
            }

            if(isProportional)
            {
                AResidues(i,indexProportional) = sFreq[i];
            }
        }



        // 1. 获取矩阵维度
        rows = AResidues.rows();
        cols = AResidues.cols();

        // 2. 创建结果矩阵（行数为原矩阵2倍，列数相同）
        Eigen::MatrixXd AResiduesCombineRealImag(2 * rows, cols);

        // 3. 将实部放在上半部分
        AResiduesCombineRealImag.topRows(rows) = AResidues.real();

        // 4. 将虚部放在下半部分
        AResiduesCombineRealImag.bottomRows(rows) = AResidues.imag();


        Eigen::MatrixXd freqResponseRI(frequencySize*2,freqResponse[0].size());
        for(int i=0;i<freqResponse.size();i++)
        {
            for(int j=0;j<freqResponse[0].size();j++)
            {
                freqResponseRI(i,j) = freqResponse[i][j].real();
                freqResponseRI(i+frequencySize,j) = freqResponse[i][j].imag();
            }
        }

        Eigen::MatrixXd xResiduse = AResiduesCombineRealImag.colPivHouseholderQr().solve(freqResponseRI);


        Eigen::MatrixXcd residuseTemp(portNum*portNum,index_real_poles.size() + 2*index_complex_poles.size());

        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<(int)(idx_A_real.size());j++)
            {
                residuseTemp(i,index_real_poles(j)) = xResiduse(idx_A_real[j],i);
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                residuseTemp(i,index_complex_poles(j)) = std::complex<double>(xResiduse(idx_A_complex_re[j],i),xResiduse(idx_A_complex_im[j],i));
            }
        }

        if(isConstant)
        {
            Eigen::VectorXd constant_coeff = xResiduse.row(indexConstant);

            constant.clear();
            constant.resize(constant_coeff.size());
            for(int i=0;i<constant_coeff.size();i++)
            {
                constant[i] = constant_coeff(i);
            }
        }

        if(isProportional)
        {
            Eigen::VectorXd proportional_coeff = xResiduse.row(indexProportional);

            if(isFrequencyNormalized)
            {
                for(int i=0;i<proportional_coeff.size();i++)
                {
                    proportional_coeff(i) /= frequencyAverage;
                }
            }

            proportional.clear();
            proportional.resize(proportional_coeff.size());
            for(int i=0;i<proportional_coeff.size();i++)
            {
                proportional[i] = proportional_coeff(i);
            }
        }



        residuse.clear();
        residuse = std::vector<std::vector<std::complex<double>>>(residuseTemp.rows(),std::vector<std::complex<double>>(residuseTemp.cols(),0));
        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<residuseTemp.cols();j++)
            {
                if(isFrequencyNormalized)
                {
                    residuse[i][j] = residuseTemp(i,j)*frequencyAverage;
                }
                else
                {
                    residuse[i][j] = residuseTemp(i,j);
                }
            }
        }



        poles1.clear();
        if(isFrequencyNormalized)
        {
            for(int i=0;i<init_poles.size();i++)
            {
                poles1.push_back(init_poles[i]*frequencyAverage);
            }
        }
        else
        {
            for(int i=0;i<init_poles.size();i++)
            {
                poles1.push_back(init_poles[i]);
            }
        }




        qDebug()<<QString("The %1th Iteration.").arg(10 - iterations + 1);


        std::vector<std::vector<std::vector<std::complex<double>>>> sparameterDataVF;
        getSparameter(frequency,poles1,residuse,constant,proportional,sparameterDataVF);

        double currentRms = getRMS(sparameterData,sparameterDataVF);
        finalRms = currentRms;
        qDebug()<<QString("The %1th Iteration VF error(rms): %2").arg(10 - iterations + 1).arg(currentRms);

        for(int i=0;i<init_poles.size();i++)
        {
            qDebug()<<QString("Pole real: %1   Pole imag: %2").arg(init_poles[i].real()).arg(init_poles[i].imag());
        }


        if(currentRms<rms)
        {
            break;
        }
        else
        {
            iterations--;
        }
    }

    return finalRms;
}


double VFRelaxed(int realPolesNum, int complexPolesNum, bool isProportional, bool isConstant,
                 int portNum,
                 std::vector<double> frequency,
                 std::vector<std::vector<std::vector<std::complex<double> > > > sparameterData,
                 int iterations, double rms,bool isFrequencyNormalized,bool isWeight,
                 std::vector<std::complex<double> > &poles1,
                 std::vector<std::vector<std::complex<double> > > &residuse,
                 std::vector<double> &constant,
                 std::vector<double> &proportional)
{

    if(frequency.empty() ||
        sparameterData.empty())
    {
        StLogger::instance()->log(QString("Input VF data is empty!"));
        qDebug() << QString("Input VF data is empty!");
        return 0;
    }
    int frequencySize = frequency.size();
    std::vector<std::vector<std::complex<double>>> freqResponse(frequencySize,
                                                                std::vector<std::complex<double>>(portNum*portNum,0));


    for(int freq = 0;freq<frequencySize;freq++)
    {
        for(int row = 0;row<portNum;row++)
        {
            for(int col=0;col<portNum;col++)
            {
                freqResponse[freq][row*portNum+col] = sparameterData[freq][row][col];
            }
        }
    }


    std::vector<double> weightResponse(frequencySize,0);
    double weightExtra = 0;
    if(isWeight)
    {
        for(int freq = 0;freq<frequencySize;freq++)
        {
            double allPortSparameterSum = 0;
            for(int row = 0;row<freqResponse[0].size();row++)
            {
                allPortSparameterSum += std::pow(std::abs(freqResponse[freq][row]),2);
            }

            weightResponse[freq] = std::sqrt(allPortSparameterSum);
        }

        double allSparameterSum = 0;
        for(int freq = 0;freq<frequencySize;freq++)
        {
            for(int row = 0;row<freqResponse[0].size();row++)
            {
                allSparameterSum += std::pow(std::abs(freqResponse[freq][row]) * weightResponse[freq],2);
            }
        }
        weightExtra = std::sqrt(allSparameterSum)/frequencySize/portNum/portNum;
    }


    double frequencyAverage = 0;
    for(int freq = 0;freq<frequencySize;freq++)
    {
        frequencyAverage += frequency[freq];
    }
    frequencyAverage /= frequencySize;


    double fMin = frequency[0];
    double fMax = frequency[frequencySize-1];


    // equalZero = 1/(2*M_PI*fMax*100);

    double interval = (fMax - fMin)/(realPolesNum-1);
    std::vector<std::complex<double>> init_poles(realPolesNum + complexPolesNum,0);
    for(int i=0;i<realPolesNum;i++)
    {
        init_poles[i] = std::complex<double>(-2*M_PI*(fMin+i*interval),0);
    }

    if(complexPolesNum == 1)
    {
        init_poles[realPolesNum] = std::complex<double>(-0.01*2*M_PI*(fMin+fMax)/2,2*M_PI*(fMin+fMax)/2);
    }
    else
    {
        interval = (fMax - fMin)/(complexPolesNum-1);
        for(int i=0;i<complexPolesNum;i++)
        {
            init_poles[i + realPolesNum] = std::complex<double>(-0.01*2*M_PI*(fMin+i*interval),2*M_PI*(fMin+i*interval));
        }
    }


    if(isFrequencyNormalized)
    {
        for(int i=0;i<init_poles.size();i++)
        {
            init_poles[i] = init_poles[i] / frequencyAverage;
        }
    }


    for(int i=0;i<init_poles.size();i++)
    {
        init_poles[i] = init_poles[i];
        qDebug()<<QString("Pole real: %1").arg(init_poles[i].real());
        qDebug()<<QString("Pole imag: %1").arg(init_poles[i].imag());
    }


    std::vector<std::complex<double>> sFreq(frequencySize,0);
    for(int i=0;i<frequencySize;i++)
    {
        sFreq[i] = std::complex<double>(0,2*M_PI*frequency[i]);
    }


    if(isFrequencyNormalized)
    {
        for(int i=0;i<frequencySize;i++)
        {
            sFreq[i]  = sFreq[i]  / frequencyAverage;
        }
    }

    double finalRms = 0;
    int iterationNum = iterations;
    while(iterations>0)
    {

        StLogger::instance()->log(QString("The number of iterations: %1").arg(iterationNum - iterations +1));

        int colOrigin = 0;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                colOrigin += 1;
            }
            else
            {
                colOrigin += 2;
            }
        }

        int colSigma = colOrigin + 1;

        int indexConstant = colOrigin;
        if(isConstant)
        {
            colOrigin++;
        }

        int indexProportional = colOrigin;
        if(isProportional)
        {
            colOrigin++;
        }

        std::vector<int> indexRealPolse;
        std::vector<int> indexComplexPoles;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                indexRealPolse.push_back(i);
            }
            else
            {
                indexComplexPoles.push_back(i);
            }
        }

        int realPolesNum = indexRealPolse.size();
        int complexPolesNum = indexComplexPoles.size();

        std::vector<int> idx_A_real(realPolesNum);
        for(int i=0;i<realPolesNum;i++)
        {
            idx_A_real[i] = i;
        }
        std::vector<int> idx_A_complex_re(complexPolesNum);
        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_re[i] = realPolesNum + 2*i;
        }
        std::vector<int> idx_A_complex_im(complexPolesNum);
        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_im[i] = idx_A_complex_re[i] + 1;
        }

        // 假设 s, poles 是 Eigen 列向量，index_real_poles 是 std::vector<int> 或 Eigen::VectorXi
        Eigen::VectorXcd s(sFreq.size()); // 可能是复数，例如在频域分析中
        for(int i=0;i<frequencySize;i++)
        {
            s[i]= sFreq[i];
        }

        Eigen::VectorXcd poles(init_poles.size());
        for(int i=0;i<init_poles.size();i++)
        {
            poles[i]= init_poles[i];
        }

        Eigen::VectorXi index_real_poles(indexRealPolse.size()); // 或 std::vector<int>
        for(int i=0;i<realPolesNum;i++)
        {
            index_real_poles[i]= indexRealPolse[i];
        }

        // 1. 从 poles 中选出指定的极点，并转换为行向量
        Eigen::RowVectorXcd selected_poles = poles(index_real_poles).transpose();
        // poles(index_real_poles) 是一个列向量，.transpose() 将其转为行向量

        // 2. 获取维度
        int m = s.size();
        int k = selected_poles.size(); // 也等于 index_real_poles.size()

        // 3. 通过复制扩展创建两个 (m, k) 的矩阵
        // s 作为列向量，将其复制扩展为 m x k 矩阵，每一列都是 s
        Eigen::MatrixXcd s_expanded = s.replicate(1, k);
        // selected_poles 作为行向量，将其复制扩展为 m x k 矩阵，每一行都是 selected_poles
        Eigen::MatrixXcd p_expanded = selected_poles.replicate(m, 1);

        // 4. 执行逐元素减法，得到差值矩阵
        Eigen::MatrixXcd diff_matrix = s_expanded - p_expanded;
        // 此时 diff_matrix(i, j) = s(i) - selected_poles(j)

        // 5. 对差值矩阵的每个元素求倒数（使用 .array() 进行逐元素操作）
        Eigen::MatrixXcd coeff_real = diff_matrix.array().inverse();





        Eigen::VectorXi index_complex_poles(complexPolesNum); // 复数极点的索引，k x 1
        for(int i=0;i<complexPolesNum;i++)
        {
            index_complex_poles[i]= indexComplexPoles[i];
        }

        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cR = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj = selected_poles_cR.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cR.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cR = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cR = selected_poles_cR.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded = selected_poles_conj.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数
        Eigen::MatrixXcd term1_matrix = (s_expanded_cR - p_expanded_cR).array().inverse();
        Eigen::MatrixXcd term2_matrix = (s_expanded_cR - p_conj_expanded).array().inverse();

        // 7. 两个项相加得到最终结果
        Eigen::MatrixXcd coeff_complex_re = term1_matrix + term2_matrix;






        // 定义虚数单位 1j（Python中为1j，C++中为complex(0,1)）
        std::complex<double> I(0.0, 1.0);

        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cI = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj_I = selected_poles_cI.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cI.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cI = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cI = selected_poles_cI.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded_cI = selected_poles_conj_I.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数，然后乘以虚数单位I
        Eigen::MatrixXcd term1_matrix_cI = I * (s_expanded_cI - p_expanded_cI).array().inverse();
        Eigen::MatrixXcd term2_matrix_cI = I * (s_expanded_cI - p_conj_expanded_cI).array().inverse();

        // 7. 两个项相减得到最终结果（term1 - term2）
        Eigen::MatrixXcd coeff_complex_im = term1_matrix_cI - term2_matrix_cI;




        Eigen::Tensor<std::complex<double>, 3> A(portNum*portNum, frequencySize, colOrigin+colSigma);
        A.setZero();


        // coeff_real 是二维矩阵（需要转换为Tensor）
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> coeff_real_tensor(
            coeff_real.data(), coeff_real.rows(), coeff_real.cols());

        for (int index = 0;index<idx_A_real.size();index++) {
            int idx = idx_A_real[index];
            if (idx >= 0 && idx < A.dimension(2)) {
                // 对每个第三维索引进行赋值
                for (int i = 0; i < A.dimension(0); i++) {
                    for (int j = 0; j < A.dimension(1); j++) {
                        A(i, j, idx) = coeff_real(j, index);
                    }
                }
            }
        }


        // coeff_real 是二维矩阵（需要转换为Tensor）
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> coeff_complexReal_tensor(
            coeff_complex_re.data(), coeff_complex_re.rows(), coeff_complex_re.cols());

        for (int index = 0;index<idx_A_complex_re.size();index++) {
            int idx = idx_A_complex_re[index];
            if (idx >= 0 && idx < A.dimension(2)) {
                // 对每个第三维索引进行赋值
                for (int i = 0; i < A.dimension(0); i++) {
                    for (int j = 0; j < A.dimension(1); j++) {
                        A(i, j, idx) = coeff_complexReal_tensor(j, index);
                    }
                }
            }
        }



        // coeff_real 是二维矩阵（需要转换为Tensor）
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> coeff_complexImag_tensor(
            coeff_complex_im.data(), coeff_complex_im.rows(), coeff_complex_im.cols());


        for (int index = 0;index<idx_A_complex_im.size();index++) {
            int idx = idx_A_complex_im[index];
            if (idx >= 0 && idx < A.dimension(2)) {
                // 对每个第三维索引进行赋值
                for (int i = 0; i < A.dimension(0); i++) {
                    for (int j = 0; j < A.dimension(1); j++) {
                        A(i, j, idx) = coeff_complexImag_tensor(j, index);
                    }
                }
            }
        }


        if(isConstant)
        {
            for(int i=0;i<A.dimension(0);i++)
            {
                for(int j=0;j<A.dimension(1);j++)
                {
                    A(i,j,indexConstant) = std::complex<double>(1,0);
                }
            }
        }

        if(isProportional)
        {
            for(int i=0;i<A.dimension(0);i++)
            {
                for(int j=0;j<A.dimension(1);j++)
                {
                    A(i,j,indexProportional) = sFreq[j];
                }
            }
        }


        for(int i=0;i<A.dimension(0);i++)
        {
            for(int j=0;j<A.dimension(1);j++)
            {
                for(int k=0;k<(int)idx_A_real.size();k++)
                {
                    A(i,j,idx_A_real[k]+colOrigin) = std::complex<double>(-1,0)*freqResponse[j][i]*coeff_real(j,k);
                }

                for(int k=0;k<(int)idx_A_complex_re.size();k++)
                {
                    A(i,j,idx_A_complex_re[k]+colOrigin) = std::complex<double>(-1,0)*freqResponse[j][i]*coeff_complex_re(j,k);
                }

                for(int k=0;k<idx_A_complex_im.size();k++)
                {
                    A(i,j,idx_A_complex_im[k]+colOrigin) = std::complex<double>(-1,0)*freqResponse[j][i]*coeff_complex_im(j,k);
                }


                A(i,j,A.dimension(2)-1) = std::complex<double>(-1,0)*freqResponse[j][i];

                if(isWeight)
                {
                    for(int k=0;k<A.dimension(2);k++)
                    {
                        A(i,j,k) *= weightResponse[j];
                    }
                }
            }
        }

        // A.reshape(portNum*portNum*frequencySize,colOrigin+colSigma);

        Eigen::Map<Eigen::MatrixXcd> AReCombine(A.data(), portNum*portNum*frequencySize,colOrigin+colSigma);

        // Eigen::MatrixXcd AReCombine;

        // 1. 获取矩阵维度
        int rows = AReCombine.rows();
        int cols = AReCombine.cols();

        int AReCombineRealImagRows = 2 * rows+1;


        // 2. 创建结果矩阵（行数为原矩阵2倍，列数相同）
        Eigen::MatrixXd AReCombineRealImag(AReCombineRealImagRows, cols);

        // 3. 将实部放在上半部分
        // AReCombineRealImag.topRows(rows) = AReCombine.real();
        AReCombineRealImag.block(0,0,rows,cols) = AReCombine.real();

        // 4. 将虚部放在下半部分
        // AReCombineRealImag.bottomRows(rows) = AReCombine.imag();
        AReCombineRealImag.block(rows,0,rows,cols) = AReCombine.imag();


        Eigen::RowVectorXd row1(cols);
        row1.setConstant(frequencySize);

        AReCombineRealImag.row(2*rows) = row1;




        int bRows = 0;

        bRows = 2*portNum*portNum*frequencySize + 1;

        Eigen::VectorXd b(bRows);
        b.setZero();

        double sumFrequencyReal = 0;
        double freqResponseNorm2 = 0;
        for(int i=0;i<frequencySize;i++)
        {
            for(int j=0;j<portNum;j++)
            {
                for(int k=0;k<portNum;k++)
                {
                    sumFrequencyReal += sparameterData[i][j][k].real();

                    freqResponseNorm2 += std::pow(std::abs(sparameterData[i][j][k]),2);
                }
            }
        }
        freqResponseNorm2 /= frequencySize;
        b(bRows - 1) = sumFrequencyReal;

        freqResponseNorm2 = std::sqrt(freqResponseNorm2);


        AReCombineRealImag.row(2*rows) = AReCombineRealImag.row(2*rows) * freqResponseNorm2;
        b(bRows - 1) *= freqResponseNorm2;

        if(isWeight)
        {
            b(bRows - 1) *= weightExtra;
        }


        Eigen::VectorXd x;

        x = AReCombineRealImag.colPivHouseholderQr().solve(b);



        std::vector<double> c_res;
        for(int i=colOrigin;i<x.size();i++)
        {
            c_res.push_back(x[i]);
        }

        double d_res = 0;

        d_res = x[x.size()-1];
        if(std::abs(d_res)<1e-8)
        {
            if(std::abs(d_res) < equalZero)
            {
                d_res = 1e-8;
            }
            else
            {
                d_res = 1e-8*std::abs(d_res)/d_res;
            }
        }






        Eigen::MatrixXcd H(c_res.size(),c_res.size());
        H.setZero();

        std::vector<std::complex<double>> realPoles;
        std::vector<std::complex<double>> complexPoles;

        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                realPoles.push_back(init_poles[i]);
            }
            else
            {
                complexPoles.push_back(init_poles[i]);
            }
        }

        for(int i=0;i<idx_A_real.size();i++)
        {
            H(idx_A_real[i],idx_A_real[i]) = realPoles[i].real();
        }

        for(int i=0;i<idx_A_real.size();i++)
        {
            for(int j=0;j<H.cols();j++)
            {
                H(idx_A_real[i],j) -= c_res[j]/d_res;
            }
        }

        for(int i=0;i<(int)idx_A_complex_re.size();i++)
        {
            for(int j=0;j<idx_A_complex_re.size();j++)
            {
                H(idx_A_complex_re[i],idx_A_complex_re[j]) = complexPoles[j].real();
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                H(idx_A_complex_re[i],idx_A_complex_im[j]) = complexPoles[j].imag();
            }
        }

        for(int i=0;i<idx_A_complex_im.size();i++)
        {
            for(int j=0;j<idx_A_complex_re.size();j++)
            {
                H(idx_A_complex_im[i],idx_A_complex_re[j]) = -1*complexPoles[j].imag();
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                H(idx_A_complex_im[i],idx_A_complex_im[j]) = complexPoles[j].real();
            }
        }

        for(int i=0;i<idx_A_complex_re.size();i++)
        {
            for(int j=0;j<H.cols();j++)
            {
                // H(idx_A_complex_re[i],j) -= 2*c_res[j]/d_res;
                H(idx_A_complex_re[i],j) -= c_res[j]/d_res;
            }
        }

        Eigen::VectorXcd polesNew = H.eigenvalues();

        for(int i=0;i<polesNew.size();i++)
        {
            if(polesNew(i).real()>0)
            {
                polesNew(i) = -polesNew(i);
            }
        }

        init_poles.clear();
        for(int i=0;i<polesNew.size();i++)
        {
            if(std::abs(polesNew[i].real())<equalZero &&
                std::abs(polesNew[i].imag())<equalZero)
            {

            }
            else
            {
                init_poles.push_back(polesNew(i));
                // qDebug()<< polesNew(i).real() << polesNew(i).imag();
            }
        }




        int n_cols = 0;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                n_cols +=1;
            }
            else
            {
                n_cols+=2;
            }
        }



        std::vector<int> realIndex;
        std::vector<int> complexIndex;
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                realIndex.push_back(i);
            }
            else
            {
                complexIndex.push_back(i);
            }
        }

        index_real_poles.resize(realIndex.size());
        index_complex_poles.resize(complexIndex.size());
        for(int i=0;i<index_real_poles.size();i++)
        {
            index_real_poles(i) = realIndex[i];
        }
        for(int i=0;i<index_complex_poles.size();i++)
        {
            index_complex_poles(i) = complexIndex[i];
        }

        if(isConstant)
        {
            indexConstant = n_cols;
            n_cols++;
        }

        if(isProportional)
        {
            indexProportional = n_cols;
            n_cols++;
        }


        int index = 0;
        idx_A_real.clear();
        idx_A_complex_re.clear();
        idx_A_complex_im.clear();
        for(int i=0;i<init_poles.size();i++)
        {
            if(std::abs(init_poles[i].imag())<equalZero)
            {
                idx_A_real.push_back(index);
                index++;
            }
            else
            {
                idx_A_complex_re.push_back(index);
                idx_A_complex_im.push_back(index+1);
                index+=2;
            }
        }

        Eigen::MatrixXcd AResidues(frequencySize,n_cols);


        poles.resize(init_poles.size());
        for(int i=0;i<init_poles.size();i++)
        {
            poles[i]= init_poles[i];
        }

        // 1. 从 poles 中选出指定的极点，并转换为行向量
        Eigen::RowVectorXcd selected_poles_res = poles(index_real_poles).transpose();
        // poles(index_real_poles) 是一个列向量，.transpose() 将其转为行向量

        // 2. 获取维度
        m = s.size();
        k = selected_poles_res.size(); // 也等于 index_real_poles.size()

        // 3. 通过复制扩展创建两个 (m, k) 的矩阵
        // s 作为列向量，将其复制扩展为 m x k 矩阵，每一列都是 s
        Eigen::MatrixXcd s_expanded_res = s.replicate(1, k);
        // selected_poles 作为行向量，将其复制扩展为 m x k 矩阵，每一行都是 selected_poles
        Eigen::MatrixXcd p_expanded_res = selected_poles_res.replicate(m, 1);

        // 4. 执行逐元素减法，得到差值矩阵
        Eigen::MatrixXcd diff_matrix_res = s_expanded_res - p_expanded_res;
        // 此时 diff_matrix(i, j) = s(i) - selected_poles_res(j)

        // 5. 对差值矩阵的每个元素求倒数（使用 .array() 进行逐元素操作）
        Eigen::MatrixXcd coeff_real_res = diff_matrix_res.array().inverse();





        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cR_res = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj_res = selected_poles_cR_res.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cR_res.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cR_res = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cR_res = selected_poles_cR_res.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded_res = selected_poles_conj_res.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数
        Eigen::MatrixXcd term1_matrix_res = (s_expanded_cR_res - p_expanded_cR_res).array().inverse();
        Eigen::MatrixXcd term2_matrix_res = (s_expanded_cR_res - p_conj_expanded_res).array().inverse();

        // 7. 两个项相加得到最终结果
        Eigen::MatrixXcd coeff_complex_re_res = term1_matrix_res + term2_matrix_res;







        // 1. 从 poles 中选出指定的复数极点，并转为行向量
        Eigen::RowVectorXcd selected_poles_cI_res = poles(index_complex_poles).transpose();

        // 2. 获取其共轭行向量
        Eigen::RowVectorXcd selected_poles_conj_I_res = selected_poles_cI_res.conjugate();

        // 3. 获取维度
        m = s.size();      // s 的长度
        k = selected_poles_cI_res.size();  // 复数极点的数量

        // 4. 将 s 复制扩展为 m x k 矩阵
        Eigen::MatrixXcd s_expanded_cI_res = s.replicate(1, k);

        // 5. 分别将两个极点行向量复制扩展为 m x k 矩阵
        Eigen::MatrixXcd p_expanded_cI_res = selected_poles_cI_res.replicate(m, 1);
        Eigen::MatrixXcd p_conj_expanded_cI_res = selected_poles_conj_I_res.replicate(m, 1);

        // 6. 计算两个差值矩阵并分别求倒数，然后乘以虚数单位I
        Eigen::MatrixXcd term1_matrix_cI_res = I * (s_expanded_cI_res - p_expanded_cI_res).array().inverse();
        Eigen::MatrixXcd term2_matrix_cI_res = I * (s_expanded_cI_res - p_conj_expanded_cI_res).array().inverse();

        // 7. 两个项相减得到最终结果（term1 - term2）
        Eigen::MatrixXcd coeff_complex_im_res = term1_matrix_cI_res - term2_matrix_cI_res;



        for(int i=0;i<AResidues.rows();i++)
        {
            for(int j=0;j<idx_A_real.size();j++)
            {
                AResidues(i,idx_A_real[j]) = coeff_real_res(i,j);
            }

            for(int j=0;j<idx_A_complex_re.size();j++)
            {
                AResidues(i,idx_A_complex_re[j]) = coeff_complex_re_res(i,j);
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                AResidues(i,idx_A_complex_im[j]) = coeff_complex_im_res(i,j);
            }

            if(isConstant)
            {
                AResidues(i,indexConstant) = std::complex<double>(1,0);
            }

            if(isProportional)
            {
                AResidues(i,indexProportional) = sFreq[i];
            }
        }



        // 1. 获取矩阵维度
        rows = AResidues.rows();
        cols = AResidues.cols();

        // 2. 创建结果矩阵（行数为原矩阵2倍，列数相同）
        Eigen::MatrixXd AResiduesCombineRealImag(2 * rows, cols);

        // 3. 将实部放在上半部分
        AResiduesCombineRealImag.topRows(rows) = AResidues.real();

        // 4. 将虚部放在下半部分
        AResiduesCombineRealImag.bottomRows(rows) = AResidues.imag();


        Eigen::MatrixXd freqResponseRI(frequencySize*2,freqResponse[0].size());
        for(int i=0;i<freqResponse.size();i++)
        {
            for(int j=0;j<freqResponse[0].size();j++)
            {
                freqResponseRI(i,j) = freqResponse[i][j].real();
                freqResponseRI(i+frequencySize,j) = freqResponse[i][j].imag();
            }
        }

        Eigen::MatrixXd xResiduse = AResiduesCombineRealImag.colPivHouseholderQr().solve(freqResponseRI);


        Eigen::MatrixXcd residuseTemp(portNum*portNum,index_real_poles.size() + 2*index_complex_poles.size());

        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<(int)(idx_A_real.size());j++)
            {
                residuseTemp(i,index_real_poles(j)) = xResiduse(idx_A_real[j],i);
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                residuseTemp(i,index_complex_poles(j)) = std::complex<double>(xResiduse(idx_A_complex_re[j],i),xResiduse(idx_A_complex_im[j],i));
            }
        }

        if(isConstant)
        {
            Eigen::VectorXd constant_coeff = xResiduse.row(indexConstant);

            constant.clear();
            constant.resize(constant_coeff.size());
            for(int i=0;i<constant_coeff.size();i++)
            {
                constant[i] = constant_coeff(i);
            }
        }

        if(isProportional)
        {
            Eigen::VectorXd proportional_coeff = xResiduse.row(indexProportional);


            if(isFrequencyNormalized)
            {
                for(int i=0;i<proportional_coeff.size();i++)
                {
                    proportional_coeff(i) /= frequencyAverage;
                }
            }

            proportional.clear();
            proportional.resize(proportional_coeff.size());
            for(int i=0;i<proportional_coeff.size();i++)
            {
                proportional[i] = proportional_coeff(i);
            }
        }



        residuse.clear();
        residuse = std::vector<std::vector<std::complex<double>>>(residuseTemp.rows(),std::vector<std::complex<double>>(residuseTemp.cols(),0));
        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<residuseTemp.cols();j++)
            {
                if(isFrequencyNormalized)
                {
                    residuse[i][j] = residuseTemp(i,j)*frequencyAverage;
                }
                else
                {
                    residuse[i][j] = residuseTemp(i,j);
                }
            }
        }


        poles1.clear();
        if(isFrequencyNormalized)
        {
            for(int i=0;i<init_poles.size();i++)
            {
                poles1.push_back(init_poles[i]*frequencyAverage);
            }
        }
        else
        {
            for(int i=0;i<init_poles.size();i++)
            {
                poles1.push_back(init_poles[i]);
            }
        }



        qDebug()<<QString("The %1th Iteration.").arg(10 - iterations + 1);


        std::vector<std::vector<std::vector<std::complex<double>>>> sparameterDataVF;
        getSparameter(frequency,poles1,residuse,constant,proportional,sparameterDataVF);

        double currentRms = getRMS(sparameterData,sparameterDataVF);
        finalRms = currentRms;
        qDebug()<<QString("The %1th Iteration VF error(rms): %2").arg(10 - iterations + 1).arg(currentRms);

        for(int i=0;i<init_poles.size();i++)
        {
            qDebug()<<QString("Pole real: %1   Pole imag: %2").arg(init_poles[i].real()).arg(init_poles[i].imag());
        }


        if(currentRms<rms)
        {
            break;
        }
        else
        {
            iterations--;
        }

    }

    return finalRms;
}


double vfImprove(VectorFittingInputData& inputData)
{
    int portNum = inputData.portNum;
    std::vector<double> frequency = inputData.frequency;
    std::vector<std::vector<std::vector<std::complex<double>>>> sparameterData = inputData.sparametrRI;
    int realPolesNum = inputData.realPolesNum;
    int complexPolesNum = inputData.complexPolesNum;
    double complexPolesRealCoefficient = inputData.complexPolesRealCoefficient;
    bool isComputeComplexPoles = inputData.isComputeComplexPoles;
    double dataPrecision = inputData.dataPrecision;
    bool isProportional = inputData.isProportional;
    bool isConstant = inputData.isConstant;
    bool isFrequencyNormalized = inputData.isNormalFrequency;
    int iterations = inputData.iterationNum;
    double rms = inputData.fittingError;
    bool isWeight = inputData.isWeight;
    bool isDirectComputePoles = inputData.isDirectComputePoles;
    bool isEnforcePolesStable = inputData.isEnforcePolesStable;
    bool isRelax = inputData.isRelax;
    std::vector<std::complex<double>> poles = inputData.poles1;
    std::vector<std::vector<std::complex<double>>> residuse = inputData.residuse;
    std::vector<double> constant = inputData.constant;
    std::vector<double> proportional = inputData.proportional;


    int maxThreadNum = omp_get_num_procs();

    #ifdef QT_DEBUG
    qDebug() << "计算机的最大线程数"<< maxThreadNum;
    #endif

    int threadNum = 1;
    Eigen::setNbThreads(1);


    // omp_set_num_threads(threadNum);

    if(frequency.empty() ||
        sparameterData.empty())
    {
        StLogger::instance()->log(QString("Input VF data is empty!"));
#ifdef QT_DEBUG
        qDebug() << QString("Input VF data is empty!");
#endif
        return 0;
    }

    int frequencySize = frequency.size();
    double fMin = frequency[0];
    double fMax = frequency[frequencySize-1];
    double frequencyAverage = 0;
    for(int freq = 0;freq<frequencySize;freq++)
    {
        frequencyAverage += frequency[freq];
    }
    frequencyAverage /= frequencySize;

    std::vector<std::complex<double>> sFreq(frequencySize,0);


    #pragma omp parallel for num_threads(threadNum)
    for(int i=0;i<frequencySize;i++)
    {
        sFreq[i] = std::complex<double>(0,2*M_PI*frequency[i]);
    }


    if(isFrequencyNormalized)
    {

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<frequencySize;i++)
        {
            sFreq[i]  = sFreq[i]  / frequencyAverage;
        }
    }

    std::vector<std::vector<std::complex<double>>> freqResponse(frequencySize,
                                                                std::vector<std::complex<double>>(portNum*portNum,0));


    #pragma omp parallel for num_threads(threadNum)
    for(int freq = 0;freq<frequencySize;freq++)
    {
        for(int row = 0;row<portNum;row++)
        {
            for(int col=0;col<portNum;col++)
            {
                freqResponse[freq][row*portNum+col] = sparameterData[freq][row][col];
            }
        }
    }


    std::vector<std::vector<double>> weightResponse(frequencySize,
                                                                std::vector<double>(portNum*portNum,0));
    double weightExtra = 0;
    if(isWeight)
    {
        for(int row = 0;row<portNum*portNum;row++)
        {
            double allPortSparameterSum = 0;
            for(int freq = 0;freq<frequencySize;freq++)
            {
                allPortSparameterSum += std::pow(std::abs(freqResponse[freq][row]),2);
            }
            for(int freq = 0;freq<frequencySize;freq++)
            {

                weightResponse[freq][row] = std::sqrt(std::sqrt(allPortSparameterSum));
                weightResponse[freq][row] = weightResponse[freq][row] * std::sqrt(frequency[freq]);
            }
        }


        double allSparameterSum = 0;
        for(int freq = 0;freq<frequencySize;freq++)
        {
            for(int row = 0;row<freqResponse[0].size();row++)
            {
                allSparameterSum += std::pow(std::abs(freqResponse[freq][row] * weightResponse[freq][row]),2);
            }
        }
        weightExtra = std::sqrt(allSparameterSum)/frequencySize;
    }




    // equalZero = 1/(2*M_PI*fMax*100);

    double interval = (fMax - fMin)/(realPolesNum-1);
    poles = std::vector<std::complex<double>>(realPolesNum + complexPolesNum,0);
    for(int i=0;i<realPolesNum;i++)
    {
        poles[i] = std::complex<double>(-2*M_PI*(fMin+i*interval),0);
    }

    if(complexPolesNum == 1)
    {
        poles[realPolesNum] = std::complex<double>(-complexPolesRealCoefficient*2*M_PI*(fMin+fMax)/2,2*M_PI*(fMin+fMax)/2);
    }
    else
    {
        interval = (fMax - fMin)/(complexPolesNum-1);
        for(int i=0;i<complexPolesNum;i++)
        {
            poles[i + realPolesNum] = std::complex<double>(-complexPolesRealCoefficient*2*M_PI*(fMin+i*interval),2*M_PI*(fMin+i*interval));
        }
    }


    if(isFrequencyNormalized)
    {

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<poles.size();i++)
        {
            poles[i] = poles[i] / frequencyAverage;
        }
    }



#ifdef QT_DEBUG
    qDebug()<<QString("Pole num: %1").arg(poles.size());
#endif




    double finalRms = 0;
    int iterationNum = iterations;
    while(iterations>0)
    {

        StLogger::instance()->log(QString("The number of iterations: %1").arg(iterationNum - iterations +1));

        int colOrigin = 0;
        for(int i=0;i<poles.size();i++)
        {
            if(std::abs(poles[i].imag())<equalZero)
            {
                colOrigin += 1;
            }
            else
            {
                colOrigin += 2;
            }
        }

        int colSigma = colOrigin + 1;

        int indexConstant = colOrigin;
        if(isConstant)
        {
            colOrigin++;
        }

        int indexProportional = colOrigin;
        if(isProportional)
        {
            colOrigin++;
        }

        std::vector<int> indexRealPolse;
        std::vector<int> indexComplexPoles;
        for(int i=0;i<poles.size();i++)
        {
            if(std::abs(poles[i].imag())<equalZero)
            {
                indexRealPolse.push_back(i);
            }
            else
            {
                indexComplexPoles.push_back(i);
            }
        }

        int realPolesNum = indexRealPolse.size();
        int complexPolesNum = indexComplexPoles.size();

        std::vector<int> idx_A_real(realPolesNum);

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<realPolesNum;i++)
        {
            idx_A_real[i] = i;
        }
        std::vector<int> idx_A_complex_re(complexPolesNum);

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_re[i] = realPolesNum + 2*i;
        }
        std::vector<int> idx_A_complex_im(complexPolesNum);

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_im[i] = idx_A_complex_re[i] + 1;
        }

        // 假设 s, poles 是 Eigen 列向量，index_real_poles 是 std::vector<int> 或 Eigen::VectorXi
        std::vector<std::complex<double>> s(sFreq.size()); // 可能是复数，例如在频域分析中

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<frequencySize;i++)
        {
            s[i]= sFreq[i];
        }



        Eigen::Tensor<std::complex<double>, 3> A(portNum*portNum, frequencySize, colOrigin+colSigma);
        A.setZero();



        // realPolesCoefficient = zeros(frequencyNum,length(realPolesIndexInPoles));
        std::vector<std::vector<std::complex<double>>> realPolesCoefficient(frequencySize,
            std::vector<std::complex<double>>(indexRealPolse.size(), std::complex<double>(0.0, 0.0)));

        // for freqIndex = 1:frequencyNum
        #pragma omp parallel for num_threads(threadNum)
        for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
            // for index = 1:length(realPolesIndexInPoles)
            for (int index = 0; index < (int)indexRealPolse.size(); ++index) {
                // realPolesCoefficient(freqIndex,index) = 1 / (s(freqIndex) - poles(realPolesIndexInPoles(index)));
                // Note: Subtracting 1 from Matlab indices for 0-based C++ access
                realPolesCoefficient[freqIndex][index] = 1.0 / (s[freqIndex] - poles[indexRealPolse[index]]);
            }
        }

        // complexPolesCoefficientReal = zeros(frequencyNum,length(complexPolesIndexInPoles));
        std::vector<std::vector<std::complex<double>>> complexPolesCoefficientReal(frequencySize,
            std::vector<std::complex<double>>(indexComplexPoles.size(), std::complex<double>(0.0, 0.0)));
        // complexPolesCoefficientImag = zeros(frequencyNum,length(complexPolesIndexInPoles));
        std::vector<std::vector<std::complex<double>>> complexPolesCoefficientImag(
            frequencySize,
            std::vector<std::complex<double>>(indexComplexPoles.size(), std::complex<double>(0.0, 0.0))
            );

        const std::complex<double> ej(0.0, 1.0); // 1j in Matlab

        // for freqIndex = 1:frequencyNum
        #pragma omp parallel for num_threads(threadNum)
        for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
            // for index = 1:length(complexPolesIndexInPoles)
            for (int index = 0; index < (int)indexComplexPoles.size(); ++index) {
                // complexPolesCoefficientReal(freqIndex,index) = 1 / (s(freqIndex) - poles(complexPolesIndexInPoles(index))) + 1 / (s(freqIndex) - conj(poles(complexPolesIndexInPoles(index)))) ;
                complexPolesCoefficientReal[freqIndex][index] = 1.0 / (s[freqIndex] - poles[indexComplexPoles[index]]) +
                                                                1.0 / (s[freqIndex] - std::conj(poles[indexComplexPoles[index]]));

                // complexPolesCoefficientImag(freqIndex,index) = 1j / (s(freqIndex) - poles(complexPolesIndexInPoles(index))) - 1j / (s(freqIndex) - conj(poles(complexPolesIndexInPoles(index)))) ;
                complexPolesCoefficientImag[freqIndex][index] = ej / (s[freqIndex] - poles[indexComplexPoles[index]]) -
                                                                ej / (s[freqIndex] - std::conj(poles[indexComplexPoles[index]]));
            }
        }



        // for portIndex = 1:numPorts*numPorts
        #pragma omp parallel for num_threads(threadNum)
        for (int portIndex = 0; portIndex < portNum * portNum; ++portIndex)
        {
            // for freqIndex = 1:frequencyNum
            for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {

                // for coeffIndex = 1:length(realPolesIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_real.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,realPolesIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * realPolesCoefficient(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,idx_A_real[coeffIndex]) = weightResponse[freqIndex][portIndex] * realPolesCoefficient[freqIndex][coeffIndex];
                }

                // for coeffIndex = 1:length(complexPolesRealIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_re.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,complexPolesRealIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * complexPolesCoefficientReal(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,idx_A_complex_re[coeffIndex]) = weightResponse[freqIndex][portIndex] * complexPolesCoefficientReal[freqIndex][coeffIndex];
                }

                // for coeffIndex = 1:length(complexPolesImagIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_im.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,complexPolesImagIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * complexPolesCoefficientImag(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,idx_A_complex_im[coeffIndex]) = weightResponse[freqIndex][portIndex] * complexPolesCoefficientImag[freqIndex][coeffIndex];
                }

                // if isProportional
                if (isProportional) {
                    // A(portIndex,freqIndex,proportionalIndexInA) = weight(portIndex, freqIndex) * s(freqIndex);
                    A(portIndex,freqIndex,indexProportional) = weightResponse[freqIndex][portIndex] * s[freqIndex];
                }

                // if isConstant
                if (isConstant) {
                    // A(portIndex,freqIndex,constantIndexInA) = weight(portIndex, freqIndex) * 1;
                    A(portIndex,freqIndex,indexConstant) = weightResponse[freqIndex][portIndex] * 1.0;
                }

                // for coeffIndex = 1:length(realPolesIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_real.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,AColNum + realPolesIndexInA(coeffIndex)) = -1*weight(portIndex, freqIndex) * sparameterAllPort(portIndex,freqIndex)*realPolesCoefficient(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,colOrigin + idx_A_real[coeffIndex]) = -1.0 * weightResponse[freqIndex][portIndex] * freqResponse[freqIndex][portIndex] * realPolesCoefficient[freqIndex][coeffIndex];
                }

                // for coeffIndex = 1:length(complexPolesRealIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_re.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,AColNum + complexPolesRealIndexInA(coeffIndex)) =  -1*weight(portIndex, freqIndex) * sparameterAllPort(portIndex,freqIndex)*complexPolesCoefficientReal(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,colOrigin + idx_A_complex_re[coeffIndex]) = -1.0 * weightResponse[freqIndex][portIndex] * freqResponse[freqIndex][portIndex] * complexPolesCoefficientReal[freqIndex][coeffIndex];
                }

                // for coeffIndex = 1:length(complexPolesImagIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_im.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,AColNum + complexPolesImagIndexInA(coeffIndex)) =  -1*weight(portIndex, freqIndex) * sparameterAllPort(portIndex,freqIndex)*complexPolesCoefficientImag(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,colOrigin + idx_A_complex_im[coeffIndex]) = -1.0 * weightResponse[freqIndex][portIndex] * freqResponse[freqIndex][portIndex] * complexPolesCoefficientImag[freqIndex][coeffIndex];
                }

                // A(portIndex,freqIndex,end) = -1 * weight(portIndex, freqIndex) * sparameterAllPort(portIndex,freqIndex);
                // In Matlab, 'end' refers to the last index of the dimension.
                // Here, the 3rd dimension size is (AColNum + ADeltaColNum).
                A(portIndex,freqIndex,colOrigin + colSigma-1) = -1.0 * weightResponse[freqIndex][portIndex] * freqResponse[freqIndex][portIndex];

            }
        }






        // A.reshape(portNum*portNum*frequencySize,colOrigin+colSigma);

        Eigen::Map<Eigen::MatrixXcd> AReCombine(A.data(), portNum*portNum*frequencySize,colOrigin+colSigma);

        // Eigen::MatrixXcd AReCombine;

        // 1. 获取矩阵维度
        int rows = AReCombine.rows();
        int cols = AReCombine.cols();

        int AReCombineRealImagRows = 2 * rows + 2;


        // 2. 创建结果矩阵（行数为原矩阵2倍，列数相同）
        Eigen::MatrixXd AReCombineRealImag = Eigen::MatrixXd::Zero(AReCombineRealImagRows, cols);

        // 3. 将实部放在上半部分
        // AReCombineRealImag.topRows(rows) = AReCombine.real();
        AReCombineRealImag.block(0,0,rows,cols) = AReCombine.real();

        // 4. 将虚部放在下半部分
        // AReCombineRealImag.bottomRows(rows) = AReCombine.imag();
        AReCombineRealImag.block(rows,0,rows,cols) = AReCombine.imag();


        Eigen::VectorXd realPolesCoeffRealSum = Eigen::VectorXd::Zero((int)idx_A_real.size());
        Eigen::VectorXd realPolesCoeffImagSum = Eigen::VectorXd::Zero((int)idx_A_real.size());
        Eigen::VectorXd compelxPolesCoeffRealRealSum = Eigen::VectorXd::Zero((int)idx_A_complex_re.size());
        Eigen::VectorXd compelxPolesCoeffRealImagSum = Eigen::VectorXd::Zero((int)idx_A_complex_re.size());
        Eigen::VectorXd compelxPolesCoeffImagRealSum = Eigen::VectorXd::Zero((int)idx_A_complex_im.size());
        Eigen::VectorXd compelxPolesCoeffImagImagSum = Eigen::VectorXd::Zero((int)idx_A_complex_im.size());




        // for coeffIndex = 1:length(realPolesIndexInA)
        #pragma omp parallel for num_threads(threadNum)
        for (int coeffIndex = 0; coeffIndex < (int)idx_A_real.size(); ++coeffIndex) {
            for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
                // A(portIndex,freqIndex,realPolesIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * realPolesCoefficient(freqIndex,coeffIndex);
                realPolesCoeffRealSum(coeffIndex) += realPolesCoefficient[freqIndex][coeffIndex].real();
                realPolesCoeffImagSum(coeffIndex) += realPolesCoefficient[freqIndex][coeffIndex].imag();
            }
        }

        // for coeffIndex = 1:length(complexPolesRealIndexInA)
        #pragma omp parallel for num_threads(threadNum)
        for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_re.size(); ++coeffIndex) {
            for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
                // A(portIndex,freqIndex,complexPolesRealIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * complexPolesCoefficientReal(freqIndex,coeffIndex);
                compelxPolesCoeffRealRealSum(coeffIndex) += complexPolesCoefficientReal[freqIndex][coeffIndex].real();
                compelxPolesCoeffRealImagSum(coeffIndex) += complexPolesCoefficientReal[freqIndex][coeffIndex].imag();
            }
        }

        // for coeffIndex = 1:length(complexPolesImagIndexInA)
        #pragma omp parallel for num_threads(threadNum)
        for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_im.size(); ++coeffIndex) {
            for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
                // A(portIndex,freqIndex,complexPolesImagIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * complexPolesCoefficientImag(freqIndex,coeffIndex);
                compelxPolesCoeffImagRealSum(coeffIndex) += complexPolesCoefficientImag[freqIndex][coeffIndex].real();
                compelxPolesCoeffImagImagSum(coeffIndex) += complexPolesCoefficientImag[freqIndex][coeffIndex].imag();
            }
        }


        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<idx_A_real.size();i++)
        {
            AReCombineRealImag(2*rows, colOrigin + idx_A_real[i]) = realPolesCoeffRealSum(i);
            AReCombineRealImag(2*rows+1, colOrigin + idx_A_real[i]) = realPolesCoeffImagSum(i);
        }

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<idx_A_complex_re.size();i++)
        {
            AReCombineRealImag(2*rows, colOrigin + idx_A_complex_re[i]) = compelxPolesCoeffRealRealSum(i);
            AReCombineRealImag(2*rows+1, colOrigin + idx_A_complex_re[i]) = compelxPolesCoeffRealImagSum(i);
        }

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<idx_A_complex_im.size();i++)
        {
            AReCombineRealImag(2*rows, colOrigin + idx_A_complex_im[i]) = compelxPolesCoeffImagRealSum(i);
            AReCombineRealImag(2*rows+1, colOrigin + idx_A_complex_im[i]) = compelxPolesCoeffImagImagSum(i);
        }





        int bRows = AReCombineRealImagRows;
        Eigen::VectorXd b(bRows);
        b.setZero();
        b(bRows - 2) = weightExtra*frequencySize*portNum*portNum;
        b(bRows-1) = 0;



        Eigen::VectorXd scale(AReCombineRealImag.cols());

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<AReCombineRealImag.cols();i++)
        {
            scale(i) = (AReCombineRealImag.col(i).norm());
            AReCombineRealImag.col(i) = AReCombineRealImag.col(i) / scale(i);
        }

        Eigen::VectorXd x;

        // x = AReCombineRealImag.colPivHouseholderQr().solve(b);
        x = AReCombineRealImag.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

        for(int i=0;i<x.size();i++)
        {
            x(i) = x(i)/scale(i);
        }




        std::vector<double> c_res;
        for(int i=colOrigin;i<x.size()-1;i++)
        {
            c_res.push_back(x[i]);
        }

        double d_res = 0;

        d_res = x(x.size()-1);
        if(std::abs(d_res)<1e-8)
        {
            if(std::abs(d_res) < equalZero)
            {
                d_res = 1e-8;
            }
            else
            {
                d_res = 1e-8*std::abs(d_res)/d_res;
            }
        }








        Eigen::MatrixXd AH = Eigen::MatrixXd::Zero(c_res.size(), c_res.size());
        Eigen::VectorXd bH= Eigen::VectorXd::Zero(c_res.size());

        #pragma omp parallel for num_threads(threadNum)
        for (std::size_t index = 0; index < idx_A_real.size(); ++index) {
            AH(idx_A_real[index],idx_A_real[index]) = std::real(poles[indexRealPolse[index]]);
            bH(idx_A_real[index]) = 1.0;
        }

        #pragma omp parallel for num_threads(threadNum)
        for (std::size_t index = 0; index < indexComplexPoles.size(); ++index) {
            AH(idx_A_complex_re[index],idx_A_complex_re[index]) = std::real(poles[indexComplexPoles[index]]);
            AH(idx_A_complex_re[index],idx_A_complex_im[index]) = std::imag(poles[indexComplexPoles[index]]);
            AH(idx_A_complex_im[index],idx_A_complex_re[index]) = -std::imag(poles[indexComplexPoles[index]]);
            AH(idx_A_complex_im[index],idx_A_complex_im[index]) = std::real(poles[indexComplexPoles[index]]);
        }


        #pragma omp parallel for num_threads(threadNum)
        for (std::size_t index1 = 0; index1 < idx_A_complex_re.size(); ++index1) {
            bH(idx_A_complex_re[index1]) = 2.0;
        }

        #pragma omp parallel for num_threads(threadNum)
        for (std::size_t index2 = 0; index2 < idx_A_complex_im.size(); ++index2) {
            bH(idx_A_complex_im[index2]) = 0.0;
        }

        Eigen::VectorXd c_resV= Eigen::VectorXd::Zero(c_res.size());

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<c_res.size();i++)
        {
            c_resV(i) = c_res[i];
        }

        Eigen::MatrixXd H = AH-bH * (c_resV/d_res).transpose();

        Eigen::VectorXcd polesNew = H.eigenvalues();

        if(isEnforcePolesStable)
        {
            for(int i=0;i<polesNew.size();i++)
            {
                if(polesNew(i).real()>0)
                {
                    polesNew(i) = std::complex<double>(-polesNew(i).real(), polesNew(i).imag());
                }
            }
        }

        #ifdef QT_DEBUG
        qDebug()<<QString("deal before");
        for(int i=0;i<polesNew.size();i++)
        {
            qDebug()<<QString("Pole real: %1   Pole imag: %2").arg(polesNew[i].real()).arg(polesNew[i].imag());
        }
        #endif


        poles.clear();
        for(int i=0;i<polesNew.size();i++)
        {
            if(std::abs(polesNew[i].real())<equalZero &&
                std::abs(polesNew[i].imag())<equalZero)
            {

            }
            else
            {
                if(std::abs(polesNew[i].imag()) > equalZero)
                {
                    poles.push_back(std::complex<double>(polesNew(i).real(), std::abs(polesNew(i).imag())));
                    i++;
                }
                else
                {
                    poles.push_back(std::complex<double>(polesNew(i).real(), polesNew(i).imag()));
                }
                // qDebug()<< polesNew(i).real() << polesNew(i).imag();
            }
            // poles.push_back(std::complex<double>(polesNew(i).real(), polesNew(i).imag()));
        }


        double tolImagEqualZero = 1e-12;  //设置容差判断是否小于0
        poles.clear();
        for(int i=0;i<polesNew.size();i++)
        {
            if(std::abs(polesNew[i].real())<=tolImagEqualZero &&
                std::abs(polesNew[i].imag())<=tolImagEqualZero)
            {

            }
            else
            {
                // if(std::abs(polesNew[i].imag()) > tolImagEqualZero)
                // {
                //     poles.push_back(std::complex<double>(polesNew(i).real(), std::abs(polesNew(i).imag())));
                //     i++;
                // }
                // else
                // {
                //     poles.push_back(std::complex<double>(polesNew(i).real(), polesNew(i).imag()));
                // }

                poles.push_back(std::complex<double>(polesNew(i).real(), polesNew(i).imag()));
            }
        }

        poles.erase(
            std::remove_if(poles.begin(), poles.end(),
                           [tolImagEqualZero](const std::complex<double>& c) {
                               return c.imag() < -tolImagEqualZero;   // 只有明显为负才删
                           }),
            poles.end()
            );

        std::vector<std::complex<double>>poleTemp;

        for(int index = 1;index<poles.size();index++)
        {
            if(abs(imag(poles[index])) <= tolImagEqualZero)
            {
                if(abs(real(poles[index])) > 2*M_PI*frequency[frequency.size()-1]/frequencyAverage*10
                    || abs(real(poles[index])) < 2*M_PI*frequency[0]/frequencyAverage*0.1)
                {
                }
                else
                {
                    poleTemp.push_back(poles[index]);
                }
            }
            else
            {
                if(abs(imag(poles[index])) > 2*M_PI*frequency[frequency.size()-1]/frequencyAverage*10
                    || abs(real(poles[index])) > 2*M_PI*frequency[frequency.size()-1]/frequencyAverage*complexPolesRealCoefficient*10
                    || abs(imag(poles[index])) < 2*M_PI*frequency[0]/frequencyAverage*0.01
                    ||abs(real(poles[index])) < 2*M_PI*frequency[0]/frequencyAverage*complexPolesRealCoefficient*0.01)
                {

                }
                else
                {
                    poleTemp.push_back(poles[index]);
                }
            }
        }
        if(poleTemp.empty())
        {
            break;
        }

        poles.clear();
        poles = poleTemp;


        #ifdef QT_DEBUG
        qDebug()<<QString("deal after");
        for(int i=0;i<poles.size();i++)
        {
            qDebug()<<QString("Pole real: %1   Pole imag: %2").arg(poles[i].real()).arg(poles[i].imag());
        }
        #endif




        int n_cols = 0;
        for(int i=0;i<poles.size();i++)
        {
            if(std::abs(poles[i].imag())<equalZero)
            {
                n_cols +=1;
            }
            else
            {
                n_cols+=2;
            }
        }



        std::vector<int> realIndex;
        std::vector<int> complexIndex;
        for(int i=0;i<poles.size();i++)
        {
            if(std::abs(poles[i].imag())<equalZero)
            {
                realIndex.push_back(i);
            }
            else
            {
                complexIndex.push_back(i);
            }
        }

        indexRealPolse.resize(realIndex.size());
        indexComplexPoles.resize(complexIndex.size());

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<indexRealPolse.size();i++)
        {
            indexRealPolse[i] = realIndex[i];
        }

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<indexComplexPoles.size();i++)
        {
            indexComplexPoles[i] = complexIndex[i];
        }

        if(isConstant)
        {
            indexConstant = n_cols;
            n_cols++;
        }

        if(isProportional)
        {
            indexProportional = n_cols;
            n_cols++;
        }


        int index = 0;
        idx_A_real.clear();
        idx_A_complex_re.clear();
        idx_A_complex_im.clear();
        for(int i=0;i<poles.size();i++)
        {
            if(std::abs(poles[i].imag())<equalZero)
            {
                idx_A_real.push_back(index);
                index++;
            }
            else
            {
                idx_A_complex_re.push_back(index);
                idx_A_complex_im.push_back(index+1);
                index+=2;
            }
        }

        Eigen::MatrixXcd AResidues = Eigen::MatrixXcd::Zero(frequencySize,n_cols);





        // realPolesCoefficient = zeros(frequencyNum,length(realPolesIndexInPoles));
        realPolesCoefficient.clear();
        realPolesCoefficient = std::vector<std::vector<std::complex<double>>>(frequencySize,
                                std::vector<std::complex<double>>(indexRealPolse.size(), std::complex<double>(0.0, 0.0)));

        #pragma omp parallel for num_threads(threadNum)
        for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
            for (int index = 0; index < (int)indexRealPolse.size(); ++index) {
                // realPolesCoefficient(freqIndex,index) = 1 / (s(freqIndex) - poles(realPolesIndexInPoles(index)));
                realPolesCoefficient[freqIndex][index] = 1.0 / (s[freqIndex] - poles[indexRealPolse[index]]);
            }
        }

        // complexPolesCoefficientReal = zeros(frequencyNum,length(complexPolesIndexInPoles));
        complexPolesCoefficientReal.clear();
        complexPolesCoefficientReal = std::vector<std::vector<std::complex<double>>>(frequencySize,
                                        std::vector<std::complex<double>>(indexComplexPoles.size(), std::complex<double>(0.0, 0.0)));
        // complexPolesCoefficientImag = zeros(frequencyNum,length(complexPolesIndexInPoles));
        complexPolesCoefficientImag = std::vector<std::vector<std::complex<double>>>(frequencySize,
                                        std::vector<std::complex<double>>(indexComplexPoles.size(), std::complex<double>(0.0, 0.0)));

        #pragma omp parallel for num_threads(threadNum)
        for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
            for (int index = 0; index < (int)indexComplexPoles.size(); ++index) {
                std::complex<double> currentPole = poles[indexComplexPoles[index]];
                std::complex<double> currentS = s[freqIndex];
                std::complex<double> j_unit(0.0, 1.0);

                // complexPolesCoefficientReal(freqIndex,index) = 1 / (s(freqIndex) - poles(complexPolesIndexInPoles(index))) + 1 / (s(freqIndex) - conj(poles(complexPolesIndexInPoles(index)))) ;
                complexPolesCoefficientReal[freqIndex][index] = 1.0 / (currentS - currentPole) + 1.0 / (currentS - std::conj(currentPole));
                // complexPolesCoefficientImag(freqIndex,index) = 1j / (s(freqIndex) - poles(complexPolesIndexInPoles(index))) - 1j / (s(freqIndex) - conj(poles(complexPolesIndexInPoles(index)))) ;
                complexPolesCoefficientImag[freqIndex][index] = j_unit / (currentS - currentPole) - j_unit / (currentS - std::conj(currentPole));
            }
        }


        #pragma omp parallel for num_threads(threadNum)
        for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
            for (int coeffIndex = 0; coeffIndex < (int)idx_A_real.size(); ++coeffIndex) {
                // A(freqIndex,realPolesIndexInA(coeffIndex)) = realPolesCoefficient(freqIndex,coeffIndex);
                AResidues(freqIndex,idx_A_real[coeffIndex]) = realPolesCoefficient[freqIndex][coeffIndex];
            }
            for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_re.size(); ++coeffIndex) {
                // A(freqIndex,complexPolesRealIndexInA(coeffIndex)) = complexPolesCoefficientReal(freqIndex,coeffIndex);
                AResidues(freqIndex,idx_A_complex_re[coeffIndex]) = complexPolesCoefficientReal[freqIndex][coeffIndex];
            }
            for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_im.size(); ++coeffIndex) {
                // A(freqIndex,complexPolesImagIndexInA(coeffIndex)) = complexPolesCoefficientImag(freqIndex,coeffIndex);
                AResidues(freqIndex,idx_A_complex_im[coeffIndex]) = complexPolesCoefficientImag[freqIndex][coeffIndex];
            }

            if (isProportional) {
                // A(freqIndex,proportionalIndexInA) = s(freqIndex);
                AResidues(freqIndex,indexProportional) = s[freqIndex];
            }
            if (isConstant) {
                // A(freqIndex,constantIndexInA) = 1;
                AResidues(freqIndex,indexConstant) = std::complex<double>(1.0, 0.0);
            }
        }









        // 1. 获取矩阵维度
        rows = AResidues.rows();
        cols = AResidues.cols();

        // 2. 创建结果矩阵（行数为原矩阵2倍，列数相同）
        Eigen::MatrixXd AResiduesCombineRealImag = Eigen::MatrixXd::Zero(2 * rows, cols);

        // 3. 将实部放在上半部分
        AResiduesCombineRealImag.topRows(rows) = AResidues.real();

        // 4. 将虚部放在下半部分
        AResiduesCombineRealImag.bottomRows(rows) = AResidues.imag();


        Eigen::MatrixXd freqResponseRI = Eigen::MatrixXd::Zero(frequencySize*2,freqResponse[0].size());

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<freqResponse.size();i++)
        {
            for(int j=0;j<freqResponse[0].size();j++)
            {
                freqResponseRI(i,j) = freqResponse[i][j].real();
                freqResponseRI(i+frequencySize,j) = freqResponse[i][j].imag();
            }
        }




        Eigen::VectorXd scale1 = Eigen::VectorXd::Zero(AResiduesCombineRealImag.cols());

        #pragma omp parallel for num_threads(threadNum)
        for(int i=0;i<AResiduesCombineRealImag.cols();i++)
        {
            scale1(i) = (AResiduesCombineRealImag.col(i).norm());
            AResiduesCombineRealImag.col(i) = AResiduesCombineRealImag.col(i) / scale1(i);
        }






        // Eigen::MatrixXd xResiduse = AResiduesCombineRealImag.colPivHouseholderQr().solve(freqResponseRI);
        Eigen::MatrixXd xResiduse = AResiduesCombineRealImag.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(freqResponseRI);



        // xResiduse = xResiduse.array() / scale1.array();
        for(int i=0;i<xResiduse.rows();i++)
        {
            for(int j=0;j<xResiduse.cols();j++)
            {
                xResiduse(i,j) = xResiduse(i,j) / scale1(i);
            }
        }


        Eigen::MatrixXcd residuseTemp = Eigen::MatrixXcd::Zero(portNum*portNum,indexRealPolse.size() + indexComplexPoles.size());

        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<(int)(idx_A_real.size());j++)
            {
                residuseTemp(i,indexRealPolse[j]) = xResiduse(idx_A_real[j],i);
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                residuseTemp(i,indexComplexPoles[j]) = std::complex<double>(xResiduse(idx_A_complex_re[j],i),xResiduse(idx_A_complex_im[j],i));
            }
        }

        if(isConstant)
        {
            Eigen::VectorXd constant_coeff = xResiduse.row(indexConstant);

            constant.clear();
            constant.resize(constant_coeff.size());
            for(int i=0;i<constant_coeff.size();i++)
            {
                constant[i] = constant_coeff(i);
            }
        }

        if(isProportional)
        {
            Eigen::VectorXd proportional_coeff = xResiduse.row(indexProportional);


            if(isFrequencyNormalized)
            {
                for(int i=0;i<proportional_coeff.size();i++)
                {
                    proportional_coeff(i) /= frequencyAverage;
                }
            }

            proportional.clear();
            proportional.resize(proportional_coeff.size());
            for(int i=0;i<proportional_coeff.size();i++)
            {
                proportional[i] = proportional_coeff(i);
            }
        }



        residuse.clear();
        residuse = std::vector<std::vector<std::complex<double>>>(residuseTemp.rows(),std::vector<std::complex<double>>(residuseTemp.cols(),0));
        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<residuseTemp.cols();j++)
            {
                if(isFrequencyNormalized)
                {
                    residuse[i][j] = residuseTemp(i,j)*frequencyAverage;
                }
                else
                {
                    residuse[i][j] = residuseTemp(i,j);
                }
            }
        }


        if(isFrequencyNormalized)
        {
            for(int i=0;i<poles.size();i++)
            {
                poles[i] = (poles[i]*frequencyAverage);
            }
        }




        #ifdef QT_DEBUG
        qDebug()<<QString("The %1th Iteration.").arg(10 - iterations + 1);
        #endif

        std::vector<std::vector<std::vector<std::complex<double>>>> sparameterDataVF;
        getSparameter(frequency,poles,residuse,constant,proportional,sparameterDataVF);

        double currentRms = getRMS(sparameterData,sparameterDataVF);
        finalRms = currentRms;

        #ifdef QT_DEBUG
        qDebug()<<QString("The %1th Iteration VF error(rms): %2").arg(10 - iterations + 1).arg(currentRms);


        for(int i=0;i<poles.size();i++)
        {
            qDebug()<<QString("Pole real: %1   Pole imag: %2").arg(poles[i].real()).arg(poles[i].imag());
        }
        #endif

        if(isProportional)
        {
            if(isFrequencyNormalized)
            {
                for(int i=0;i<proportional.size();i++)
                {
                    proportional[i] *= frequencyAverage;
                }
            }
        }

        if(isFrequencyNormalized)
        {
            for(int i=0;i<poles.size();i++)
            {
                poles[i] = poles[i]/frequencyAverage;
            }
        }



        if(currentRms<rms)
        {
            break;
        }
        else
        {
            iterations--;
        }

    }


    if(isProportional)
    {
        if(isFrequencyNormalized)
        {
            for(int i=0;i<proportional.size();i++)
            {
                proportional[i] /= frequencyAverage;
            }
        }
    }

    if(isFrequencyNormalized)
    {
        for(int i=0;i<poles.size();i++)
        {
            poles[i] = (poles[i]*frequencyAverage);
        }
    }

    inputData.poles1.clear();
    inputData.poles1 = poles;
    inputData.residuse.clear();
    inputData.residuse = residuse;
    inputData.constant.clear();
    inputData.constant = constant;
    inputData.proportional.clear();
    inputData.proportional = proportional;


    return finalRms;
}

/**
 * @brief Computes the number of complex poles and the poles themselves based on S-parameter data.
 *
 * @param frequency Vector of frequency points.
 * @param sparameter 3D vector representing S-parameters [port_i][port_j][freq_index].
 * @param magLimit Magnitude limit for pole detection (default 1e-2).
 * @param proportionalLimit Proportional limit for peak detection (default 1e-2).
 * @param polesNum Output parameter for the number of poles found.
 * @param poles Output vector for the detected complex poles.
 */
void computeComplexPolesNum(const std::vector<double>& frequency,
                            const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameter,
                            double magLimit,
                            double proportionalLimit,
                            int& polesNum,
                            std::vector<std::complex<double>>& poles) {

    poles.clear();

    if (sparameter.empty() || sparameter[0].empty()) {
        polesNum = 0;
        return;
    }

    size_t numPorts = sparameter[0].size();
    size_t frequencyNum = frequency.size();

    for (size_t i = 0; i < numPorts; ++i) {
        for (size_t j = 0; j < numPorts; ++j) {
            std::vector<double> real;
            std::vector<double> imag;
            for (size_t k = 0; k < frequencyNum; ++k) {
                real.push_back(sparameter[k][i][j].real());
                imag.push_back(sparameter[k][i][j].imag());
            }
            double realMax = *std::max_element(real.begin(),real.end());
            double realMin = *std::min_element(real.begin(),real.end());
            double imagMax = *std::max_element(imag.begin(),imag.end());
            double imagMin = *std::min_element(imag.begin(),imag.end());




            for (size_t k = 1; k < frequencyNum-1; ++k) {
                if ( (abs(sparameter[k][i][j].real() - sparameter[k-1][i][j].real()) + abs(sparameter[k][i][j].real() - sparameter[k+1][i][j].real()))/2 / abs(realMax - realMin) > proportionalLimit &&
                    (sparameter[k][i][j].real() - sparameter[k-1][i][j].real())/(frequency[k] - frequency[k-1]) * (sparameter[k][i][j].real() - sparameter[k+1][i][j].real())/(frequency[k] - frequency[k+1]) < 0 &&
                           abs(sparameter[k][i][j].real()) > magLimit)
                {
                    poles.push_back(std::complex<double>(sparameter[k][i][j].real(), 2*M_PI*frequency[k]));
                }
            }

            for (size_t k = 1; k < frequencyNum-1; ++k) {
                if ( (abs(sparameter[k][i][j].imag() - sparameter[k-1][i][j].imag()) + abs(sparameter[k][i][j].imag() - sparameter[k+1][i][j].imag()))/2 / abs(imagMax - imagMin) > proportionalLimit &&
                    (sparameter[k][i][j].imag() - sparameter[k-1][i][j].imag())/(frequency[k] - frequency[k-1]) * (sparameter[k][i][j].imag() - sparameter[k+1][i][j].imag())/(frequency[k] - frequency[k+1]) < 0 &&
                    abs(sparameter[k][i][j].imag()) > magLimit)
                {
                    poles.push_back(std::complex<double>(sparameter[k][i][j].imag(), 2*M_PI*frequency[k]));
                }
            }


        }
    }

    polesNum = static_cast<int>(poles.size());
}

// Overload to handle default arguments as in Matlab
void computeComplexPolesNum(const std::vector<double>& frequency,
                            const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameter,
                            int& polesNum,
                            std::vector<std::complex<double>>& poles) {
    computeComplexPolesNum(frequency, sparameter, 1e-2, 1e-2, polesNum, poles);
}

void computeComplexPolesNum(const std::vector<double>& frequency,
                            const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameter,
                            double magLimit,
                            int& polesNum,
                            std::vector<std::complex<double>>& poles) {
    computeComplexPolesNum(frequency, sparameter, magLimit, 1e-2, polesNum, poles);
}


double VFImprove(VectorFittingInputData &inputData)
{

    int portNum = inputData.portNum;
    std::vector<double> frequency = inputData.frequency;
    std::vector<std::vector<std::vector<std::complex<double>>>> sparameterData = inputData.sparametrRI;
    int realPolesNum = inputData.realPolesNum;
    int complexPolesNum = inputData.complexPolesNum;
    double complexPolesRealCoefficient = inputData.complexPolesRealCoefficient;
    // bool isComputeComplexPoles = inputData.isComputeComplexPoles;
    // double dataPrecision = inputData.dataPrecision;
    bool isProportional = inputData.isProportional;
    bool isConstant = inputData.isConstant;
    // bool isFrequencyNormalized = inputData.isNormalFrequency;
    int iterations = inputData.iterationNum;
    double rms = inputData.fittingError;
    // bool isWeight = inputData.isWeight;
    // bool isDirectComputePoles = inputData.isDirectComputePoles;
    bool isEnforcePolesStable = inputData.isEnforcePolesStable;
    // bool isRelax = inputData.isRelax;
    std::vector<std::complex<double>> poles = inputData.poles1;
    std::vector<std::vector<std::complex<double>>> residuse = inputData.residuse;
    std::vector<double> constant = inputData.constant;
    std::vector<double> proportional = inputData.proportional;










    if(frequency.empty() ||
        sparameterData.empty())
    {
        StLogger::instance()->log(QString("Input VF data is empty!"));
#ifdef QT_DEBUG
        qDebug() << QString("Input VF data is empty!");
#endif
        return 0;
    }


    int frequencySize = frequency.size();
    double fMin = frequency[0];
    double fMax = frequency[frequencySize-1];
    double frequencyAverage = 0;
    for(int freq = 0;freq<frequencySize;freq++)
    {
        frequencyAverage += frequency[freq];
    }
    frequencyAverage /= frequencySize;

    std::vector<std::complex<double>> sFreq(frequencySize,0);



    for(int i=0;i<frequencySize;i++)
    {
        sFreq[i] = std::complex<double>(0,2*M_PI*frequency[i]/ frequencyAverage);
    }

    double tolImagEqualZero = 1e-12;  //设置容差判断是否小于0

    std::vector<std::vector<std::complex<double>>> freqResponse(frequencySize,
                                                                std::vector<std::complex<double>>(portNum*portNum,0));


    for(int freq = 0;freq<frequencySize;freq++)
    {
        for(int row = 0;row<portNum;row++)
        {
            for(int col=0;col<portNum;col++)
            {
                freqResponse[freq][row*portNum+col] = sparameterData[freq][row][col];
            }
        }
    }


    std::vector<std::vector<double>> weightResponse(frequencySize,
                                                    std::vector<double>(portNum*portNum,0));
    double weightExtra = 0;
    // if(isWeight)
    {
        for(int row = 0;row<portNum*portNum;row++)
        {
            for(int freq = 0;freq<frequencySize;freq++)
            {

                weightResponse[freq][row] = 1/(std::abs(freqResponse[freq][row])+1e-12);
            }
        }


        double allSparameterSum = 0;
        for(int freq = 0;freq<frequencySize;freq++)
        {
            for(int row = 0;row<freqResponse[0].size();row++)
            {
                allSparameterSum += std::pow(std::abs(freqResponse[freq][row] * weightResponse[freq][row]),2);
            }
        }
        weightExtra = std::sqrt(allSparameterSum)/frequencySize;
    }




    // equalZero = 1/(2*M_PI*fMax*100);

    double interval = (fMax - fMin)/(realPolesNum-1);
    poles = std::vector<std::complex<double>>(realPolesNum + complexPolesNum,0);
    for(int i=0;i<realPolesNum;i++)
    {
        poles[i] = std::complex<double>(-2*M_PI*(fMin+i*interval),0);
    }

    if(complexPolesNum == 1)
    {
        poles[realPolesNum] = std::complex<double>(-std::abs(complexPolesRealCoefficient)*2*M_PI*(fMin+fMax)/2,2*M_PI*(fMin+fMax)/2);
    }
    else
    {
        interval = (fMax - fMin)/(complexPolesNum-1);
        for(int i=0;i<complexPolesNum;i++)
        {
            poles[i + realPolesNum] = std::complex<double>(-std::abs(complexPolesRealCoefficient)*2*M_PI*(fMin+i*interval),2*M_PI*(fMin+i*interval));
        }
    }


    // if(isFrequencyNormalized)
    {

        for(int i=0;i<poles.size();i++)
        {
            poles[i] = poles[i] / frequencyAverage;
        }
    }



#ifdef QT_DEBUG
    qDebug()<<QString("Pole num: %1").arg(poles.size());
#endif




    double finalRms = 0;
    int iterationNum = iterations;
    while(iterations>0)
    {

        // StLogger::instance()->log(QString("The number of iterations: %1").arg(iterationNum - iterations +1));

        std::vector<int> indexRealPolse;
        std::vector<int> indexComplexPoles;
        for(int i=0;i<poles.size();i++)
        {
            if(std::abs(poles[i].imag())<=tolImagEqualZero)
            {
                indexRealPolse.push_back(i);
            }
            else
            {
                indexComplexPoles.push_back(i);
            }
        }

        int colOrigin = indexRealPolse.size() + 2*indexComplexPoles.size();
        int colSigma = colOrigin + 1;

        int indexConstant = colOrigin;
        if(isConstant)
        {
            colOrigin++;
        }

        int indexProportional = colOrigin;
        if(isProportional)
        {
            colOrigin++;
        }



        int realPolesNum = indexRealPolse.size();
        int complexPolesNum = indexComplexPoles.size();

        std::vector<int> idx_A_real(realPolesNum);

        for(int i=0;i<realPolesNum;i++)
        {
            idx_A_real[i] = i;
        }
        std::vector<int> idx_A_complex_re(complexPolesNum);

        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_re[i] = realPolesNum + 2*i;
        }
        std::vector<int> idx_A_complex_im(complexPolesNum);

        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_im[i] = idx_A_complex_re[i] + 1;
        }

        // 假设 s, poles 是 Eigen 列向量，index_real_poles 是 std::vector<int> 或 Eigen::VectorXi
        std::vector<std::complex<double>> s(sFreq.size()); // 可能是复数，例如在频域分析中

        for(int i=0;i<frequencySize;i++)
        {
            s[i]= sFreq[i];
        }



        Eigen::Tensor<std::complex<double>, 3> A(portNum*portNum, frequencySize, colOrigin+colSigma);
        A.setZero();



        // realPolesCoefficient = zeros(frequencyNum,length(realPolesIndexInPoles));
        std::vector<std::vector<std::complex<double>>> realPolesCoefficient(frequencySize,
                                                                            std::vector<std::complex<double>>(indexRealPolse.size(), std::complex<double>(0.0, 0.0)));

// for freqIndex = 1:frequencyNum
        for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
            // for index = 1:length(realPolesIndexInPoles)
            for (int index = 0; index < (int)indexRealPolse.size(); ++index) {
                // realPolesCoefficient(freqIndex,index) = 1 / (s(freqIndex) - poles(realPolesIndexInPoles(index)));
                // Note: Subtracting 1 from Matlab indices for 0-based C++ access
                realPolesCoefficient[freqIndex][index] = 1.0 / (s[freqIndex] - poles[indexRealPolse[index]]);
            }
        }

        // complexPolesCoefficientReal = zeros(frequencyNum,length(complexPolesIndexInPoles));
        std::vector<std::vector<std::complex<double>>> complexPolesCoefficientReal(
            frequencySize,
            std::vector<std::complex<double>>(indexComplexPoles.size(), std::complex<double>(0.0, 0.0)));
        // complexPolesCoefficientImag = zeros(frequencyNum,length(complexPolesIndexInPoles));
        std::vector<std::vector<std::complex<double>>> complexPolesCoefficientImag(
            frequencySize,
            std::vector<std::complex<double>>(indexComplexPoles.size(), std::complex<double>(0.0, 0.0))
            );

        const std::complex<double> ej(0.0, 1.0); // 1j in Matlab

// for freqIndex = 1:frequencyNum
        for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
            // for index = 1:length(complexPolesIndexInPoles)
            for (int index = 0; index < (int)indexComplexPoles.size(); ++index) {
                // complexPolesCoefficientReal(freqIndex,index) = 1 / (s(freqIndex) - poles(complexPolesIndexInPoles(index))) + 1 / (s(freqIndex) - conj(poles(complexPolesIndexInPoles(index)))) ;
                complexPolesCoefficientReal[freqIndex][index] = 1.0 / (s[freqIndex] - poles[indexComplexPoles[index]]) +
                                                                1.0 / (s[freqIndex] - std::conj(poles[indexComplexPoles[index]]));

                // complexPolesCoefficientImag(freqIndex,index) = 1j / (s(freqIndex) - poles(complexPolesIndexInPoles(index))) - 1j / (s(freqIndex) - conj(poles(complexPolesIndexInPoles(index)))) ;
                complexPolesCoefficientImag[freqIndex][index] = ej / (s[freqIndex] - poles[indexComplexPoles[index]]) -
                                                                ej / (s[freqIndex] - std::conj(poles[indexComplexPoles[index]]));
            }
        }



// for portIndex = 1:numPorts*numPorts
        for (int portIndex = 0; portIndex < portNum * portNum; ++portIndex)
        {
            // for freqIndex = 1:frequencyNum
            for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {

                // for coeffIndex = 1:length(realPolesIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_real.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,realPolesIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * realPolesCoefficient(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,idx_A_real[coeffIndex]) = weightResponse[freqIndex][portIndex] * realPolesCoefficient[freqIndex][coeffIndex];
                }

                // for coeffIndex = 1:length(complexPolesRealIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_re.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,complexPolesRealIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * complexPolesCoefficientReal(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,idx_A_complex_re[coeffIndex]) = weightResponse[freqIndex][portIndex] * complexPolesCoefficientReal[freqIndex][coeffIndex];
                }

                // for coeffIndex = 1:length(complexPolesImagIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_im.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,complexPolesImagIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * complexPolesCoefficientImag(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,idx_A_complex_im[coeffIndex]) = weightResponse[freqIndex][portIndex] * complexPolesCoefficientImag[freqIndex][coeffIndex];
                }

                // if isProportional
                if (isProportional) {
                    // A(portIndex,freqIndex,proportionalIndexInA) = weight(portIndex, freqIndex) * s(freqIndex);
                    A(portIndex,freqIndex,indexProportional) = weightResponse[freqIndex][portIndex] * s[freqIndex];
                }

                // if isConstant
                if (isConstant) {
                    // A(portIndex,freqIndex,constantIndexInA) = weight(portIndex, freqIndex) * 1;
                    A(portIndex,freqIndex,indexConstant) = weightResponse[freqIndex][portIndex] * 1.0;
                }

                // for coeffIndex = 1:length(realPolesIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_real.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,AColNum + realPolesIndexInA(coeffIndex)) = -1*weight(portIndex, freqIndex) * sparameterAllPort(portIndex,freqIndex)*realPolesCoefficient(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,colOrigin + idx_A_real[coeffIndex]) = -1.0 * weightResponse[freqIndex][portIndex] * freqResponse[freqIndex][portIndex] * realPolesCoefficient[freqIndex][coeffIndex];
                }

                // for coeffIndex = 1:length(complexPolesRealIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_re.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,AColNum + complexPolesRealIndexInA(coeffIndex)) =  -1*weight(portIndex, freqIndex) * sparameterAllPort(portIndex,freqIndex)*complexPolesCoefficientReal(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,colOrigin + idx_A_complex_re[coeffIndex]) = -1.0 * weightResponse[freqIndex][portIndex] * freqResponse[freqIndex][portIndex] * complexPolesCoefficientReal[freqIndex][coeffIndex];
                }

                // for coeffIndex = 1:length(complexPolesImagIndexInA)
                for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_im.size(); ++coeffIndex) {
                    // A(portIndex,freqIndex,AColNum + complexPolesImagIndexInA(coeffIndex)) =  -1*weight(portIndex, freqIndex) * sparameterAllPort(portIndex,freqIndex)*complexPolesCoefficientImag(freqIndex,coeffIndex);
                    A(portIndex,freqIndex,colOrigin + idx_A_complex_im[coeffIndex]) = -1.0 * weightResponse[freqIndex][portIndex] * freqResponse[freqIndex][portIndex] * complexPolesCoefficientImag[freqIndex][coeffIndex];
                }

                // A(portIndex,freqIndex,end) = -1 * weight(portIndex, freqIndex) * sparameterAllPort(portIndex,freqIndex);
                // In Matlab, 'end' refers to the last index of the dimension.
                // Here, the 3rd dimension size is (AColNum + ADeltaColNum).
                A(portIndex,freqIndex,colOrigin + colSigma-1) = -1.0 * weightResponse[freqIndex][portIndex] * freqResponse[freqIndex][portIndex];

            }
        }



        int dim_m = 2 * frequencySize;
        int dim_n = colOrigin + colSigma;
        int dim_k = std::min(dim_m, dim_n);

        Eigen::Tensor<double, 3> R(portNum*portNum, dim_k, dim_n);
        R.setZero();

        Eigen::Tensor<double, 3> realA = A.real();  // 强制求值
        Eigen::Tensor<double, 3> imagA = A.imag();

        Eigen::Tensor<double, 3> result = realA.concatenate(imagA, 1);


        for(int i=0;i<result.dimension(0);i++)
        {
            Eigen::MatrixXd resultTemp(result.dimension(1),result.dimension(2));
            for(int j=0;j<result.dimension(1);j++)
            {
                for(int k=0;k<result.dimension(2);k++)
                {
                    resultTemp(j,k) = result(i,j,k);
                }
            }

            // 3. 进行 QR 分解
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(resultTemp);

            // 获取 Q 和 R
            Eigen::MatrixXd Q = qr.householderQ();
            Eigen::MatrixXd R1 = qr.matrixR().triangularView<Eigen::Upper>();

            for(int j=0;j<dim_k;j++)
            {
                for(int k=0;k<result.dimension(2);k++)
                {
                    R(i,j,k) = R1(j,k);
                }
            }
        }

        int n_rows_r12 = 0;
        int n_rows_r22 =0;
        if (dim_k == dim_m)
        {
            n_rows_r12 = frequencySize;
            n_rows_r22 = frequencySize;
        }
        else
        {
            n_rows_r12 = colOrigin;
            n_rows_r22 = colSigma;
        }

        Eigen::Tensor<double, 3> R22(portNum*portNum, dim_k-n_rows_r12, colSigma);

        for(int i=0;i<result.dimension(0);i++)
        {
            for(int j=0;j<R22.dimension(1);j++)
            {
                for(int k=0;k<R22.dimension(2);k++)
                {
                    R22(i,j,k) = R(i,n_rows_r12 + j,dim_n - colSigma + k);
                }
            }
        }



        int dim0 = portNum*portNum * n_rows_r22;

        Eigen::MatrixXd A_fast(dim0, colSigma);
        A_fast.setZero();
        Eigen::Map<Eigen::MatrixXd> B(R22.data(), dim0, colSigma);
        A_fast.block(0,0,dim0,colSigma) = B;





        Eigen::MatrixXd AReCombineRealImag(dim0+2,colSigma);
        AReCombineRealImag.setZero();
        AReCombineRealImag.block(0,0,dim0,colSigma) = A_fast;


        Eigen::VectorXd realPolesCoeffRealSum = Eigen::VectorXd::Zero((int)idx_A_real.size());
        Eigen::VectorXd realPolesCoeffImagSum = Eigen::VectorXd::Zero((int)idx_A_real.size());
        Eigen::VectorXd compelxPolesCoeffRealRealSum = Eigen::VectorXd::Zero((int)idx_A_complex_re.size());
        Eigen::VectorXd compelxPolesCoeffRealImagSum = Eigen::VectorXd::Zero((int)idx_A_complex_re.size());
        Eigen::VectorXd compelxPolesCoeffImagRealSum = Eigen::VectorXd::Zero((int)idx_A_complex_im.size());
        Eigen::VectorXd compelxPolesCoeffImagImagSum = Eigen::VectorXd::Zero((int)idx_A_complex_im.size());




        for (int coeffIndex = 0; coeffIndex < (int)idx_A_real.size(); ++coeffIndex) {
            for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
                // A(portIndex,freqIndex,realPolesIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * realPolesCoefficient(freqIndex,coeffIndex);
                realPolesCoeffRealSum(coeffIndex) += realPolesCoefficient[freqIndex][coeffIndex].real();
                realPolesCoeffImagSum(coeffIndex) += realPolesCoefficient[freqIndex][coeffIndex].imag();
            }
        }

        for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_re.size(); ++coeffIndex) {
            for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
                // A(portIndex,freqIndex,complexPolesRealIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * complexPolesCoefficientReal(freqIndex,coeffIndex);
                compelxPolesCoeffRealRealSum(coeffIndex) += complexPolesCoefficientReal[freqIndex][coeffIndex].real();
                compelxPolesCoeffRealImagSum(coeffIndex) += complexPolesCoefficientReal[freqIndex][coeffIndex].imag();
            }
        }

        for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_im.size(); ++coeffIndex) {
            for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
                // A(portIndex,freqIndex,complexPolesImagIndexInA(coeffIndex)) = weight(portIndex, freqIndex) * complexPolesCoefficientImag(freqIndex,coeffIndex);
                compelxPolesCoeffImagRealSum(coeffIndex) += complexPolesCoefficientImag[freqIndex][coeffIndex].real();
                compelxPolesCoeffImagImagSum(coeffIndex) += complexPolesCoefficientImag[freqIndex][coeffIndex].imag();
            }
        }


        for(int i=0;i<idx_A_real.size();i++)
        {
            AReCombineRealImag(dim0,   idx_A_real[i]) = realPolesCoeffRealSum(i);
            AReCombineRealImag(dim0+1,   idx_A_real[i]) = realPolesCoeffImagSum(i);
        }

        for(int i=0;i<idx_A_complex_re.size();i++)
        {
            AReCombineRealImag(dim0,  idx_A_complex_re[i]) = compelxPolesCoeffRealRealSum(i);
            AReCombineRealImag(dim0+1,   idx_A_complex_re[i]) = compelxPolesCoeffRealImagSum(i);
        }

        for(int i=0;i<idx_A_complex_im.size();i++)
        {
            AReCombineRealImag(dim0,   idx_A_complex_im[i]) = compelxPolesCoeffImagRealSum(i);
            AReCombineRealImag(dim0+1,   idx_A_complex_im[i]) = compelxPolesCoeffImagImagSum(i);
        }
        AReCombineRealImag(dim0,   colSigma- 1) = 1;
        AReCombineRealImag(dim0+1,   colSigma- 1) = 1;





        int bRows = dim0+2;
        Eigen::VectorXd b(bRows);
        b.setZero();
        b(bRows - 2) = weightExtra*frequencySize*weightExtra*weightExtra*weightExtra;
        b(bRows-1) = 0;

        AReCombineRealImag.row(bRows-2) = AReCombineRealImag.row(bRows-2)*weightExtra;
        AReCombineRealImag.row(bRows-1) = AReCombineRealImag.row(bRows-1)*weightExtra;



        Eigen::VectorXd scale(AReCombineRealImag.cols());

        for(int i=0;i<AReCombineRealImag.cols();i++)
        {
            scale(i) = (AReCombineRealImag.col(i).norm());
            AReCombineRealImag.col(i) = AReCombineRealImag.col(i) / scale(i);
        }

        Eigen::VectorXd x;

        x = AReCombineRealImag.colPivHouseholderQr().solve(b);
        // x = AReCombineRealImag.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

        for(int i=0;i<x.size();i++)
        {
            if(std::isnan(x(i)) || std::isinf(x(i)))
            {
                break;
            }
            x(i) = x(i)/scale(i);
        }




        std::vector<double> c_res;
        for(int i=0;i<x.size()-1;i++)
        {
            c_res.push_back(x[i]);
        }

        double d_res = 0;

        d_res = x(x.size()-1);
        if(std::abs(d_res)<1e-8)
        {
            if(std::abs(d_res) < tolImagEqualZero)
            {
                d_res = 1e-8;
            }
            else
            {
                d_res = 1e-8*std::abs(d_res)/d_res;
            }
        }








        Eigen::MatrixXd AH = Eigen::MatrixXd::Zero(c_res.size(), c_res.size());
        Eigen::VectorXd bH= Eigen::VectorXd::Zero(c_res.size());

        for (std::size_t index = 0; index < idx_A_real.size(); ++index) {
            AH(idx_A_real[index],idx_A_real[index]) = poles[indexRealPolse[index]].real();
            bH(idx_A_real[index]) = 1.0;
        }

        for (std::size_t index = 0; index < indexComplexPoles.size(); ++index) {
            AH(idx_A_complex_re[index],idx_A_complex_re[index]) = poles[indexComplexPoles[index]].real();
            AH(idx_A_complex_re[index],idx_A_complex_im[index]) = poles[indexComplexPoles[index]].imag();
            AH(idx_A_complex_im[index],idx_A_complex_re[index]) = -poles[indexComplexPoles[index]].imag();
            AH(idx_A_complex_im[index],idx_A_complex_im[index]) = poles[indexComplexPoles[index]].real();
        }


        for (std::size_t index1 = 0; index1 < idx_A_complex_re.size(); ++index1) {
            bH(idx_A_complex_re[index1]) = 2.0;
        }

        for (std::size_t index2 = 0; index2 < idx_A_complex_im.size(); ++index2) {
            bH(idx_A_complex_im[index2]) = 0.0;
        }
        if(poles.empty())
        {
            break;
        }

        Eigen::VectorXd c_resV= Eigen::VectorXd::Zero(c_res.size());

        for(int i=0;i<c_res.size();i++)
        {
            c_resV(i) = c_res[i];
        }

        Eigen::MatrixXd H = AH-bH * (c_resV/d_res).transpose();

        Eigen::VectorXcd polesNew = H.eigenvalues();

        if(isEnforcePolesStable)
        {
            for(int i=0;i<polesNew.size();i++)
            {
                if(polesNew(i).real()>0)
                {
                    polesNew(i) = std::complex<double>(-polesNew(i).real(), polesNew(i).imag());
                }
            }
        }

// #ifdef QT_DEBUG
        qDebug()<<QString("deal before");
        for(int i=0;i<polesNew.size();i++)
        {
            qDebug()<<QString("Pole real: %1   Pole imag: %2").arg(polesNew[i].real()).arg(polesNew[i].imag());
        }
// #endif

//筛选和matlab不一致
        /////
        /// /////
        /// /////
        /// /////
        poles.clear();
        for(int i=0;i<polesNew.size();i++)
        {
            if(std::abs(polesNew[i].real())<=tolImagEqualZero &&
                std::abs(polesNew[i].imag())<=tolImagEqualZero)
            {

            }
            else
            {

                poles.push_back(std::complex<double>(polesNew(i).real(), polesNew(i).imag()));
            }
        }

        poles.erase(
            std::remove_if(poles.begin(), poles.end(),
                           [tolImagEqualZero](const std::complex<double>& c) {
                               return c.imag() < -tolImagEqualZero;   // 只有明显为负才删
                           }),
            poles.end()
            );

        std::vector<std::complex<double>>poleTemp;

        for(int index = 1;index<poles.size();index++)
        {
            if(abs(imag(poles[index])) <= tolImagEqualZero)
            {
                if(abs(real(poles[index])) > 2*M_PI*frequency[frequency.size()-1]/frequencyAverage*10
                    || abs(real(poles[index])) < 2*M_PI*frequency[0]/frequencyAverage*0.1)
                {
                }
                else
                {
                    poleTemp.push_back(poles[index]);
                }
            }
            else
            {
                if(abs(imag(poles[index])) > 2*M_PI*frequency[frequency.size()-1]/frequencyAverage*10
                        || abs(real(poles[index])) > 2*M_PI*frequency[frequency.size()-1]/frequencyAverage*complexPolesRealCoefficient*10
                        || abs(imag(poles[index])) < 2*M_PI*frequency[0]/frequencyAverage*0.01
                        ||abs(real(poles[index])) < 2*M_PI*frequency[0]/frequencyAverage*complexPolesRealCoefficient*0.01)
                {

                }
                else
                {
                    poleTemp.push_back(poles[index]);
                }
            }
        }
        if(poleTemp.empty())
        {
            break;
        }

        poles.clear();
        poles = poleTemp;


// #ifdef QT_DEBUG
        qDebug()<<QString("deal after");
        for(int i=0;i<poles.size();i++)
        {
            qDebug()<<QString("Pole real: %1   Pole imag: %2").arg(poles[i].real()).arg(poles[i].imag());
        }
// #endif




        std::vector<int> realIndex;
        std::vector<int> complexIndex;
        for(int i=0;i<poles.size();i++)
        {
            if(std::abs(poles[i].imag())<=tolImagEqualZero)
            {
                realIndex.push_back(i);
            }
            else
            {
                complexIndex.push_back(i);
            }
        }

        int n_cols = realIndex.size() + 2*complexIndex.size();
        if(isConstant)
        {
            indexConstant = n_cols;
            n_cols++;
        }

        if(isProportional)
        {
            indexProportional = n_cols;
            n_cols++;
        }

        indexRealPolse.resize(realIndex.size());
        indexComplexPoles.resize(complexIndex.size());

        for(int i=0;i<indexRealPolse.size();i++)
        {
            indexRealPolse[i] = realIndex[i];
        }

        for(int i=0;i<indexComplexPoles.size();i++)
        {
            indexComplexPoles[i] = complexIndex[i];
        }




        int index = 0;
        idx_A_real.clear();
        idx_A_complex_re.clear();
        idx_A_complex_im.clear();
        realPolesNum = indexRealPolse.size();
        complexPolesNum = indexComplexPoles.size();
        idx_A_real.resize(realPolesNum,0);

        for(int i=0;i<realPolesNum;i++)
        {
            idx_A_real[i] = i;
        }
        idx_A_complex_re.resize(complexPolesNum, 0);

        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_re[i] = realPolesNum + 2*i;
        }
        idx_A_complex_im.resize(complexPolesNum, 0);

        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_im[i] = idx_A_complex_re[i] + 1;
        }

        Eigen::MatrixXcd AResidues = Eigen::MatrixXcd::Zero(frequencySize,n_cols);





        // realPolesCoefficient = zeros(frequencyNum,length(realPolesIndexInPoles));
        realPolesCoefficient.clear();
        realPolesCoefficient = std::vector<std::vector<std::complex<double>>>(frequencySize,
                                                                              std::vector<std::complex<double>>(indexRealPolse.size(), std::complex<double>(0.0, 0.0)));

        for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
            for (int index = 0; index < (int)indexRealPolse.size(); ++index) {
                // realPolesCoefficient(freqIndex,index) = 1 / (s(freqIndex) - poles(realPolesIndexInPoles(index)));
                realPolesCoefficient[freqIndex][index] = 1.0 / (s[freqIndex] - poles[indexRealPolse[index]]);
            }
        }

        // complexPolesCoefficientReal = zeros(frequencyNum,length(complexPolesIndexInPoles));
        complexPolesCoefficientReal.clear();
        complexPolesCoefficientReal = std::vector<std::vector<std::complex<double>>>(frequencySize,
                                                                                     std::vector<std::complex<double>>(indexComplexPoles.size(), std::complex<double>(0.0, 0.0)));
        // complexPolesCoefficientImag = zeros(frequencyNum,length(complexPolesIndexInPoles));
        complexPolesCoefficientImag.clear();
        complexPolesCoefficientImag = std::vector<std::vector<std::complex<double>>>(frequencySize,
                                                                                     std::vector<std::complex<double>>(indexComplexPoles.size(), std::complex<double>(0.0, 0.0)));

        for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
            for (int index = 0; index < (int)indexComplexPoles.size(); ++index) {
                std::complex<double> currentPole = poles[indexComplexPoles[index]];
                std::complex<double> currentS = s[freqIndex];
                std::complex<double> j_unit(0.0, 1.0);

                // complexPolesCoefficientReal(freqIndex,index) = 1 / (s(freqIndex) - poles(complexPolesIndexInPoles(index))) + 1 / (s(freqIndex) - conj(poles(complexPolesIndexInPoles(index)))) ;
                complexPolesCoefficientReal[freqIndex][index] = 1.0 / (currentS - currentPole) + 1.0 / (currentS - std::conj(currentPole));
                // complexPolesCoefficientImag(freqIndex,index) = 1j / (s(freqIndex) - poles(complexPolesIndexInPoles(index))) - 1j / (s(freqIndex) - conj(poles(complexPolesIndexInPoles(index)))) ;
                complexPolesCoefficientImag[freqIndex][index] = j_unit / (currentS - currentPole) - j_unit / (currentS - std::conj(currentPole));
            }
        }


        for (int freqIndex = 0; freqIndex < frequencySize; ++freqIndex) {
            for (int coeffIndex = 0; coeffIndex < (int)idx_A_real.size(); ++coeffIndex) {
                // A(freqIndex,realPolesIndexInA(coeffIndex)) = realPolesCoefficient(freqIndex,coeffIndex);
                AResidues(freqIndex,idx_A_real[coeffIndex]) = realPolesCoefficient[freqIndex][coeffIndex];
            }
            for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_re.size(); ++coeffIndex) {
                // A(freqIndex,complexPolesRealIndexInA(coeffIndex)) = complexPolesCoefficientReal(freqIndex,coeffIndex);
                AResidues(freqIndex,idx_A_complex_re[coeffIndex]) = complexPolesCoefficientReal[freqIndex][coeffIndex];
            }
            for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_im.size(); ++coeffIndex) {
                // A(freqIndex,complexPolesImagIndexInA(coeffIndex)) = complexPolesCoefficientImag(freqIndex,coeffIndex);
                AResidues(freqIndex,idx_A_complex_im[coeffIndex]) = complexPolesCoefficientImag[freqIndex][coeffIndex];
            }

            if (isProportional) {
                // A(freqIndex,proportionalIndexInA) = s(freqIndex);
                AResidues(freqIndex,indexProportional) = s[freqIndex];
            }
            if (isConstant) {
                // A(freqIndex,constantIndexInA) = 1;
                AResidues(freqIndex,indexConstant) = std::complex<double>(1.0, 0.0);
            }
        }









        // 1. 获取矩阵维度
        int rows = AResidues.rows();
        int cols = AResidues.cols();

        // 2. 创建结果矩阵（行数为原矩阵2倍，列数相同）
        Eigen::MatrixXd AResiduesCombineRealImag = Eigen::MatrixXd::Zero(2 * rows, cols);

        // 3. 将实部放在上半部分
        AResiduesCombineRealImag.topRows(rows) = AResidues.real();

        // 4. 将虚部放在下半部分
        AResiduesCombineRealImag.bottomRows(rows) = AResidues.imag();


        Eigen::MatrixXd freqResponseRI = Eigen::MatrixXd::Zero(frequencySize*2,freqResponse[0].size());

        for(int i=0;i<freqResponse.size();i++)
        {
            for(int j=0;j<freqResponse[0].size();j++)
            {
                freqResponseRI(i,j) = freqResponse[i][j].real();
                freqResponseRI(i+frequencySize,j) = freqResponse[i][j].imag();
            }
        }




        Eigen::VectorXd scale1 = Eigen::VectorXd::Zero(AResiduesCombineRealImag.cols());

        for(int i=0;i<AResiduesCombineRealImag.cols();i++)
        {
            scale1(i) = (AResiduesCombineRealImag.col(i).norm());
            AResiduesCombineRealImag.col(i) = AResiduesCombineRealImag.col(i) / scale1(i);
        }






        Eigen::MatrixXd xResiduse = AResiduesCombineRealImag.colPivHouseholderQr().solve(freqResponseRI);
        // Eigen::MatrixXd xResiduse = AResiduesCombineRealImag.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(freqResponseRI);



        // xResiduse = xResiduse.array() / scale1.array();
        for(int i=0;i<xResiduse.rows();i++)
        {
            for(int j=0;j<xResiduse.cols();j++)
            {
                xResiduse(i,j) = xResiduse(i,j) / scale1(i);
            }
        }


        Eigen::MatrixXcd residuseTemp = Eigen::MatrixXcd::Zero(portNum*portNum,indexRealPolse.size() + indexComplexPoles.size());

        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<(int)(idx_A_real.size());j++)
            {
                residuseTemp(i,indexRealPolse[j]) = xResiduse(idx_A_real[j],i);
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                residuseTemp(i,indexComplexPoles[j]) = std::complex<double>(xResiduse(idx_A_complex_re[j],i),xResiduse(idx_A_complex_im[j],i));
            }
        }

        if(isConstant)
        {
            Eigen::VectorXd constant_coeff = xResiduse.row(indexConstant);

            constant.clear();
            constant.resize(constant_coeff.size());
            for(int i=0;i<constant_coeff.size();i++)
            {
                constant[i] = constant_coeff(i);
            }
        }

        if(isProportional)
        {
            Eigen::VectorXd proportional_coeff = xResiduse.row(indexProportional);


            // if(isFrequencyNormalized)
            {
                for(int i=0;i<proportional_coeff.size();i++)
                {
                    proportional_coeff(i) /= frequencyAverage;
                }
            }

            proportional.clear();
            proportional.resize(proportional_coeff.size());
            for(int i=0;i<proportional_coeff.size();i++)
            {
                proportional[i] = proportional_coeff(i);
            }
        }



        residuse.clear();
        residuse = std::vector<std::vector<std::complex<double>>>(residuseTemp.rows(),std::vector<std::complex<double>>(residuseTemp.cols(),0));
        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<residuseTemp.cols();j++)
            {
                // if(isFrequencyNormalized)
                {
                    residuse[i][j] = residuseTemp(i,j)*frequencyAverage;
                }
                // else
                // {
                //     residuse[i][j] = residuseTemp(i,j);
                // }
            }
        }


        // if(isFrequencyNormalized)
        {
            for(int i=0;i<poles.size();i++)
            {
                poles[i] = (poles[i]*frequencyAverage);
            }
        }




#ifdef QT_DEBUG
        qDebug()<<QString("The %1th Iteration.").arg(10 - iterations + 1);
#endif

        std::vector<std::vector<std::vector<std::complex<double>>>> sparameterDataVF;
        getSparameter(frequency,poles,residuse,constant,proportional,sparameterDataVF);

        double currentRms = getRMS(sparameterData,sparameterDataVF);
        finalRms = currentRms;

#ifdef QT_DEBUG
        qDebug()<<QString("The %1th Iteration VF error(rms): %2").arg(10 - iterations + 1).arg(currentRms);


        for(int i=0;i<poles.size();i++)
        {
            qDebug()<<QString("Pole real: %1   Pole imag: %2").arg(poles[i].real()).arg(poles[i].imag());
        }
#endif

        if(isProportional)
        {
            // if(isFrequencyNormalized)
            {
                for(int i=0;i<proportional.size();i++)
                {
                    proportional[i] *= frequencyAverage;
                }
            }
        }

        {
            for(int i=0;i<poles.size();i++)
            {
                poles[i] = poles[i]/frequencyAverage;
            }
        }



        if(currentRms<rms)
        {
            break;
        }
        else
        {
            iterations--;
        }

    }


    if(isProportional)
    {
        // if(isFrequencyNormalized)
        {
            for(int i=0;i<proportional.size();i++)
            {
                proportional[i] /= frequencyAverage;
            }
        }
    }

    // if(isFrequencyNormalized)
    {
        for(int i=0;i<poles.size();i++)
        {
            poles[i] = (poles[i]*frequencyAverage);
        }
    }

    inputData.poles1.clear();
    inputData.poles1 = poles;
    inputData.residuse.clear();
    inputData.residuse = residuse;
    inputData.constant.clear();
    inputData.constant = constant;
    inputData.proportional.clear();
    inputData.proportional = proportional;


    return finalRms;
}

void passivityEnforceModel(int numPorts, std::vector<double> frequency,
                           std::vector<std::complex<double> > poles,
                           std::vector<std::vector<std::complex<double> > > residue,
                           std::vector<double> proportional,
                           std::vector<double> constant, bool isProportional, bool isConstant, int iterations,
                           std::vector<std::vector<std::complex<double>>>& residueRes,
                           std::vector<double>& constantRes)
{
    double tolImagEqualZero = 1e-10;  // 设置容差判断是否小于0

    std::vector<int> frequencyVioIndex;
    std::vector<std::vector<std::complex<double> > >  residueOrigin = residue;
    std::vector<double> proportionalOrigin = proportional;
    std::vector<double> constantOrigin = constant;


    double delta = 0.999;
    int currentIterationNum = 1;
    int portNum = numPorts;

    while(iterations >= 0)
    {
        std::vector<std::vector<std::vector<std::complex<double>>>> sParameterFit;
        getSparameter(frequency,poles,residue,constant,proportional,sParameterFit);

        for(int index = 0;index<sParameterFit.size();index++)
        {
            Eigen::MatrixXcd SMatrix(portNum,portNum);
            for(int i=0;i<portNum;i++)
            {
                for(int j=0;j<portNum;j++)
                {
                    SMatrix(i,j) = sParameterFit[index][i][j];
                }
            }
            Eigen::JacobiSVD<Eigen::MatrixXcd> svd(SMatrix,Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::MatrixXd S = svd.singularValues().asDiagonal();
            // 3. 获取分解因子
            Eigen::MatrixXcd U = svd.matrixU();          // 左奇异向量
            Eigen::VectorXd S_vec = svd.singularValues(); // 奇异值（降序排列的向量）
            Eigen::MatrixXcd V = svd.matrixV();          // 右奇异向量

            Eigen::MatrixXd tempS = Eigen::MatrixXd::Zero(portNum,portNum);
            if(S.maxCoeff() > delta)
            {
                for(int i=0;i<S.rows();i++)
                {
                    for(int j=0;j<S.rows();j++)
                    {
                        if(S(i,j) > delta)
                        {
                            tempS(i,j) = S(i,j) - delta;
                        }
                    }
                }
                Eigen::MatrixXcd SMatrixVio = U * tempS * V.transpose();
                for(int i=0;i<portNum;i++)
                {
                    for(int j=0;j<portNum;j++)
                    {
                        sParameterFit[index][i][j] = SMatrixVio(i,j);
                    }
                }
            frequencyVioIndex.push_back(index);
            }
            else
            {
                sParameterFit[index] = std::vector<std::vector<std::complex<double>>>(numPorts, std::vector<std::complex<double>>(numPorts,0));
            }
        }

        std::sort(frequencyVioIndex.begin(), frequencyVioIndex.end());                       // 先排序
        frequencyVioIndex.erase(std::unique(frequencyVioIndex.begin(), frequencyVioIndex.end()), frequencyVioIndex.end());   // 再去重


        sParameterFit.clear();
        std::vector<double>frequencyVio;
        for(int i=0;i<frequencyVioIndex.size();i++)
        {
            frequencyVio.push_back(frequency[frequencyVioIndex[i]]);
        }
        getSparameter(frequencyVio,poles,residue,constant,proportional,sParameterFit);
        std::vector<std::complex<double>> s;
        for(int i=0;i<frequencyVio.size();i++)
        {
            s.push_back(std::complex<double>(0,2*M_PI*frequencyVio[i]));
        }
        int frequencyNum = frequencyVioIndex.size();

        for(int index = 0;index<sParameterFit.size();index++)
        {
            Eigen::MatrixXcd SMatrix(portNum,portNum);
            for(int i=0;i<portNum;i++)
            {
                for(int j=0;j<portNum;j++)
                {
                    SMatrix(i,j) = sParameterFit[index][i][j];
                }
            }
            Eigen::JacobiSVD<Eigen::MatrixXcd> svd(SMatrix,Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::MatrixXd S = svd.singularValues().asDiagonal();
            // 3. 获取分解因子
            Eigen::MatrixXcd U = svd.matrixU();          // 左奇异向量
            Eigen::VectorXd S_vec = svd.singularValues(); // 奇异值（降序排列的向量）
            Eigen::MatrixXcd V = svd.matrixV();          // 右奇异向量

            Eigen::MatrixXd tempS = Eigen::MatrixXd::Zero(portNum,portNum);
            if(S.maxCoeff() > delta)
            {
                for(int i=0;i<S.rows();i++)
                {
                    for(int j=0;j<S.rows();j++)
                    {
                        if(S(i,j) > delta)
                        {
                            tempS(i,j) = S(i,j) - delta;
                        }
                    }
                }
                Eigen::MatrixXcd SMatrixVio = U * tempS * V.transpose();
                for(int i=0;i<portNum;i++)
                {
                    for(int j=0;j<portNum;j++)
                    {
                        sParameterFit[index][i][j] = SMatrixVio(i,j);
                    }
                }
            }
            else
            {
                sParameterFit[index] = std::vector<std::vector<std::complex<double>>>(numPorts, std::vector<std::complex<double>>(numPorts,0));
            }
        }


        std::vector<int> realIndex;
        std::vector<int> complexIndex;
        for(int i=0;i<poles.size();i++)
        {
            if(std::abs(poles[i].imag())<=tolImagEqualZero)
            {
                realIndex.push_back(i);
            }
            else
            {
                complexIndex.push_back(i);
            }
        }

        int n_cols = realIndex.size() + 2*complexIndex.size();
        int indexConstant = 0;
        if(isConstant)
        {
            indexConstant = n_cols;
            n_cols++;
        }

        int indexProportional = 0;
        if(isProportional)
        {
            indexProportional = n_cols;
            n_cols++;
        }

        std::vector<int> indexRealPolse(realIndex.size(),0);
        std::vector<int> indexComplexPoles(complexIndex.size(),0);

        for(int i=0;i<indexRealPolse.size();i++)
        {
            indexRealPolse[i] = realIndex[i];
        }

        for(int i=0;i<indexComplexPoles.size();i++)
        {
            indexComplexPoles[i] = complexIndex[i];
        }




        int realPolesNum = indexRealPolse.size();
        int complexPolesNum = indexComplexPoles.size();
        std::vector<int> idx_A_real(realPolesNum,0);

        for(int i=0;i<realPolesNum;i++)
        {
            idx_A_real[i] = i;
        }
        std::vector<int> idx_A_complex_re(complexPolesNum, 0);

        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_re[i] = realPolesNum + 2*i;
        }
        std::vector<int> idx_A_complex_im(complexPolesNum, 0);

        for(int i=0;i<complexPolesNum;i++)
        {
            idx_A_complex_im[i] = idx_A_complex_re[i] + 1;
        }

        Eigen::MatrixXcd AResidues = Eigen::MatrixXcd::Zero(frequencyNum,n_cols);



        std::vector<std::vector<std::complex<double>>> realPolesCoefficient(frequencyNum,
                                                                              std::vector<std::complex<double>>(indexRealPolse.size(), std::complex<double>(0.0, 0.0)));

        for (int freqIndex = 0; freqIndex < frequencyNum; ++freqIndex) {
            for (int index = 0; index < (int)indexRealPolse.size(); ++index) {
                // realPolesCoefficient(freqIndex,index) = 1 / (s(freqIndex) - poles(realPolesIndexInPoles(index)));
                realPolesCoefficient[freqIndex][index] = 1.0 / (s[freqIndex] - poles[indexRealPolse[index]]);
            }
        }


        std::vector<std::vector<std::complex<double>>> complexPolesCoefficientReal(frequencyNum,
                                                                                     std::vector<std::complex<double>>(indexComplexPoles.size(), std::complex<double>(0.0, 0.0)));
        std::vector<std::vector<std::complex<double>>> complexPolesCoefficientImag(frequencyNum,
                                                                                     std::vector<std::complex<double>>(indexComplexPoles.size(), std::complex<double>(0.0, 0.0)));

        for (int freqIndex = 0; freqIndex < frequencyNum; ++freqIndex) {
            for (int index = 0; index < (int)indexComplexPoles.size(); ++index) {
                std::complex<double> currentPole = poles[indexComplexPoles[index]];
                std::complex<double> currentS = s[freqIndex];
                std::complex<double> j_unit(0.0, 1.0);

                // complexPolesCoefficientReal(freqIndex,index) = 1 / (s(freqIndex) - poles(complexPolesIndexInPoles(index))) + 1 / (s(freqIndex) - conj(poles(complexPolesIndexInPoles(index)))) ;
                complexPolesCoefficientReal[freqIndex][index] = 1.0 / (currentS - currentPole) + 1.0 / (currentS - std::conj(currentPole));
                // complexPolesCoefficientImag(freqIndex,index) = 1j / (s(freqIndex) - poles(complexPolesIndexInPoles(index))) - 1j / (s(freqIndex) - conj(poles(complexPolesIndexInPoles(index)))) ;
                complexPolesCoefficientImag[freqIndex][index] = j_unit / (currentS - currentPole) - j_unit / (currentS - std::conj(currentPole));
            }
        }


        for (int freqIndex = 0; freqIndex < frequencyNum; ++freqIndex) {
            for (int coeffIndex = 0; coeffIndex < (int)idx_A_real.size(); ++coeffIndex) {
                // A(freqIndex,realPolesIndexInA(coeffIndex)) = realPolesCoefficient(freqIndex,coeffIndex);
                AResidues(freqIndex,idx_A_real[coeffIndex]) = realPolesCoefficient[freqIndex][coeffIndex];
            }
            for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_re.size(); ++coeffIndex) {
                // A(freqIndex,complexPolesRealIndexInA(coeffIndex)) = complexPolesCoefficientReal(freqIndex,coeffIndex);
                AResidues(freqIndex,idx_A_complex_re[coeffIndex]) = complexPolesCoefficientReal[freqIndex][coeffIndex];
            }
            for (int coeffIndex = 0; coeffIndex < (int)idx_A_complex_im.size(); ++coeffIndex) {
                // A(freqIndex,complexPolesImagIndexInA(coeffIndex)) = complexPolesCoefficientImag(freqIndex,coeffIndex);
                AResidues(freqIndex,idx_A_complex_im[coeffIndex]) = complexPolesCoefficientImag[freqIndex][coeffIndex];
            }

            if (isProportional) {
                // A(freqIndex,proportionalIndexInA) = s(freqIndex);
                AResidues(freqIndex,indexProportional) = s[freqIndex];
            }
            if (isConstant) {
                // A(freqIndex,constantIndexInA) = 1;
                AResidues(freqIndex,indexConstant) = std::complex<double>(1.0, 0.0);
            }
        }









        // 1. 获取矩阵维度
        int rows = AResidues.rows();
        int cols = AResidues.cols();

        // 2. 创建结果矩阵（行数为原矩阵2倍，列数相同）
        Eigen::MatrixXd AResiduesCombineRealImag = Eigen::MatrixXd::Zero(2 * rows, cols);

        // 3. 将实部放在上半部分
        AResiduesCombineRealImag.topRows(rows) = AResidues.real();

        // 4. 将虚部放在下半部分
        AResiduesCombineRealImag.bottomRows(rows) = AResidues.imag();



        std::vector<std::vector<std::complex<double>>> freqResponse(frequencyNum,
                                                                    std::vector<std::complex<double>>(portNum*portNum,0));


        for(int freq = 0;freq<frequencyNum;freq++)
        {
            for(int row = 0;row<portNum;row++)
            {
                for(int col=0;col<portNum;col++)
                {
                    freqResponse[freq][row*portNum+col] = sParameterFit[freq][row][col];
                }
            }
        }

        Eigen::MatrixXd freqResponseRI = Eigen::MatrixXd::Zero(frequencyNum*2,freqResponse[0].size());

        for(int i=0;i<freqResponse.size();i++)
        {
            for(int j=0;j<freqResponse[0].size();j++)
            {
                freqResponseRI(i,j) = freqResponse[i][j].real();
                freqResponseRI(i+frequencyNum,j) = freqResponse[i][j].imag();
            }
        }




        Eigen::VectorXd scale1 = Eigen::VectorXd::Zero(AResiduesCombineRealImag.cols());

        for(int i=0;i<AResiduesCombineRealImag.cols();i++)
        {
            scale1(i) = (AResiduesCombineRealImag.col(i).norm());
            if(std::abs(scale1(i)) <= tolImagEqualZero || std::isnan(scale1(i)) || std::isinf(scale1(i)))
            {
                scale1(i) = 1;
            }
            AResiduesCombineRealImag.col(i) = AResiduesCombineRealImag.col(i) / scale1(i);
        }






        Eigen::MatrixXd xResiduse = AResiduesCombineRealImag.colPivHouseholderQr().solve(freqResponseRI);
        // Eigen::MatrixXd xResiduse = AResiduesCombineRealImag.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(freqResponseRI);



        // xResiduse = xResiduse.array() / scale1.array();
        for(int i=0;i<xResiduse.rows();i++)
        {
            for(int j=0;j<xResiduse.cols();j++)
            {
                xResiduse(i,j) = xResiduse(i,j) / scale1(i);
            }
        }


        Eigen::MatrixXcd residuseTemp = Eigen::MatrixXcd::Zero(portNum*portNum,indexRealPolse.size() + indexComplexPoles.size());

        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<(int)(idx_A_real.size());j++)
            {
                residuseTemp(i,indexRealPolse[j]) = xResiduse(idx_A_real[j],i);
            }

            for(int j=0;j<idx_A_complex_im.size();j++)
            {
                residuseTemp(i,indexComplexPoles[j]) = std::complex<double>(xResiduse(idx_A_complex_re[j],i),xResiduse(idx_A_complex_im[j],i));
            }
        }

        if(isConstant)
        {
            Eigen::VectorXd constant_coeff = xResiduse.row(indexConstant);

            for(int i=0;i<constant_coeff.size();i++)
            {
                constant[i] = constant[i] - constant_coeff(i);
            }
        }

        if(isProportional)
        {
            Eigen::VectorXd proportional_coeff = xResiduse.row(indexProportional);
            for(int i=0;i<proportional_coeff.size();i++)
            {
                proportional[i] = proportional[i] - proportional_coeff(i);
            }
        }



        // residue.clear();
        // residue = std::vector<std::vector<std::complex<double>>>(residuseTemp.rows(),std::vector<std::complex<double>>(residuseTemp.cols(),0));
        for(int i=0;i<residuseTemp.rows();i++)
        {
            for(int j=0;j<residuseTemp.cols();j++)
            {
                // if(isFrequencyNormalized)
                {
                    residue[i][j] = residue[i][j] - residuseTemp(i,j);
                }
            }
        }


        sParameterFit.clear();
        getSparameter(frequency,poles,residue,constant,proportional,sParameterFit);


        iterations = iterations - 1;
        currentIterationNum = currentIterationNum + 1;

        std::vector<double> maxSigular;
        std::vector<int>frequencyVioIndexTemp;
        for(int index = 0;index<sParameterFit.size();index++)
        {
            Eigen::MatrixXcd SMatrix(portNum,portNum);
            for(int i=0;i<portNum;i++)
            {
                for(int j=0;j<portNum;j++)
                {
                    SMatrix(i,j) = sParameterFit[index][i][j];
                }
            }
            Eigen::JacobiSVD<Eigen::MatrixXcd> svd(SMatrix,Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::MatrixXd S = svd.singularValues().asDiagonal();

            maxSigular.push_back(S.maxCoeff());
            if(S.maxCoeff() > delta)
            {
                frequencyVioIndexTemp.push_back(index);
            }
        }

        if(*std::max_element(maxSigular.begin(),maxSigular.end()) < 1)
        {
            break;
        }


        for(int i=0;i<frequencyVioIndexTemp.size();i++)
        {
            if (!(std::find(frequencyVioIndex.begin(), frequencyVioIndex.end(), frequencyVioIndexTemp[i]) != frequencyVioIndex.end())) {
                for(int j=0;j<frequencyVioIndexTemp.size();j++)
                {
                    frequencyVioIndex.push_back(frequencyVioIndexTemp[j]);
                }

                residue = residueOrigin;
                proportional = proportionalOrigin;
                constant = constantOrigin;
            }
        }

    }

    residueRes.clear();
    residueRes = residue;
    constantRes.clear();
    constantRes = constant;


}
