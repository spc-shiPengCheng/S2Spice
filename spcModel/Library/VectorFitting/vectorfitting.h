#ifndef VECTORFITTING_H
#define VECTORFITTING_H


#include <complex>
#include <vector>

#include "Eigen/Core"
#include <unsupported/Eigen/CXX11/Tensor>

#include "Tool/st_struct_file.h"



double VF(int realPolesNum,int complexPolesNum,bool isProportional,bool isConstant,
        int portNum,std::vector<double> frequency,
        std::vector<std::vector<std::vector<std::complex<double>>>> sparameterData,
        int iterations, double rms,bool isRelaxed,
        std::vector<std::complex<double>>& poles1,
        std::vector<std::vector<std::complex<double>>>& residuse,
        std::vector<double>& constant,
        std::vector<double>& proportional);


double VFOrigin(int realPolesNum,int complexPolesNum,bool isProportional,bool isConstant,
          int portNum,std::vector<double> frequency,
          std::vector<std::vector<std::vector<std::complex<double>>>> sparameterData,
          int iterations, double rms,bool isFrequencyNormalized,bool isWeight,
          std::vector<std::complex<double>>& poles1,
          std::vector<std::vector<std::complex<double>>>& residuse,
          std::vector<double>& constant,
          std::vector<double>& proportional);


double VFRelaxed(int realPolesNum,int complexPolesNum,bool isProportional,bool isConstant,
                int portNum,std::vector<double> frequency,
                std::vector<std::vector<std::vector<std::complex<double>>>> sparameterData,
                int iterations, double rms,bool isFrequencyNormalized,bool isWeight,
                std::vector<std::complex<double>>& poles1,
                std::vector<std::vector<std::complex<double>>>& residuse,
                std::vector<double>& constant,
                std::vector<double>& proportional);


double vfImprove(VectorFittingInputData& inputData);

double VFImprove(VectorFittingInputData& inputData);


void getSparameter(std::vector<double> frequency,
                   std::vector<std::complex<double>> poles,
                   std::vector<std::vector<std::complex<double>>> residuse,
                   std::vector<double> constant,
                   std::vector<double> proportional,
                   std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterData);

double getRMS(std::vector<std::vector<std::vector<std::complex<double>>>>& originData,
              std::vector<std::vector<std::vector<std::complex<double>>>>& VFData);

Eigen::VectorXd getNonZeroSolutionQR(const Eigen::MatrixXd& A, double threshold = 1e-13);

void vectorfittingInvoke(int realPolesNum,int complexPolesNum,bool isProportional,bool isConstant,
                         int portNum,std::vector<double> frequency,
                         std::vector<std::vector<std::vector<std::complex<double>>>> sparameterData,
                         int iterations, double rms,bool isRelaxed,
                         std::vector<std::complex<double>>& poles1,
                         std::vector<std::vector<std::complex<double>>>& residuse,
                         std::vector<double>& constant,
                         std::vector<double>& proportional);

bool isPassiveSparameter(std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterDataRI);

void passivityEnforce(std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterDataRI);

void computeComplexPolesNum(const std::vector<double>& frequency,
                            const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameter,
                            double magLimit,
                            double proportionalLimit,
                            int& polesNum,
                            std::vector<std::complex<double>>& poles);

void passivityEnforceModel(int numPorts,std::vector<double> frequency,
                           std::vector<std::complex<double>> poles,
                           std::vector<std::vector<std::complex<double>>> residue,
                           std::vector<double> proportional,
                           std::vector<double> constant,bool isProportional,bool isConstant,int iterations,
                           std::vector<std::vector<std::complex<double>>>& residueRes,
                           std::vector<double>& constantRes);

#endif // VECTORFITTING_H
