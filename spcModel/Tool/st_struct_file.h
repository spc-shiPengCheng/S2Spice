#ifndef ST_STRUCT_FILE_H
#define ST_STRUCT_FILE_H

#include <complex>
#include <QString>

const double equalZero = 1e-10;


struct snp
{
    QString filePath = "";
    QString fileName = "";
    qint64 fileSize = 0;
    std::vector<double> frequency = std::vector<double>();
    double frequencyUnit = 1;
    QString frequencyUnitLabel = "";
    std::vector<std::vector<std::vector<std::complex<double>>>> sparametrRI;
    std::vector<std::vector<std::vector<std::complex<double>>>> sparametrRIVF;
    std::vector<std::vector<std::vector<std::complex<double>>>> sparametrRINetlist;
    std::vector<std::vector<std::vector<std::complex<double>>>> yparametrRI;
    std::vector<std::vector<std::vector<std::complex<double>>>> yparametrRIVF;
    std::vector<std::vector<std::vector<std::complex<double>>>> yparametrRINetlist;
    std::vector<std::vector<std::vector<std::complex<double>>>> zparametrRI;
    std::vector<std::vector<std::vector<std::complex<double>>>> zparametrRIVF;
    std::vector<std::vector<std::vector<std::complex<double>>>> zparametrRINetlist;
    double referenceZ = 0;
    int portNum = 0;
    QString fittingError = "";
    std::vector<std::complex<double>> poles;
    std::vector<std::vector<std::complex<double>>> residuse;
    std::vector<double> constant;
    std::vector<double> proportional;

    void reset()
    {
        filePath = "";
        fileName = "";
        fileSize = 0;
        frequency.clear();
        frequencyUnit = 1;
        frequencyUnitLabel = "";
        sparametrRI.clear();
        sparametrRIVF.clear();
        sparametrRINetlist.clear();
        yparametrRI.clear();
        yparametrRIVF.clear();
        yparametrRINetlist.clear();
        zparametrRI.clear();
        zparametrRIVF.clear();
        zparametrRINetlist.clear();
        referenceZ = 0;
        portNum = 0;
        fittingError = "";
        poles.clear();
        residuse.clear();
        constant.clear();
        proportional.clear();
    }
};


struct VectorFittingInputData
{
    int portNum;
    std::vector<double> frequency;
    std::vector<std::vector<std::vector<std::complex<double>>>> sparametrRI;
    int realPolesNum = 2;
    int complexPolesNum = 2;
    double complexPolesRealCoefficient = 1e-2;
    bool isComputeComplexPoles = false;
    double dataPrecision = 1e-2;
    bool isProportional = false;
    bool isConstant = true;
    bool isNormalFrequency = true;
    int iterationNum = 10;
    double fittingError = 5e-3;
    bool isWeight = true;
    bool isDirectComputePoles = false;
    bool isEnforcePolesStable = true;
    bool isRelax = true;
    std::vector<std::complex<double>> poles1;
    std::vector<std::vector<std::complex<double>>> residuse;
    std::vector<double> constant;
    std::vector<double> proportional;
};


struct VectorFittingParameter
{
    int realPolesNum = 2;
    int complexPolesNum = 2;
    bool isProportional = false;
    bool isConstant = true;
    int iterations = 10;
    double rms = 0.005;
    bool isRelaxed = true;
    bool isFrequencyNormalized = true;
    bool isWeigth = true;
    int minOrder = 2;
    int maxOrder = 50;
};

enum DisplayParameterType
{
    S = 0,
    Y,
    Z,
};

enum DisplayDataType
{
    REAL_IMAG = 0,
    MAG_ANGLE,
    DB_ANGLE,
    SMITH,
};

enum NetlistType
{
    ControlledSource = 0,
    RLC,
};

#endif // ST_STRUCT_FILE_H
