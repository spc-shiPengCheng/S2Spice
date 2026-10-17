#ifndef ST_CONVERT_SYZ_H
#define ST_CONVERT_SYZ_H

#include <complex>
#include <vector>



void s2z(std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterData,
         std::vector<std::vector<std::vector<std::complex<double>>>>& zparameterData,
         std::complex<double> referenceZ);

void s2y(const std::vector<std::vector<std::vector<std::complex<double>>>>& sParameterDataRI,
         std::vector<std::vector<std::vector<std::complex<double>>>>& yParameterDataRI,
         double referenceZ);

#endif // ST_CONVERT_SYZ_H
