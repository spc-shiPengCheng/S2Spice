#pragma once


#include <vector>
#include <complex>


bool isReciprocal(const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterDataRI, double threshold = 1e-6);

bool isPassive(const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterDataRI, double threshold = 1e-6);

std::vector<double> getSingular(const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterDataRI);

bool isCausal(const std::vector<std::vector<std::vector<std::complex<double>>>>& sparameterDataRI, double threshold = 1e-6);

double calculate_cqmi_matrix(
    const std::vector<std::vector<std::vector<std::complex<double>>>>& s_matrix);

double calculate_cqmi_1port(const std::vector<std::complex<double>>& s) ;
