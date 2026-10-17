#pragma once

#include <vector>
#include <complex>


void write_spice_subcircuit_s_controlled_source(int nports,
                                                std::complex<double> z0,
                                                const std::vector<std::complex<double>>& poles,
                                                const std::vector<std::vector<std::complex<double>>>& residues,
                                                const std::vector<double>& proportional_coeff,
                                                const std::vector<double>& constant_coeff,
                                                const std::string& file,
                                                std::string fitted_model_name = "s_equivalent");

void write_spice_subcircuit_s_RLC(int nports,
                                  std::complex<double> z0,
                                  const std::vector<std::complex<double>>& poles,
                                  const std::vector<std::vector<std::complex<double>>>& residues,
                                  const std::vector<double>& proportional_coeff,
                                  const std::vector<double>& constant_coeff,
                                  const std::string& file,
                                  std::string fitted_model_name = "s_equivalent");
