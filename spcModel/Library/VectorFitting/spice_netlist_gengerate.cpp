#include "spice_netlist_gengerate.h"


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>


std::string formatter(double val) {
    if (val == 0.0) {
        return "0";
    }

    // Engineering prefixes mapping
    const std::vector<int> expMap = {-15, -12, -9, -6, -3, 0, 3, 6, 9, 12};
    const std::vector<std::string> prefixMap = {"f", "p", "n", "u", "m", "", "k", "meg", "g", "t"};

    double absVal = std::abs(val);
    int exponent = static_cast<int>(std::floor(std::log10(absVal) / 3.0) * 3.0);
    double mantissa = val / std::pow(10.0, exponent);

    auto it = std::find(expMap.begin(), expMap.end(), exponent);
    if (it != expMap.end()) {
        int idx = std::distance(expMap.begin(), it);
        char buffer[64];
        std::snprintf(buffer, sizeof(buffer), "%.3f%s", mantissa, prefixMap[idx].c_str());
        return std::string(buffer);
    } else {
        char buffer[64];
        std::snprintf(buffer, sizeof(buffer), "%.3e", val);
        return std::string(buffer);
    }
}

/**
 * Creates an equivalent N-port SPICE subcircuit based on its vector fitted S parameter responses.
 *
 * Parameters
 * ----------
 * nports : int
 *     Number of ports.
 * z0 : std::complex<double>
 *     Port reference impedance (should be resistive).
 * poles : std::vector<std::complex<double>>
 *     Vector of poles.
 * residues : std::vector<std::vector<std::complex<double>>>
 *     Matrix of residues (rows: responses, columns: poles).
 * proportional_coeff : std::vector<double>
 *     Proportional coefficients (e) for each response.
 * constant_coeff : std::vector<double>
 *     Constant coefficients (d) for each response.
 * file : std::string
 *     Path and filename for the SPICE subcircuit file.
 * fitted_model_name : std::string
 *     Name of the resulting model.
 */
void write_spice_subcircuit_s_controlled_source(int nports,
                                                std::complex<double> z0,
                                                const std::vector<std::complex<double>>& poles,
                                                const std::vector<std::vector<std::complex<double>>>& residues,
                                                const std::vector<double>& proportional_coeff,
                                                const std::vector<double>& constant_coeff,
                                                const std::string& file,
                                                std::string fitted_model_name)
{

    // counter for subcircuits for the equivalent admittances
    int subckt_count = 0;

    // provides a unique SPICE subcircuit identifier (X1, X2, X3, ...)
    auto get_new_subckt_identifier = [&subckt_count]() {
        subckt_count++;
        char buffer[32];
        std::snprintf(buffer, sizeof(buffer), "X%d", subckt_count);
        return std::string(buffer);
    };

    std::ofstream f(file);
    if (!f.is_open()) {
        throw std::runtime_error("Could not open file " + file + " for writing.");
    }

    // write title line
    f << "* EQUIVALENT CIRCUIT FOR VECTOR FITTED S-MATRIX" << std::endl;
    f << "*" << std::endl;

    // define the complete equivalent circuit as a subcircuit with one input node per port
    // those port nodes are labeled p1, p2, p3, ...
    // all ports share a common node for ground reference (node 0)
    std::string str_input_nodes = "";
    for (int n = 1; n <= nports; ++n) {
        str_input_nodes += "p" + std::to_string(n) + " ";
    }

    f << ".SUBCKT " << fitted_model_name << " " << str_input_nodes << std::endl;

    for (int n = 1; n <= nports; ++n) {
        f << "*" << std::endl;
        f << "* port " << n << std::endl;
        // add port reference impedance z0 (has to be resistive, no imaginary part)
        f << "R" << n << " a" << n << " 0 " << std::real(z0) << std::endl;

        // add dummy voltage sources (V=0) to measure the input current
        f << "V" << n << " p" << n << " a" << n << " 0" << std::endl;

        // CCVS and VCVS driving the transfer admittances with a = V/2/sqrt(Z0) + I/2*sqrt(Z0)
        // In
        f << "H" << n << " nt" << n << " nts" << n << " V" << n << " " << std::real(z0) << std::endl;
        // Vn
        f << "E" << n << " nts" << n << " 0 p" << n << " 0 " << 1 << std::endl;

        for (int j = 1; j <= nports; ++j) {
            f << "* transfer network for s" << n << j << std::endl;

            // stacking order in VectorFitting class variables:
            // s11, s12, s13, ..., s21, s22, s23, ...
            int i_response = (n - 1) * nports + (j - 1);

            // add CCCS to generate the scattered current I_nj at port n
            // control current is measured by the dummy voltage source at the transfer network Y_nj
            // the scattered current is injected into the port (source positive connected to ground)
            f << "F" << n << j << " 0 a" << n << " V" << n << j << " " << formatter(1.0 / std::real(z0)) << std::endl;
            f << "F" << n << j << "_inv a" << n << " 0 V" << n << j << "_inv " << formatter(1.0 / std::real(z0)) << std::endl;

            // add dummy voltage source (V=0) in series with Y_nj to measure current through transfer admittance
            f << "V" << n << j << " nt" << j << " nt" << n << j << " 0" << std::endl;
            f << "V" << n << j << "_inv nt" << j << " nt" << n << j << "_inv 0" << std::endl;

            // add corresponding transfer admittance Y_nj, which is modulating the control current
            // the transfer admittance is a parallel circuit (sum) of individual admittances
            f << "* transfer admittances for S" << n << j << std::endl;

            // start with proportional and constant term of the model
            // H(s) = d + s * e  model
            // Y(s) = G + s * C  equivalent admittance
            double g = 0.0;
            if (!constant_coeff.empty()) {
                g = constant_coeff[i_response];
                // add R for constant term
                if (g < 0) {
                    f << "R" << n << j << " nt" << n << j << "_inv 0 " << formatter(std::abs(1.0 / g)) << std::endl;
                } else if (g > 0) {
                    f << "R" << n << j << " nt" << n << j << " 0 " << formatter(1.0 / g) << std::endl;
                }
            }
            double c = 0.0;
            if (!proportional_coeff.empty()) {
                c = proportional_coeff[i_response];
                // add C for proportional term
                if (c < 0) {
                    f << "C" << n << j << " nt" << n << j << "_inv 0 " << formatter(std::abs(c)) << std::endl;
                } else if (c > 0) {
                    f << "C" << n << j << " nt" << n << j << " 0 " << formatter(c) << std::endl;
                }
            }





            // add pairs of poles and residues
            for (size_t i_pole = 0; i_pole < poles.size(); ++i_pole) {
                std::complex<double> pole = poles[i_pole];
                std::complex<double> residue = residues[i_response][i_pole];

                char nt_node_buffer[64];
                std::snprintf(nt_node_buffer, sizeof(nt_node_buffer), " nt%d%d", n, j);
                std::string node = get_new_subckt_identifier() + std::string(nt_node_buffer);

                if (std::real(residue) < 0.0) {
                    // multiplication with -1 required, otherwise the values for RLC would be negative
                    // this gets compensated by inverting the transfer current direction for this subcircuit
                    residue = -1.0 * residue;
                    node = node + "_inv";
                }

                if (std::imag(pole) == 0.0) {
                    // real pole; add rl_admittance
                    double l_val = 1.0 / std::real(residue);
                    double r_val = -1.0 * std::real(pole) / std::real(residue);
                    f << node << " 0 rl_admittance res=" << formatter(r_val) << " ind=" << formatter(l_val) << std::endl;
                } else {
                    // complex pole of a conjugate pair; add rcl_vccs_admittance
                    double l_val = 1.0 / (2.0 * std::real(residue));
                    double b = -2.0 * (std::real(residue) * std::real(pole) + std::imag(residue) * std::imag(pole));
                    double r_val = -1.0 * std::real(pole) / std::real(residue);
                    double c_val = 2.0 * std::real(residue) / (std::pow(std::abs(pole), 2));
                    if(isinf(l_val) || isinf(r_val) || isinf(c_val))
                    {
                        continue;
                    }
                    double gm_add = b * l_val * c_val;
                    int m;
                    if (gm_add < 0) {
                        m = -1;
                    } else {
                        m = 1;
                    }
                    f << node << " 0 rcl_vccs_admittance res=" << formatter(r_val)
                      << " cap=" << formatter(c_val) << " ind=" << formatter(l_val)
                      << " gm=" << formatter(std::abs(gm_add)) << " mult=" << m << std::endl;
                }
            }
        }
    }

    f << ".ENDS " << fitted_model_name << std::endl;

    f << "*" << std::endl;

    // subcircuit for an active RCL+VCCS equivalent admittance Y(s) of a complex-conjugate pole-residue pair H(s)
    // Residue: c = c' + j * c"
    // Pole: p = p' + j * p"
    // H(s)  = c / (s - p) + conj(c) / (s - conj(p))
    //       = (2 * c' * s - 2 * (c'p' + c"p")) / (s ** 2 - 2 * p' * s + |p| ** 2)
    // Y(S)  = (1 / L * s + b) / (s ** 2 + R / L * s + 1 / (L * C))
    f << ".SUBCKT rcl_vccs_admittance n_pos n_neg res=1k cap=1n ind=100p gm=1m mult=1" << std::endl;
    f << "L1 n_pos 1 ind" << std::endl;
    f << "C1 1 2 cap" << std::endl;
    f << "R1 2 n_neg res" << std::endl;
    f << "G1 n_pos n_neg 1 2 gm * mult" << std::endl;
    f << ".ENDS rcl_vccs_admittance" << std::endl;

    f << "*" << std::endl;

    // subcircuit for a passive RL equivalent admittance Y(s) of a real pole-residue pair H(s)
    // H(s) = c / (s - p)
    // Y(s) = 1 / L / (s + s * R / L)
    f << ".SUBCKT rl_admittance n_pos n_neg res=1k ind=100p" << std::endl;
    f << "L1 n_pos 1 ind" << std::endl;
    f << "R1 1 n_neg res" << std::endl;
    f << ".ENDS rl_admittance" << std::endl;

    f.close();
}


void write_spice_subcircuit_s_RLC(int nports,
                                  std::complex<double> z0,
                                  const std::vector<std::complex<double>>& poles,
                                  const std::vector<std::vector<std::complex<double>>>& residues,
                                  const std::vector<double>& proportional_coeff,
                                  const std::vector<double>& constant_coeff,
                                  const std::string& file,
                                  std::string fitted_model_name)
{

    // counter for subcircuits for the equivalent admittances
    int subckt_count = 0;

    // provides a unique SPICE subcircuit identifier (X1, X2, X3, ...)
    auto get_new_subckt_identifier = [&subckt_count]() {
        subckt_count++;
        char buffer[32];
        std::snprintf(buffer, sizeof(buffer), "X%d", subckt_count);
        return std::string(buffer);
    };

    std::ofstream f(file);
    if (!f.is_open()) {
        throw std::runtime_error("Could not open file " + file + " for writing.");
    }

    // write title line
    f << "* EQUIVALENT CIRCUIT FOR VECTOR FITTED S-MATRIX" << std::endl;
    f << "*" << std::endl;

    // define the complete equivalent circuit as a subcircuit with one input node per port
    // those port nodes are labeled p1, p2, p3, ...
    // all ports share a common node for ground reference (node 0)
    std::string str_input_nodes = "";
    for (int n = 1; n <= nports; ++n) {
        str_input_nodes += "p" + std::to_string(n) + " ";
    }

    f << ".SUBCKT " << fitted_model_name << " " << str_input_nodes << std::endl;

    for (int n = 1; n <= nports; ++n) {
        f << "*" << std::endl;
        f << "* port " << n << std::endl;
        // add port reference impedance z0 (has to be resistive, no imaginary part)
        f << "R" << n << " a" << n << " 0 " << std::real(z0) << std::endl;

        // add dummy voltage sources (V=0) to measure the input current
        f << "V" << n << " p" << n << " a" << n << " 0" << std::endl;

        // CCVS and VCVS driving the transfer admittances with a = V/2/sqrt(Z0) + I/2*sqrt(Z0)
        // In
        f << "H" << n << " nt" << n << " nts" << n << " V" << n << " " << std::real(z0) << std::endl;
        // Vn
        f << "E" << n << " nts" << n << " 0 p" << n << " 0 " << 1 << std::endl;

        for (int j = 1; j <= nports; ++j) {
            f << "* transfer network for s" << n << j << std::endl;

            // stacking order in VectorFitting class variables:
            // s11, s12, s13, ..., s21, s22, s23, ...
            int i_response = (n - 1) * nports + (j - 1);

            // add CCCS to generate the scattered current I_nj at port n
            // control current is measured by the dummy voltage source at the transfer network Y_nj
            // the scattered current is injected into the port (source positive connected to ground)
            f << "F" << n << j << " 0 a" << n << " V" << n << j << " " << formatter(1.0 / std::real(z0)) << std::endl;
            f << "F" << n << j << "_inv a" << n << " 0 V" << n << j << "_inv " << formatter(1.0 / std::real(z0)) << std::endl;

            // add dummy voltage source (V=0) in series with Y_nj to measure current through transfer admittance
            f << "V" << n << j << " nt" << j << " nt" << n << j << " 0" << std::endl;
            f << "V" << n << j << "_inv nt" << j << " nt" << n << j << "_inv 0" << std::endl;

            // add corresponding transfer admittance Y_nj, which is modulating the control current
            // the transfer admittance is a parallel circuit (sum) of individual admittances
            f << "* transfer admittances for S" << n << j << std::endl;

            // start with proportional and constant term of the model
            // H(s) = d + s * e  model
            // Y(s) = G + s * C  equivalent admittance
            double g = 0.0;
            if (!constant_coeff.empty()) {
                g = constant_coeff[i_response];
            }
            double c = 0.0;
            if (!proportional_coeff.empty()) {
                c = proportional_coeff[i_response];
            }

            // add R for constant term
            if (g < 0) {
                f << "R" << n << j << " nt" << n << j << "_inv 0 " << formatter(std::abs(1.0 / g)) << std::endl;
            } else if (g > 0) {
                f << "R" << n << j << " nt" << n << j << " 0 " << formatter(1.0 / g) << std::endl;
            }

            // add C for proportional term
            if (c < 0) {
                f << "C" << n << j << " nt" << n << j << "_inv 0 " << formatter(std::abs(c)) << std::endl;
            } else if (c > 0) {
                f << "C" << n << j << " nt" << n << j << " 0 " << formatter(c) << std::endl;
            }

            // add pairs of poles and residues
            for (size_t i_pole = 0; i_pole < poles.size(); ++i_pole) {
                std::complex<double> pole = poles[i_pole];
                std::complex<double> residue = residues[i_response][i_pole];

                char nt_node_buffer[64];
                std::snprintf(nt_node_buffer, sizeof(nt_node_buffer), " nt%d%d", n, j);
                std::string node = get_new_subckt_identifier() + std::string(nt_node_buffer);

                if (std::real(residue) < 0.0) {
                    // multiplication with -1 required, otherwise the values for RLC would be negative
                    // this gets compensated by inverting the transfer current direction for this subcircuit
                    residue = -1.0 * residue;
                    node = node + "_inv";
                }

                if (std::imag(pole) == 0.0) {
                    // real pole; add rl_admittance
                    double l_val = 1.0 / std::real(residue);
                    double r_val = -1.0 * std::real(pole) / std::real(residue);
                    f << node << " 0 rl_admittance res=" << formatter(r_val) << " ind=" << formatter(l_val) << std::endl;
                } else {
                    // complex pole of a conjugate pair; add rcl_vccs_admittance
                    double l_val = 1.0 / (2.0 * std::real(residue));
                    double temp1 =  2 * (std::real(residue) * std::real(pole) + std::imag(residue) * std::imag(pole)) / (2 * std::real(residue));
                    double r_val1 = 1 / (2 * std::real(residue)) * (-2*std::real(pole) +temp1);
                    double c_val = 2 * std::real(residue) / (abs(pole) * abs(pole) + (-2*std::real(pole) +temp1)*temp1);
                    double r_val2 = -1/c_val/temp1;

                    if(isinf(l_val) || isinf(r_val1) || isinf(c_val) || isinf(r_val2))
                    {
                        continue;
                    }

                    f << node << " 0 rcl_admittance res1=" << formatter(r_val1)
                      << " res2=" << formatter(r_val2) << " cap=" << formatter(c_val)
                      << " ind=" << formatter(l_val) << std::endl;

                }
            }
        }
    }

    f << ".ENDS " << fitted_model_name << std::endl;

    f << "*" << std::endl;

    // subcircuit for an active RCL+VCCS equivalent admittance Y(s) of a complex-conjugate pole-residue pair H(s)
    // Residue: c = c' + j * c"
    // Pole: p = p' + j * p"
    // H(s)  = c / (s - p) + conj(c) / (s - conj(p))
    //       = (2 * c' * s - 2 * (c'p' + c"p")) / (s ** 2 - 2 * p' * s + |p| ** 2)
    // Y(S)  = (1 / L * s + b) / (s ** 2 + R / L * s + 1 / (L * C))
    f << ".SUBCKT rcl_admittance n_pos n_neg res1=1 res2=1 cap=1 ind=1" << std::endl;
    f << "L1 n_pos 1 ind" << std::endl;
    f << "R1 1 2 res1" << std::endl;
    f << "C1 2 n_neg cap" << std::endl;
    f << "R2 2 n_neg res2" << std::endl;
    f << ".ENDS rcl_admittance" << std::endl;

    f << "*" << std::endl;

    // subcircuit for a passive RL equivalent admittance Y(s) of a real pole-residue pair H(s)
    // H(s) = c / (s - p)
    // Y(s) = 1 / L / (s + s * R / L)
    f << ".SUBCKT rl_admittance n_pos n_neg res=1 ind=1" << std::endl;
    f << "L1 n_pos 1 ind" << std::endl;
    f << "R1 1 n_neg res" << std::endl;
    f << ".ENDS rl_admittance" << std::endl;

    f.close();
}
