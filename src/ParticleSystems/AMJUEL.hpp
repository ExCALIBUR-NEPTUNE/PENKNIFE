#ifndef AMJUEL_HPP
#define AMJUEL_HPP

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <reactions/reactions.hpp>
#include "../Misc/Constants.hpp"

namespace NESO::Solvers::tokamak
{
using namespace VANTAGE::Reactions;
struct norm
{
    static constexpr double time        = 1.0E-8;
    static constexpr double length      = -1;
    static constexpr double temp        = 1.0;
    static constexpr double dens        = 1e18;
    static constexpr double mass_amu    = 1.0;
    static constexpr double mass_amu_SI = constants::m_p_si * mass_amu;
    static inline const double vel      = std::sqrt((2 * constants::e) / constants::m_p_si);
    static inline const double potential_energy =
        13.6 * constants::e / (mass_amu_SI * vel * vel);
};

class AMJUEL
{
private:
    template <int n1, int n2>
    static inline const std::array<std::array<double, n2>, n1> fetch_amjuel_coeffs(
        const std::string &filename)
    {
        std::ifstream input{filename};

        std::vector<std::vector<float>> csv_rows;
        for (std::string line; std::getline(input, line);)
        {
            std::istringstream ss(std::move(line));
            std::vector<float> row;
            for (std::string value; std::getline(ss, value, ',');)
            {
                row.push_back(std::move(strtod(value.c_str(), 0)));
            }
            csv_rows.push_back(std::move(row));
        }

        std::array<std::array<double, n2>, n1> arr;
        for (int i = 0; i < n2; ++i)
        {
            for (int j = 0; j < n1; ++j)
            {
                arr[i][j] = csv_rows[i][j];
            }
        }
        return arr;
    }

    static constexpr int num_coeffs_E  = 9;
    static constexpr int num_coeffs_n  = 9;
    static constexpr int num_coeffs_np = 9;
    static constexpr int num_coeffs_T  = 9;

    // Ionisation
public:
    static inline const auto ionise_rate_data(
        const std::string &filename = "data/H.4_2.1.5.csv")
    {
        const auto h4_2_1_5_coeffs =
            fetch_amjuel_coeffs<num_coeffs_n, num_coeffs_T>(filename);

        return AMJUEL2DData<num_coeffs_T, num_coeffs_n>(
            1.0, norm::dens, norm::temp, norm::time, h4_2_1_5_coeffs);
    }

    static inline const auto ionise_energy_data(
        const std::string &filename = "data/H.10_2.1.8.csv")
    {
        const auto h10_2_1_5_coeffs =
            fetch_amjuel_coeffs<num_coeffs_np, num_coeffs_T>(filename);

        return AMJUEL2DData<num_coeffs_T, num_coeffs_np>(
            norm::mass_amu * norm::vel * norm::vel, norm::dens, norm::temp,
            norm::time, h10_2_1_5_coeffs);
    }

    // Charge Exchange
private:
    static constexpr int h1_3_1_8_num_coeffs   = 9;
    static constexpr int h1_3_1_8_num_l_coeffs = 3;
    static constexpr int h1_3_1_8_num_r_coeffs = 3;

    static constexpr std::array<double, h1_3_1_8_num_coeffs> h1_3_1_8_coeffs{
        {-3.274123792568e+01, -8.916456579806e-02, -3.016990732025e-02,
         9.205482406462e-03, 2.400266568315e-03, -1.927122311323e-03,
         3.654750340106e-04, -2.788866460622e-05, 7.422296363524e-07}};

    static constexpr std::array<double, h1_3_1_8_num_l_coeffs>
        h1_3_1_8_l_coeffs{
            {-3.294589355000e+01, -1.713112000000e-01, 0.000000000000e+00}};

    static constexpr std::array<double, h1_3_1_8_num_r_coeffs>
        h1_3_1_8_r_coeffs{
            {-5.787734011000e+01, 7.671416829000e+00, -5.208376804000e-01}};

    static constexpr double h1_3_1_8_elabmin = 0.10000e0;
    static constexpr double h1_3_1_8_elabmax = 5.00000e0;
    static constexpr double h1_3_1_8_emax    = 1.0e6;

public:
    static inline const auto cx_rate_data(
        double parent_mass, double child_mass,
        const std::string &filename = "data/H.3_3.1.8.csv")
    {
        const auto h3_3_1_8_coeffs =
            fetch_amjuel_coeffs<num_coeffs_E, num_coeffs_T>(filename);
        return AMJUEL2DDataH3<num_coeffs_T, num_coeffs_E, 2>(
            1.0, norm::dens, norm::temp / child_mass, norm::time, norm::vel,
            parent_mass, h3_3_1_8_coeffs);
    }

    static inline const auto amjuel_fit_cross_section(double reduced_mass)
    {
        static constexpr double amjuel_cs_norm = 1E-4;

        return AMJUELFitCrossSection<h1_3_1_8_num_coeffs, h1_3_1_8_num_l_coeffs,
                                     h1_3_1_8_num_r_coeffs>(
            norm::vel, amjuel_cs_norm, reduced_mass, h1_3_1_8_coeffs,
            h1_3_1_8_l_coeffs, h1_3_1_8_r_coeffs, h1_3_1_8_elabmin,
            h1_3_1_8_elabmax, h1_3_1_8_emax);
    }

    // Recombination
    static inline const auto recomb_rate_data(
        const std::string &filename = "data/H.4_2.1.8.csv")
    {
        const auto h4_2_1_8_coeffs =
            fetch_amjuel_coeffs<num_coeffs_T, num_coeffs_n>(filename);
        return AMJUEL2DData<num_coeffs_T, num_coeffs_n>(
            1.0, norm::dens, norm::temp, norm::time, h4_2_1_8_coeffs);
    }

    static inline const auto recomb_energy_data(
        const std::string &filename = "data/H.10_2.1.8.csv")
    {
        const auto h10_2_1_8_coeffs =
            fetch_amjuel_coeffs<num_coeffs_T, num_coeffs_np>(filename);
        return AMJUEL2DData<num_coeffs_T, num_coeffs_np>(
            norm::mass_amu * norm::vel * norm::vel, norm::dens, norm::temp,
            norm::time, h10_2_1_8_coeffs);
    }
};
} // namespace NESO::Solvers::tokamak
#endif