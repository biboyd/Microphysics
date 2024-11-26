#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

const std::string network_name = "pynucastro-cxx";

namespace network
{

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real bion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the binding energy of the element

        if constexpr (spec == N) {
            return 0.0_rt;
        }
        else if constexpr (spec == H1) {
            return 0.0_rt;
        }
        else if constexpr (spec == He4) {
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == C12) {
            return 92.16173498399803_rt;
        }
        else if constexpr (spec == C13) {
            return 97.10804378399916_rt;
        }
        else if constexpr (spec == N13) {
            return 94.10522604799917_rt;
        }
        else if constexpr (spec == N14) {
            return 104.65860734799753_rt;
        }
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == O17) {
            return 131.76239561199873_rt;
        }
        else if constexpr (spec == O18) {
            return 139.8077658120019_rt;
        }
        else if constexpr (spec == F18) {
            return 137.36950247599816_rt;
        }
        else if constexpr (spec == F21) {
            return 162.50415677599813_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
        }
        else if constexpr (spec == Ne21) {
            return 167.40598973999658_rt;
        }
        else if constexpr (spec == Ne22) {
            return 177.77024384000106_rt;
        }
        else if constexpr (spec == Ne23) {
            return 182.97089593999772_rt;
        }
        else if constexpr (spec == Na23) {
            return 186.56435240400242_rt;
        }
        else if constexpr (spec == Na25) {
            return 202.53493510400222_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.2570479679962_rt;
        }
        else if constexpr (spec == Mg25) {
            return 205.5875680680001_rt;
        }
        else if constexpr (spec == Mg27) {
            return 223.12401426799624_rt;
        }
        else if constexpr (spec == Al27) {
            return 224.95193723199915_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real mion () {
        static_assert(spec >= 1 && spec <= NumSpec);

        constexpr amrex::Real A = NetworkProperties::aion(spec);
        constexpr amrex::Real Z = NetworkProperties::zion(spec);

        return (A - Z) * C::Legacy::m_n + Z * (C::Legacy::m_p + C::Legacy::m_e) - bion<spec>() * C::Legacy::MeV2gr;
    }

    // Legacy (non-templated) interfaces

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real bion (int spec) {
        using namespace Species;

        amrex::Real b = 0.0_rt;

        // Set the binding energy of the element
        constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                b = bion<n>();
            }
        });

        return b;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real mion (int spec) {
        using namespace Species;

        amrex::Real m = 0.0_rt;

        constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                m = mion<n>();
            }
        });

        return m;
    }
}

namespace Rates
{

    enum NetworkRates
    {
        k_n_C12_to_C13 = 1,
        k_p_C12_to_N13 = 2,
        k_He4_C12_to_O16 = 3,
        k_p_C13_to_N14 = 4,
        k_He4_N14_to_F18 = 5,
        k_n_O16_to_O17 = 6,
        k_He4_O16_to_Ne20 = 7,
        k_p_O17_to_F18 = 8,
        k_He4_O17_to_Ne21 = 9,
        k_He4_O18_to_Ne22 = 10,
        k_He4_Ne20_to_Mg24 = 11,
        k_n_Ne22_to_Ne23 = 12,
        k_p_Ne22_to_Na23 = 13,
        k_p_Na23_to_Mg24 = 14,
        k_C12_C12_to_p_Na23 = 15,
        k_C12_C12_to_He4_Ne20 = 16,
        k_He4_C13_to_n_O16 = 17,
        k_C12_O16_to_p_Al27 = 18,
        k_C12_O16_to_He4_Mg24 = 19,
        k_p_O17_to_He4_N14 = 20,
        k_He4_Ne21_to_n_Mg24 = 21,
        k_p_Ne23_to_n_Na23 = 22,
        k_p_Na23_to_He4_Ne20 = 23,
        k_p_Na25_to_n_Mg25 = 24,
        k_p_Na25_to_He4_Ne22 = 25,
        k_p_Mg27_to_n_Al27 = 26,
        k_N13_to_C13 = 27,
        k_Al27_to_Mg27 = 28,
        k_F18_to_O18 = 29,
        k_F21_to_Ne21 = 30,
        k_Mg25_to_Na25 = 31,
        k_Na23_to_Ne23 = 32,
        k_Na25_to_Mg25 = 33,
        k_Ne21_to_F21 = 34,
        k_Ne23_to_Na23 = 35,
        NumRates = k_Ne23_to_Na23
    };

    // number of reaclib rates

    const int NrateReaclib = 26;

    // number of tabular rates

    const int NrateTabular = 9;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "n_C12_to_C13",  // 1,
        "p_C12_to_N13",  // 2,
        "He4_C12_to_O16",  // 3,
        "p_C13_to_N14",  // 4,
        "He4_N14_to_F18",  // 5,
        "n_O16_to_O17",  // 6,
        "He4_O16_to_Ne20",  // 7,
        "p_O17_to_F18",  // 8,
        "He4_O17_to_Ne21",  // 9,
        "He4_O18_to_Ne22",  // 10,
        "He4_Ne20_to_Mg24",  // 11,
        "n_Ne22_to_Ne23",  // 12,
        "p_Ne22_to_Na23",  // 13,
        "p_Na23_to_Mg24",  // 14,
        "C12_C12_to_p_Na23",  // 15,
        "C12_C12_to_He4_Ne20",  // 16,
        "He4_C13_to_n_O16",  // 17,
        "C12_O16_to_p_Al27",  // 18,
        "C12_O16_to_He4_Mg24",  // 19,
        "p_O17_to_He4_N14",  // 20,
        "He4_Ne21_to_n_Mg24",  // 21,
        "p_Ne23_to_n_Na23",  // 22,
        "p_Na23_to_He4_Ne20",  // 23,
        "p_Na25_to_n_Mg25",  // 24,
        "p_Na25_to_He4_Ne22",  // 25,
        "p_Mg27_to_n_Al27",  // 26,
        "N13_to_C13",  // 27,
        "Al27_to_Mg27",  // 28,
        "F18_to_O18",  // 29,
        "F21_to_Ne21",  // 30,
        "Mg25_to_Na25",  // 31,
        "Na23_to_Ne23",  // 32,
        "Na25_to_Mg25",  // 33,
        "Ne21_to_F21",  // 34,
        "Ne23_to_Na23"  // 35,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 1;
    constexpr int N_index = 0;
    constexpr int He4_index = 2;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
