#ifndef SCREEN_H
#define SCREEN_H

// these need to be defined before screen_data.H is included
#define SCREEN_METHOD_null 0
#define SCREEN_METHOD_screen5 1
#define SCREEN_METHOD_chugunov2007 2
#define SCREEN_METHOD_chugunov2009 3
#define SCREEN_METHOD_chabrier1998 4

#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <network_properties.H>
#include <microphysics_math.H>
#include <fundamental_constants.H>
#include <cmath>
#include <screen_data.H>
#include <extern_parameters.H>

using namespace amrex::literals;

#if SCREEN_METHOD == SCREEN_METHOD_null
const std::string screen_name = "null";
#elif SCREEN_METHOD == SCREEN_METHOD_screen5
const std::string screen_name = "screen5";
#elif SCREEN_METHOD == SCREEN_METHOD_chugunov2007
const std::string screen_name = "chugunov2007";
#elif SCREEN_METHOD == SCREEN_METHOD_chugunov2009
const std::string screen_name = "chugunov2009";
#elif SCREEN_METHOD == SCREEN_METHOD_chabrier1998
const std::string screen_name = "chabrier1998";
#endif

struct plasma_state_t {

    amrex::Real qlam0z;
    amrex::Real qlam0zdt;
    //amrex::Real qlam0zdd;

    amrex::Real taufac;
    amrex::Real taufacdt;

    amrex::Real aa;
    amrex::Real daadt;
    //amrex::Real daadd;

    amrex::Real temp;
    amrex::Real zbar;
    amrex::Real z2bar;
    amrex::Real n_e;
    amrex::Real gamma_e_fac;
};

inline
std::ostream& operator<< (std::ostream& o, plasma_state_t const& pstate)
{
    o << "qlam0z   = " << pstate.qlam0z << std::endl;
    o << "qlam0zdt = " << pstate.qlam0zdt << std::endl;
    o << "taufac   = " << pstate.taufac << std::endl;
    o << "taufacdt = " << pstate.taufacdt << std::endl;
    o << "aa       = " << pstate.aa << std::endl;
    o << "daadt    = " << pstate.daadt << std::endl;
    o << "temp     = " << pstate.temp << std::endl;
    o << "zbar     = " << pstate.zbar << std::endl;
    o << "z2bar    = " << pstate.z2bar << std::endl;
    o << "n_e      = " << pstate.n_e << std::endl;
    o << "gamma_e_fac = " << pstate.gamma_e_fac << std::endl;

    return o;
}

AMREX_FORCE_INLINE
void
screening_init() {

}

AMREX_FORCE_INLINE
void
screening_finalize() {

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
fill_plasma_state(plasma_state_t& state, const amrex::Real temp,
                  const amrex::Real dens, amrex::Array1D<amrex::Real, 1, NumSpec> const& y) {

    amrex::Real sum = 0.0_rt;
    for (int n = 1; n <= NumSpec; n++) {
        sum += y(n);
    }
    amrex::Real abar = 1.0_rt / sum;
    amrex::Real ytot = sum;

    sum = 0.0_rt;
    amrex::Real sum2 = 0.0_rt;
    for (int n = 1; n <= NumSpec; n++) {
        sum += zion[n-1]*y(n);
        sum2 += zion[n-1]*zion[n-1]*y(n);
    }

    // Part of Eq.6 in Itoh:1979
    // 4.248719e3 = (27*pi^2*e^4*m_u/(2*k_B*hbar^2))^(1/3)
    // the extra (1/3) to make tau -> tau/3
    const amrex::Real co2 = (1.0_rt/3.0_rt) * 4.248719e3_rt;

    amrex::Real zbar = sum * abar;
    amrex::Real z2bar = sum2 * abar;

    // ntot
    amrex::Real rr = dens * ytot;
    amrex::Real tempi = 1.0_rt / temp;
    [[maybe_unused]] amrex::Real dtempi;
    if constexpr (do_T_derivatives) {
        dtempi = -tempi * tempi;
    }

    // Part of Eq. 19 in Graboske:1973
    // pp = sqrt( \tilde{z}*(rho/u_I/T) )
    amrex::Real pp = std::sqrt(rr*tempi*(z2bar + zbar));
    [[maybe_unused]] amrex::Real dppdt;
    if constexpr (do_T_derivatives) {
        amrex::Real qq = 0.5_rt/pp *(z2bar + zbar);
        dppdt = qq*rr*dtempi;
    }

    // Part version of Eq. 19 in Graboske:1973
    state.qlam0z = 1.88e8_rt * tempi * pp;
    if constexpr (do_T_derivatives) {
        state.qlam0zdt = 1.88e8_rt * (dtempi*pp + tempi*dppdt);
    }

    // Part of Eq.6 in Itoh:1979
    state.taufac = co2 * std::cbrt(tempi);
    if constexpr (do_T_derivatives) {
        state.taufacdt = -(1.0_rt/3.0_rt) * state.taufac * tempi;
    }

    amrex::Real xni = std::cbrt(rr * zbar);

    // Part of Eq.4 in Itoh:1979
    // 2.27493e5 = e^2 / ( (3*m_u/(4pi))^(1/3) *k_B )
    state.aa = 2.27493e5_rt * tempi * xni;
    if constexpr (do_T_derivatives) {
        state.daadt = 2.27493e5_rt * dtempi * xni;
    }

    state.temp = temp;
    state.zbar = zbar;
    state.z2bar = z2bar;

    // Electron number density
    // zbar * ntot works out to sum(z[i] * n[i]), after cancelling terms
    state.n_e = zbar * rr * C::n_A;

    // precomputed part of Gamma_e, from Chugunov 2009 eq. 6
    constexpr amrex::Real gamma_e_constants =
        C::q_e*C::q_e/C::k_B * gcem::pow(4.0_rt/3.0_rt*M_PI, 1.0_rt/3.0_rt);
    state.gamma_e_fac = gamma_e_constants * std::cbrt(state.n_e);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
fill_plasma_state(plasma_state_t& state, const amrex::Real temp,
                  const amrex::Real dens, amrex::Array1D<amrex::Real, 1, NumSpec> const& y) {
    if (integrator_rp::jacobian == 1) {
        constexpr int do_T_derivatives = 1;
        fill_plasma_state<do_T_derivatives>(state, temp, dens, y);
    } else {
        constexpr int do_T_derivatives = 0;
        fill_plasma_state<do_T_derivatives>(state, temp, dens, y);
    }
}

#if SCREEN_METHOD == SCREEN_METHOD_screen5
template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_screen5 (const plasma_state_t& state,
                     const scrn::screen_factors_t& scn_fac,
                     amrex::Real& scor, amrex::Real& scordt)
{
    // this subroutine calculates screening factors and their derivatives
    // for nuclear reaction rates in the weak, intermediate and strong regimes.
    // based on graboske, dewit, grossman and cooper apj 181 457 1973 for
    // weak screening. based on alastuey and jancovici apj 226 1034 1978,
    // with plasma parameters from itoh et al apj 234 1079 1979, for strong
    // screening.

    // input:
    // state   = plasma state (T, rho, abar, zbar, etc.)
    // scn_fac = screening factors for A and Z

    // output:
    // scor    = screening correction
    // scordt  = derivative of screening correction with temperature


    // fact = 2^(1/3)
    const amrex::Real fact    = 1.25992104989487e0_rt;
    const amrex::Real gamefx  = 0.3e0_rt;          // lower gamma limit for intermediate screening
    const amrex::Real gamefs  = 0.8e0_rt;          // upper gamma limit for intermediate screening
    const amrex::Real h12_max = 300.e0_rt;

    // Get the ion data based on the input index
    amrex::Real z1 = scn_fac.z1;
    amrex::Real z2 = scn_fac.z2;

    // calculate individual screening factors
    amrex::Real bb = z1 * z2;
    amrex::Real gamp = state.aa;
    [[maybe_unused]] amrex::Real gampdt;
    if constexpr (do_T_derivatives) {
        gampdt = state.daadt;
    }

    // In Eq.4 in Itoh:1979, this term is 2*Z_1*Z_2/(Z_1^(1/3) + Z_2^(1/3))
    // However here we follow Wallace:1982 Eq. A13, which is Z_1*Z_2*(2/(Z_1+Z_2))^(1/3)

    amrex::Real qq = fact * bb * scn_fac.zs13inv;

    // Full Equation of Wallace:1982 Eq. A13

    amrex::Real gamef = qq * gamp;
    [[maybe_unused]] amrex::Real gamefdt;
    if constexpr (do_T_derivatives) {
        gamefdt = qq * gampdt;
    }

    // Full version of Eq.6 in Itoh:1979 with extra 1/3 factor
    // the extra 1/3 factor is there for convenience.
    // tau12 = Eq.6 / 3

    amrex::Real tau12 = state.taufac * scn_fac.aznut;
    [[maybe_unused]] amrex::Real tau12dt;
    if constexpr (do_T_derivatives) {
        tau12dt = state.taufacdt * scn_fac.aznut;
    }

    qq = 1.0_rt/tau12;

    // alph12 = 3*gamma_ij/tau_ij

    amrex::Real alph12 = gamef * qq;
    [[maybe_unused]] amrex::Real alph12dt;
    if constexpr (do_T_derivatives) {
        alph12dt = (gamefdt - alph12*tau12dt) * qq;
    }


    // limit alph12 to 1.6 to prevent unphysical behavior.
    // See Introduction in Alastuey:1978

    // this should really be replaced by a pycnonuclear reaction rate formula
    if (alph12 > 1.6_rt) {
        alph12   = 1.6e0_rt;
        if constexpr (do_T_derivatives) {
            alph12dt = 0.0_rt;
        }

        // redetermine previous factors if 3*gamma_ij/tau_ij > 1.6

        gamef    = 1.6e0_rt * tau12;
        if constexpr (do_T_derivatives) {
            gamefdt  = 1.6e0_rt * tau12dt;
        }

        qq = scn_fac.zs13/(fact * bb);
        gamp = gamef * qq;
        if constexpr (do_T_derivatives) {
            gampdt = gamefdt * qq;
        }
    }

    // weak screening regime
    // Full version of Eq. 19 in Graboske:1973 by considering weak regime
    // and Wallace:1982 Eq. A14. Here the degeneracy factor is assumed to be 1.

    amrex::Real h12w = bb * state.qlam0z;
    [[maybe_unused]] amrex::Real dh12wdt;
    if constexpr (do_T_derivatives) {
        dh12wdt = bb * state.qlam0zdt;
    }

    amrex::Real h12 = h12w;
    [[maybe_unused]] amrex::Real dh12dt;
    if constexpr (do_T_derivatives) {
        dh12dt = dh12wdt;
    }

    // intermediate and strong sceening regime

    if (gamef > gamefx) {

      // gamma_ij^(1/4)

        amrex::Real gamp14 = std::pow(gamp, 0.25_rt);
        amrex::Real rr = 1.0_rt/gamp;

        // Here we follow Eq. A9 in Wallace:1982
        // See Eq. 25 Alastuey:1978, Eq. 16 and 17 in Jancovici:1977 for reference
        amrex::Real cc = 0.896434e0_rt * gamp * scn_fac.zhat
            - 3.44740e0_rt * gamp14 * scn_fac.zhat2
            - 0.5551e0_rt * (std::log(gamp) + scn_fac.lzav)
            - 2.996e0_rt;

        [[maybe_unused]] amrex::Real dccdt;
        if constexpr (do_T_derivatives) {
            qq = 0.25_rt * gamp14 * rr;
            amrex::Real gamp14dt = qq * gampdt;
            dccdt = 0.896434e0_rt * gampdt * scn_fac.zhat
                - 3.44740e0_rt * gamp14dt * scn_fac.zhat2
                - 0.5551e0_rt *rr * gampdt;
        }

        // (3gamma_ij/tau_ij)^3
        amrex::Real a3 = alph12 * alph12 * alph12;
        amrex::Real da3 = 3.0e0_rt * alph12 * alph12;

        // Part of Eq. 28 in Alastuey:1978
        qq = 0.014e0_rt + 0.0128e0_rt*alph12;

        // Part of Eq. 28 in Alastuey:1978
        rr = (5.0_rt/32.0_rt) - alph12*qq;
        [[maybe_unused]] amrex::Real drrdt;
        if constexpr (do_T_derivatives) {
            amrex::Real dqqdt  = 0.0128e0_rt*alph12dt;
            drrdt  = -(alph12dt*qq + alph12*dqqdt);
        }

        // Part of Eq. 28 in Alastuey:1978
        amrex::Real ss = tau12*rr;

        // Part of Eq. 31 in Alastuey:1978
        amrex::Real tt = -0.0098e0_rt + 0.0048e0_rt*alph12;

        // Part of Eq. 31 in Alastuey:1978
        amrex::Real uu = 0.0055e0_rt + alph12*tt;

        // Part of Eq. 31 in Alastuey:1978
        amrex::Real vv = gamef * alph12 * uu;

        // Exponent of Eq. 32 in Alastuey:1978, which uses Eq.28 and Eq.31
        // Strong screening factor
        h12 = cc - a3 * (ss + vv);
        if constexpr (do_T_derivatives) {
            amrex::Real dssdt  = tau12dt*rr + tau12*drrdt;
            amrex::Real dttdt  = 0.0048e0_rt*alph12dt;
            amrex::Real duudt  = alph12dt*tt + alph12*dttdt;
            amrex::Real dvvdt = gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt;
            rr = da3 * (ss + vv);
            dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt);
        }

        // See conclusion and Eq. 34 in Alastuey:1978
        // This is an extra factor to account for quantum effects
        rr = 1.0_rt - 0.0562e0_rt*a3;
        if constexpr (do_T_derivatives) {
            ss = -0.0562e0_rt*da3;
            drrdt = ss*alph12dt;
        }

        amrex::Real xlgfac;
        [[maybe_unused]] amrex::Real dxlgfacdt;

        // In extreme case, rr is 0.77, see conclusion in Alastuey:1978
        if (rr >= 0.77e0_rt) {
            xlgfac = rr;
            if constexpr (do_T_derivatives) {
                dxlgfacdt = drrdt;
            }
        } else {
            xlgfac = 0.77e0_rt;
            if constexpr (do_T_derivatives) {
                dxlgfacdt = 0.0_rt;
            }
        }

        // Include the extra factor that accounts for quantum effects
        h12 = std::log(xlgfac) + h12;
        rr = 1.0_rt/xlgfac;
        if constexpr (do_T_derivatives) {
            dh12dt = rr*dxlgfacdt + dh12dt;
        }

        // If gamma_ij < upper limit of intermediate regime
        // then it is in the intermediate regime, else strong screening.
        if (gamef <= gamefs) {
            amrex::Real dgamma = 1.0e0_rt/(gamefs - gamefx);

            rr =  dgamma*(gamefs - gamef);

            ss = dgamma*(gamef - gamefx);

            vv = h12;

            // Then the screening factor is a combination
            // of the strong and weak screening factor.
            h12 = h12w*rr + vv*ss;
            if constexpr (do_T_derivatives) {
                drrdt  = -dgamma*gamefdt;
                amrex::Real dssdt = dgamma*gamefdt;
                dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt;
            }
        }

        // end of intermediate and strong screening
    }

    // machine limit the output
    // further limit to avoid the pycnonuclear regime
    h12 = amrex::max(amrex::min(h12, h12_max), 0.0_rt);
    scor = std::exp(h12);

    if constexpr (do_T_derivatives) {
        if (h12 == h12_max) {
            scordt = 0.0;
        } else {
            scordt = scor * dh12dt;
        }
    }
}

#elif SCREEN_METHOD == SCREEN_METHOD_chugunov2007
template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void chugunov2007 (const plasma_state_t& state,
                   const scrn::screen_factors_t& scn_fac,
                   amrex::Real& scor, amrex::Real& scordt)
{
    // Calculates screening factors based on Chugunov et al. 2007, following the
    // the approach in Yakovlev 2006 to extend to a multi-component plasma.
    //
    // References:
    //   Chugunov, DeWitt, and Yakovlev 2007, PhRvD, 76, 025028
    //   Yakovlev, Gasques, Afanasjev, Beard, and Wiescher 2006, PhRvC, 74, 035803
    //   Chugunov and DeWitt 2009, PhRvC, 80, 014611

    // input:
    // state   = plasma state (T, rho, abar, zbar, etc.)
    // scn_fac = screening factors for A and Z

    // output:
    // scor    = screening correction
    // scordt  = derivative of screening correction with temperature

    amrex::Real tmp;
    // Plasma temperature T_p
    // This formula comes from working backwards from zeta_ij (Chugunov 2009 eq. 12)
    // through Chugunov 2007 eq. 3 to Chugunov 2007 eq. 2.
    // Ultimately, the changes from the expression in Chugunov 2007 are:
    //   Z^2 -> Z1 * Z2
    //   n_i -> n_e / ztilde^3, where ztilde = (Z1^(1/3) + Z2^(1/3)) / 2
    //   m_i -> 2 mu12 (reduced mass)
    // This prescription reduces to the expressions from Chugunov 2007 in the case
    // of an OCP, and to Chugunov 2009 in the case of a binary ionic mixture.
    // This also matches Yakovlev et al. 2006, eq. 10.
    //
    // For reference, MESA r21.12.1 does:
    //   Z^2 -> Z1 * Z2
    //   n_i -> n_e / zbar (=ntot)
    //   m_i -> m_u * abar
    // Sam Jones' Fortran implementation (linked in issue #146) does:
    //   Z^2 -> zbar^2
    //   n_i -> ntot
    //   m_i -> m_u * abar
    amrex::Real mu12 = scn_fac.a1 * scn_fac.a2 / (scn_fac.a1 + scn_fac.a2);
    amrex::Real z_factor = scn_fac.z1 * scn_fac.z2;
    amrex::Real n_i = state.n_e / scn_fac.ztilde3;
    amrex::Real m_i = 2.0_rt * mu12 / C::n_A;

    constexpr amrex::Real T_p_factor = C::hbar/C::k_B*C::q_e*gcem::sqrt(4.0_rt*GCEM_PI);
    amrex::Real T_p = T_p_factor * std::sqrt(z_factor * n_i / m_i);

    // Normalized temperature
    amrex::Real inv_T_p = 1.0_rt / T_p;
    amrex::Real T_norm = state.temp * inv_T_p;
    amrex::Real dT_norm_dT = inv_T_p;

    // The fit has only been verified down to T ~ 0.1 T_p, below which the rate
    // should be nearly temperature-independent (in the pycnonuclear regime),
    // and we clip the temperature to 0.1 T_p at small T.
    // start the transition here
    constexpr amrex::Real T_norm_fade = 0.2_rt;
    constexpr amrex::Real T_norm_min = 0.1_rt;

    if (T_norm < T_norm_min) {
        // clip temperature to the minimum value
        T_norm = T_norm_min;
        dT_norm_dT = 0.0_rt;
    } else if (T_norm <= T_norm_fade) {
        // blend using a cosine, after MESA
        constexpr amrex::Real delta_T = T_norm_fade - T_norm_min;
        tmp = M_PI * (T_norm - T_norm_min) / delta_T;
        amrex::Real f = 0.5_rt * (1.0_rt - std::cos(tmp));
        if constexpr (do_T_derivatives) {
            amrex::Real df_dT = 0.5_rt * M_PI / delta_T * std::sin(tmp) * dT_norm_dT;
            dT_norm_dT = -df_dT * T_norm_min + df_dT * T_norm + f * dT_norm_dT;
        }
        T_norm = (1.0_rt - f) * T_norm_min + f * T_norm;
    }

    // Coulomb coupling parameter from Yakovlev 2006 eq. 10
    amrex::Real Gamma = state.gamma_e_fac*scn_fac.z1*scn_fac.z2 / (scn_fac.ztilde*T_norm*T_p);
    [[maybe_unused]] amrex::Real dGamma_dT;
    if constexpr (do_T_derivatives) {
        dGamma_dT = -Gamma / T_norm * dT_norm_dT;
    }

    // The fit for Gamma is only applicable up to ~600, so smoothly cap its value
    constexpr amrex::Real Gamma_fade = 590;
    constexpr amrex::Real Gamma_max = 600;

    if (Gamma > Gamma_max) {
        // clip Gamma to the max value
        Gamma = Gamma_max;
        dGamma_dT = 0.0_rt;
    } else if (Gamma >= Gamma_fade) {
        // blend using a cosine, after MESA
        constexpr amrex::Real delta_gamma = Gamma_max - Gamma_fade;
        tmp = M_PI * (Gamma - Gamma_fade) / delta_gamma;
        amrex::Real f = 0.5_rt * (1.0_rt - std::cos(tmp));
        if constexpr (do_T_derivatives) {
            amrex::Real df_dT = 0.5_rt * M_PI / delta_gamma * std::sin(tmp) * dGamma_dT;
            dGamma_dT = dGamma_dT - (df_dT * Gamma + f * dGamma_dT) + df_dT * Gamma_max;
        }
        Gamma = (1.0_rt - f) * Gamma + f * Gamma_max;
    }

    // Chugunov 2007 eq. 3
    amrex::Real zeta = std::cbrt(4.0_rt / (3.0_rt * GCEM_PI*GCEM_PI * T_norm*T_norm));
    [[maybe_unused]] amrex::Real dzeta_dT;
    if constexpr (do_T_derivatives) {
        dzeta_dT = -2.0_rt / (3.0_rt * T_norm) * zeta * dT_norm_dT;
    }

    // Gamma tilde from Chugunov 2007 eq. 21
    constexpr amrex::Real fit_alpha = 0.022_rt;
    amrex::Real fit_beta = 0.41_rt - 0.6_rt / Gamma;
    amrex::Real fit_gamma = 0.06_rt + 2.2_rt / Gamma;
    // Polynomial term in Gamma tilde
    amrex::Real poly = 1.0_rt + zeta*(fit_alpha + zeta*(fit_beta + fit_gamma*zeta));
    [[maybe_unused]] amrex::Real dpoly_dT;
    if constexpr (do_T_derivatives) {
        tmp = dGamma_dT / (Gamma * Gamma);
        amrex::Real dfit_beta_dT = 0.6_rt * tmp;
        amrex::Real dfit_gamma_dT = -2.2_rt * tmp;
        dpoly_dT = (fit_alpha + 2.0_rt*zeta*fit_beta + 3.0_rt*fit_gamma*zeta*zeta) * dzeta_dT
                   + zeta*zeta*(dfit_beta_dT + dfit_gamma_dT*zeta);
    }

    amrex::Real gamtilde = Gamma / std::cbrt(poly);
    // this is gamtilde * dlog(gamtilde)/dT
    [[maybe_unused]] amrex::Real dgamtilde_dT;
    if constexpr (do_T_derivatives) {
        dgamtilde_dT = gamtilde * (dGamma_dT / Gamma - dpoly_dT / poly / 3.0_rt);
    }

    // fit parameters just after Chugunov 2007 eq. 19
    constexpr amrex::Real A1 = 2.7822_rt;
    constexpr amrex::Real A2 = 98.34_rt;
    constexpr amrex::Real A3 = gcem::sqrt(3.0_rt) - A1 / gcem::sqrt(A2);
    const amrex::Real B1 = -1.7476_rt;
    const amrex::Real B2 = 66.07_rt;
    const amrex::Real B3 = 1.12_rt;
    const amrex::Real B4 = 65_rt;
    amrex::Real gamtilde2 = gamtilde * gamtilde;

    amrex::Real term1, term2, term3, term4;
    [[maybe_unused]] amrex::Real dterm1_dT, dterm2_dT, dterm3_dT, dterm4_dT;
    // Chugunov 2007 eq. 19
    term1 = 1.0_rt / std::sqrt(A2 + gamtilde);
    if constexpr (do_T_derivatives) {
        dterm1_dT = -0.5_rt * term1 / (A2 + gamtilde) * dgamtilde_dT;
    }

    term2 = 1.0_rt / (1.0_rt + gamtilde);
    if constexpr (do_T_derivatives) {
        dterm2_dT = -term2 * term2 * dgamtilde_dT;
    }

    tmp = B2 + gamtilde;
    term3 = gamtilde2 / tmp;
    if constexpr (do_T_derivatives) {
        dterm3_dT = gamtilde * (B2 + tmp) / (tmp * tmp) * dgamtilde_dT;
    }

    tmp = B4 + gamtilde2;
    term4 = gamtilde2 / tmp;
    if constexpr (do_T_derivatives) {
        dterm4_dT = B4 / (tmp * tmp) * 2.0_rt * gamtilde * dgamtilde_dT;
    }

    amrex::Real inner = A1 * term1 + A3 * term2;

    amrex::Real gamtilde32 = std::pow(gamtilde, 1.5_rt);
    amrex::Real h = gamtilde32 * inner + B1 * term3 + B3 * term4;
    [[maybe_unused]] amrex::Real dh_dT;
    if constexpr (do_T_derivatives) {
        amrex::Real dinner_dT = A1 * dterm1_dT + A3 * dterm2_dT;
        amrex::Real dgamtilde32_dT = 1.5_rt * std::sqrt(gamtilde) * dgamtilde_dT;
        dh_dT = dgamtilde32_dT * inner + gamtilde32 * dinner_dT
                + B1 * dterm3_dT + B3 * dterm4_dT;
    }

    // machine limit the output
    constexpr amrex::Real h_max = 300.e0_rt;
    h = amrex::min(h, h_max);
    scor = std::exp(h);

    if constexpr (do_T_derivatives) {
        if (h == h_max) {
            scordt = 0.0;
        } else {
            scordt = scor * dh_dT;
        }
    }
}

#elif SCREEN_METHOD == SCREEN_METHOD_chugunov2009
template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void chugunov2009_f0 (const amrex::Real gamma, const amrex::Real dlog_dT, amrex::Real& f, amrex::Real& df_dT)
{
    // Calculates the free energy per ion in a OCP, from Chugunov and DeWitt 2009
    // equation 24.

    // input:
    // gamma
    // dlog_dT = derivative of log(gamma) with temperature

    // output:
    // f       = free energy
    // df_dT   = derivative of free energy with temperature

    // fit parameters
    constexpr amrex::Real A1 = -0.907_rt;
    constexpr amrex::Real A2 = 0.62954_rt;
    constexpr amrex::Real A3 = -gcem::sqrt(3.0_rt) / 2.0_rt - A1 / gcem::sqrt(A2);
    constexpr amrex::Real B1 = 0.00456_rt;
    constexpr amrex::Real B2 = 211.6_rt;
    constexpr amrex::Real B3 = -1e-4_rt;
    constexpr amrex::Real B4 = 0.00462_rt;
    amrex::Real gamma_12 = std::sqrt(gamma);

    amrex::Real term1, term2, term3, term4, term5;
    amrex::Real dterm1_dgamma, dterm2_dgamma, dterm3_dgamma, dterm4_dgamma, dterm5_dgamma;

    term1 = gamma_12 * std::sqrt(A2 + gamma);
    if constexpr (do_T_derivatives) {
        dterm1_dgamma = (0.5_rt * A2 + gamma) / term1;
    }

    term2 = std::log(std::sqrt(gamma / A2) + std::sqrt(1.0_rt + gamma / A2));
    if constexpr (do_T_derivatives) {
        dterm2_dgamma = 1.0_rt / (2.0_rt * term1);
    }

    term3 = gamma_12 - std::atan(gamma_12);
    if constexpr (do_T_derivatives) {
        dterm3_dgamma = 0.5_rt * gamma_12 / (1.0_rt + gamma);
    }

    term4 = std::log1p(gamma / B2);
    if constexpr (do_T_derivatives) {
        dterm4_dgamma = 1.0_rt / (B2 + gamma);
    }

    term5 = std::log1p(gamma * gamma / B4);
    if constexpr (do_T_derivatives) {
        dterm5_dgamma = 2.0_rt * gamma / (B4 + gamma * gamma);
    }

    f = A1*(term1 - A2*term2) + 2.0_rt*A3*term3 + B1*(gamma - B2*term4) + 0.5_rt*B3*term5;
    if constexpr (do_T_derivatives) {
        df_dT = (A1*(dterm1_dgamma - A2*dterm2_dgamma) +
                 2.0_rt*A3*dterm3_dgamma +
                 B1*(1.0_rt - B2*dterm4_dgamma) +
                 0.5_rt*B3*dterm5_dgamma) * dlog_dT * gamma;
    }
}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void chugunov2009 (const plasma_state_t& state,
                   const scrn::screen_factors_t& scn_fac,
                   amrex::Real& scor, amrex::Real& scordt)
{
    // Calculates screening factors based on Chugunov and DeWitt 2009, PhRvC, 80, 014611

    // input:
    // state   = plasma state (T, rho, abar, zbar, etc.)
    // scn_fac = screening factors for A and Z

    // output:
    // scor    = screening correction
    // scordt  = derivative of screening correction with temperature

    amrex::Real z1z2 = scn_fac.z1 * scn_fac.z2;
    amrex::Real zcomp = scn_fac.z1 + scn_fac.z2;

    // Gamma_e from eq. 6
    amrex::Real Gamma_e = state.gamma_e_fac / state.temp;
    // work in terms of log derivatives, since it's much simpler
    amrex::Real dlog_Gamma_dT = -1.0_rt / state.temp;

    // Coulomb coupling parameters for ions and compound nucleus, eqs. 7 & 9
    amrex::Real Gamma_1 = Gamma_e * scn_fac.z1_53;
    amrex::Real Gamma_2 = Gamma_e * scn_fac.z2_53;
    amrex::Real Gamma_comp = Gamma_e * scn_fac.zs53;

    amrex::Real Gamma_12 = Gamma_e * z1z2 / scn_fac.ztilde;

    // Coulomb barrier penetrability, eq. 10

    constexpr amrex::Real tau_factor = gcem::pow(
        27.0_rt/2.0_rt * amrex::Math::powi<2>(M_PI*C::q_e*C::q_e/C::hbar)
        / (C::n_A*C::k_B), 1.0_rt/3.0_rt);

    amrex::Real tau_12 = tau_factor * scn_fac.aznut / std::cbrt(state.temp);
    amrex::Real dlog_tau_12_dT;
    if constexpr (do_T_derivatives) {
        dlog_tau_12_dT = -1.0_rt / state.temp / 3.0_rt;
    }

    // eq. 12
    amrex::Real zeta = 3.0_rt * Gamma_12 / tau_12;
    amrex::Real dzeta_dT;
    if constexpr (do_T_derivatives) {
        dzeta_dT = zeta * (dlog_Gamma_dT - dlog_tau_12_dT);
    }

    // additional fit parameters, eq. 25
    amrex::Real y_12 = 4.0_rt * z1z2 / (zcomp * zcomp);
    amrex::Real c1 = 0.013_rt * y_12 * y_12;
    amrex::Real c2 = 0.406_rt * std::pow(y_12, 0.14_rt);
    amrex::Real c3 = 0.062_rt * std::pow(y_12, 0.19_rt) + 1.8_rt / Gamma_12;

    amrex::Real poly = 1.0_rt + zeta*(c1 + zeta*(c2 + c3*zeta));
    amrex::Real t_12 = std::cbrt(poly);
    amrex::Real dlog_dT = 0.0_rt;
    if constexpr (do_T_derivatives) {
        amrex::Real dc3_dT = -1.8_rt / Gamma_12 * dlog_Gamma_dT;
        amrex::Real dpoly_dT = (c1 + zeta*(2.0_rt*c2 + 3.0_rt*c3*zeta))*dzeta_dT +
                        dc3_dT*zeta*zeta*zeta;
        amrex::Real dlog_t_12_dT = dpoly_dT / (3.0_rt * poly);
        dlog_dT = dlog_Gamma_dT - dlog_t_12_dT;
    }

    // strong screening enhancement factor, eq. 23, replacing tau_ij with t_ij
    // Using Gamma/tau_ij gives extremely low values, while Gamma/t_ij gives
    // values similar to those from Chugunov 2007.
    amrex::Real term1, term2, term3;
    amrex::Real dterm1_dT = 0.0_rt, dterm2_dT = 0.0_rt, dterm3_dT = 0.0_rt;
    chugunov2009_f0<do_T_derivatives>(Gamma_1 / t_12, dlog_dT, term1, dterm1_dT);
    chugunov2009_f0<do_T_derivatives>(Gamma_2 / t_12, dlog_dT, term2, dterm2_dT);
    chugunov2009_f0<do_T_derivatives>(Gamma_comp / t_12, dlog_dT, term3, dterm3_dT);
    amrex::Real h_fit = term1 + term2 - term3;
    amrex::Real dh_fit_dT;
    if constexpr (do_T_derivatives) {
        dh_fit_dT = dterm1_dT + dterm2_dT - dterm3_dT;
    }

    // weak screening correction term, eq. A3
    amrex::Real corr_C = 3.0_rt*z1z2 * std::sqrt(state.z2bar/state.zbar) /
                  (scn_fac.zs52 - scn_fac.z1_52 - scn_fac.z2_52);

    // corrected enhancement factor, eq. A4
    amrex::Real Gamma_12_2 = Gamma_12 * Gamma_12;
    amrex::Real numer = corr_C + Gamma_12_2;
    amrex::Real denom = 1.0_rt + Gamma_12_2;
    amrex::Real h12 = numer / denom * h_fit;
    amrex::Real dh12_dT;
    if constexpr (do_T_derivatives) {
        amrex::Real dGamma_12_2_dT = 2 * Gamma_12_2 * dlog_Gamma_dT;
        dh12_dT = h12 * (dGamma_12_2_dT/numer - dGamma_12_2_dT/denom + dh_fit_dT/h_fit);
    }

    // machine limit the output
    const amrex::Real h12_max = 300.e0_rt;
    h12 = amrex::min(h12, h12_max);
    scor = std::exp(h12);

    if constexpr (do_T_derivatives) {
        if (h12 == h12_max) {
            scordt = 0.0_rt;
        } else {
            scordt = scor * dh12_dT;
        }
    }
}

#elif SCREEN_METHOD == SCREEN_METHOD_chabrier1998
template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void chabrier1998_helmholtz_F(const amrex::Real gamma, const amrex::Real dgamma_dT,
                              amrex::Real& f, amrex::Real& df_dT) {
    // Helmholtz free energy, See Chabrier & Potekhin 1998 Eq. 28

    // Fitted parameters, see Chabrier & Potekhin 1998 Sec.IV
    constexpr amrex::Real A_1 = -0.9052_rt;
    constexpr amrex::Real A_2 = 0.6322_rt;
    constexpr amrex::Real A_3 = -0.5_rt * gcem::sqrt(3.0_rt) - A_1 / gcem::sqrt(A_2);

    // Precompute some expressions that are reused in the derivative
    const amrex::Real sqrt_gamma = std::sqrt(gamma);
    const amrex::Real sqrt_1_gamma_A2 = std::sqrt(1.0_rt + gamma/A_2);
    const amrex::Real sqrt_gamma_A2_gamma = std::sqrt(gamma * (A_2 + gamma));
    const amrex::Real sqrt_gamma_A2 = std::sqrt(gamma/A_2);

    f = A_1 * (sqrt_gamma_A2_gamma -
               A_2 * std::log(sqrt_gamma_A2 + sqrt_1_gamma_A2)) +
               2.0_rt * A_3 * (sqrt_gamma - std::atan(sqrt_gamma));

    if constexpr (do_T_derivatives) {
        df_dT = A_1 * ((A_2 + 2.0_rt * gamma) / (2.0_rt * sqrt_gamma_A2_gamma) -
                       0.5_rt / (sqrt_gamma_A2 + sqrt_1_gamma_A2) *
                       (1.0_rt / sqrt_gamma_A2 + 1.0_rt / sqrt_1_gamma_A2)) +
                       A_3 / sqrt_gamma * (1.0_rt - 1.0_rt / (1.0_rt + gamma));

        df_dT *= dgamma_dT;
    }
}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void chabrier1998 (const plasma_state_t& state,
                   const scrn::screen_factors_t& scn_fac,
                   amrex::Real& scor, amrex::Real& scordt)
{
    // Calculates screening factors based on Chabrier & Potekhin 1998,
    // Calder2007 and partly screen5 routine mentioned in Alastuey 1978.

    // This screening is valid for weak screening: Gamma < 0.1
    // and strong screening: 1 <= Gamma <= 160
    // Reference:
    //            Chabrier & Potekhin 1998, PhRvE, 58, 4941
    //            Calder 2007, doi:10.1086/510709
    //            Wallace & Woosley 1982
    //            Alastuey 1978

    // Eq. 2 in Chabrier & Potekhin 1998

    amrex::Real Gamma_e = state.gamma_e_fac / state.temp;

    // See Calder2007 appendix Eq. A9

    amrex::Real Gamma1 = Gamma_e * scn_fac.z1_53;
    amrex::Real Gamma2 = Gamma_e * scn_fac.z2_53;
    amrex::Real Gamma12 = Gamma_e * scn_fac.zs53;

    amrex::Real Gamma1dT{}, Gamma2dT{}, Gamma12dT{};

    if constexpr (do_T_derivatives) {
        Gamma1dT = -Gamma1 / state.temp;
        Gamma2dT = -Gamma2 / state.temp;
        Gamma12dT = -Gamma12 / state.temp;
    }

    // Helmholtz free energy

    amrex::Real f1, f2, f12;
    amrex::Real f1dT, f2dT, f12dT;

    chabrier1998_helmholtz_F<do_T_derivatives>(Gamma1, Gamma1dT, f1, f1dT);
    chabrier1998_helmholtz_F<do_T_derivatives>(Gamma2, Gamma2dT, f2, f2dT);
    chabrier1998_helmholtz_F<do_T_derivatives>(Gamma12, Gamma12dT, f12, f12dT);

    // Now we add quantum correction terms discussed in Alastuey 1978.
    // Notice in Alastuey 1978, they have a different classical term,
    // which is implemented in the strong screening limit of our screen5 routine.
    amrex::Real quantum_corr_1 = 0.0_rt;

    amrex::Real quantum_corr_2 = 0.0_rt;

    [[maybe_unused]] amrex::Real quantum_corr_1_dT = 0.0_rt;
    [[maybe_unused]] amrex::Real quantum_corr_2_dT = 0.0_rt;

    if (screening_rp::enable_chabrier1998_quantum_corr) {
        // See Wallace1982, Eq. A13

        amrex::Real Gamma_eff = std::cbrt(2.0_rt) * scn_fac.z1 * scn_fac.z2 *
                         scn_fac.zs13inv * Gamma_e;
        [[maybe_unused]] amrex::Real Gamma_eff_dT;

        if constexpr (do_T_derivatives) {
            Gamma_eff_dT = -Gamma_eff / state.temp;
        }

        // TAU/3, see Wallace1982, Eq. A2

        amrex::Real tau12 = state.taufac * scn_fac.aznut;

        [[maybe_unused]] amrex::Real tau12dT;
        if constexpr (do_T_derivatives) {
            tau12dT = state.taufacdt * scn_fac.aznut;
        }

        // see Calder 2007 Eq. A8

        amrex::Real b_fac = Gamma_eff / tau12;

        [[maybe_unused]] amrex::Real b_fac_dT;
        if constexpr (do_T_derivatives) {
            b_fac_dT = (Gamma_eff_dT - b_fac * tau12dT) / tau12;
        }

        // Quantum correction terms (same as screen5)
        //see Calder 2007 Eq.A8 and Alastuey1978, Eq. 24 and 31

        quantum_corr_1 = -tau12 * (5.0_rt/32.0_rt * amrex::Math::powi<3>(b_fac) -
                                   0.014_rt * amrex::Math::powi<4>(b_fac) -
                                   0.128_rt * amrex::Math::powi<5>(b_fac));

        quantum_corr_2 = -Gamma_eff * (0.0055_rt * amrex::Math::powi<4>(b_fac) -
                                       0.0098_rt * amrex::Math::powi<5>(b_fac) +
                                       0.0048_rt * amrex::Math::powi<6>(b_fac));

        if constexpr (do_T_derivatives) {
            quantum_corr_1_dT = tau12dT / tau12 * quantum_corr_1 - tau12 *
                b_fac_dT * (15.0_rt/32.0_rt * amrex::Math::powi<2>(b_fac) -
                            0.014_rt * 4.0_rt * amrex::Math::powi<3>(b_fac) -
                            0.128_rt * 5.0_rt * amrex::Math::powi<4>(b_fac));

            quantum_corr_2_dT = Gamma_eff_dT / Gamma_eff * quantum_corr_2 - Gamma_eff *
                b_fac_dT * (0.0055_rt * 4.0_rt * amrex::Math::powi<3>(b_fac) -
                            0.0098_rt * 5.0_rt * amrex::Math::powi<4>(b_fac) +
                            0.0048_rt * 6.0_rt * amrex::Math::powi<5>(b_fac));
        }
    }
    // See Calder2007 Appendix Eq. A8.
    // f1 + f2 - f12 gives the classical terms
    // The difference between this and strong screening of screen5
    // is that we replaced the classical term which is f1 + f2 - f12
    // using results from Chabrier&Potekhin1998.

    amrex::Real h12 = f1 + f2 - f12 + quantum_corr_1 + quantum_corr_2;

    [[maybe_unused]] amrex::Real dh12dT;
    if constexpr (do_T_derivatives) {
        dh12dT = f1dT + f2dT - f12dT + quantum_corr_1_dT + quantum_corr_2_dT;
    }

    amrex::Real h12_max = 300.0_rt;
    h12 = amrex::min(h12_max, h12);

    scor = std::exp(h12);

    if constexpr (do_T_derivatives) {
        if (h12 == h12_max) {
            scordt = 0.0_rt;
        } else {
            scordt = scor * dh12dT;
        }
    }
}
#endif

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_screen(const plasma_state_t& state,
                   const scrn::screen_factors_t& scn_fac,
                   amrex::Real& scor, amrex::Real& scordt)
{
#if SCREEN_METHOD == SCREEN_METHOD_null
    // null screening
    amrex::ignore_unused(state, scn_fac);
    scor = 1.0_rt;
    scordt = 0.0_rt;
#elif SCREEN_METHOD == SCREEN_METHOD_screen5
    actual_screen5<do_T_derivatives>(state, scn_fac, scor, scordt);
#elif SCREEN_METHOD == SCREEN_METHOD_chugunov2007
    chugunov2007<do_T_derivatives>(state, scn_fac, scor, scordt);
#elif SCREEN_METHOD == SCREEN_METHOD_chugunov2009
    chugunov2009<do_T_derivatives>(state, scn_fac, scor, scordt);
#elif SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    chabrier1998<do_T_derivatives>(state, scn_fac, scor, scordt);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_screen(const plasma_state_t& state,
                   const scrn::screen_factors_t& scn_fac,
                   amrex::Real& scor, amrex::Real& scordt)
{
    if (integrator_rp::jacobian == 1) {
        constexpr int do_T_derivatives = 1;
        actual_screen<do_T_derivatives>(state, scn_fac, scor, scordt);
    } else {
        constexpr int do_T_derivatives = 0;
        actual_screen<do_T_derivatives>(state, scn_fac, scor, scordt);
        scordt = 0.0_rt;
    }
}

#endif
