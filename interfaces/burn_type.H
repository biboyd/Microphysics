#ifndef BURN_TYPE_H
#define BURN_TYPE_H

#include <iomanip>

#include <AMReX_REAL.H>

#include <network.H>
#include <eos_type.H>
#include <eos_composition.H>
#include <extern_parameters.H>

#include <ArrayUtilities.H>

using namespace amrex::literals;
using namespace network_rp;

// A generic structure holding data necessary to do a nuclear burn.

// Set the number of independent variables -- this should be
// enuc + the number of species which participate
// in the evolution equations.

// For Strang evolution, this will be the number of equations we are
// evolving.  For simplified-SDC, we will need neqs when getting the
// reaction sources from the network's RHS.
const int neqs = 1 + NumSpec;

// Index of the energy variable in the work arrays.

const int net_ienuc = NumSpec + 1;

// this is the data type that is used to get the ydots from the actual
// RHS of the network, regardless of Strang or SDC
using YdotNetArray1D = amrex::Array1D<amrex::Real, 1, neqs>;

// these indices represent the order that the conserved state comes
// into the ODE integration from the hydro code.
//
// they also represent the order of the advective sources
//
// integrate rho*X, internal energy, total energy
// carry momentum as an unevolved variable

const int SFS   = 0;
const int SEINT = NumSpec;
// the following are not evolved
const int SFX   = SEINT + 1;
const int SRHO  = SFX + NumAux; // this is SFS + NumSpec if NumAux = 0;
const int SMX   = SRHO + 1;
const int SMY   = SRHO + 2;
const int SMZ   = SRHO + 3;
const int SEDEN = SMZ+1;

const int SVAR  = SEDEN+1;
const int SVAR_EVOLVE = SFX;

// this is the data type of the dense Jacobian that the network wants.
// it is not the same size as the Jacobian that VODE cares about when
// we are doing simplified-SDC
using  JacNetArray2D = ArrayUtil::MathArray2D<1, neqs, 1, neqs>;

struct burn_t
{

  // on exit, this will be set to the integration time we achieved

  amrex::Real time{};

  // this first group are the quantities the network RHS uses

  amrex::Real rho{};
  amrex::Real T{};
  amrex::Real e{};
  amrex::Real xn[NumSpec]{};
#if NAUX_NET > 0
  amrex::Real aux[NumAux]{};
#endif

  // scaling used to make e we integrate dimensionless
  amrex::Real e_scale{};

  // derivatives needed for calling the EOS
  amrex::Real dedr{};
  amrex::Real dedT{};
  amrex::Real dedA{};
  amrex::Real dedZ{};

  // for diagnostics / error reporting

  int i{};
  int j{};
  int k{};
#ifdef NONAKA_PLOT
  int level{};
  amrex::Real reference_time{};
#ifdef STRANG
  int strang_half{};
#endif
#endif

  // now come the bits that we need for SDC or Strang evolution

  // y is the input conserved state.  We will keep this state updated
  // in time as we integrate, such that upon output it will be the
  // final conserved state.
  amrex::Real y[SVAR]{};

  // we need to store a copy of the original state as well so we can
  // handle the non-evolved state evolution
  amrex::Real rho_orig{};

  amrex::Real rhoe_orig{};

  // ydot_a are the advective terms that will modify the state y due
  // to hydrodynamics over the timestep.
  amrex::Real ydot_a[SVAR]{};

  int sdc_iter{};
  int num_sdc_iters{};

  // for drive_initial_convection, we will fix T during the
  // integration to a passed in value.  We will interpret a positive
  // T_fixed as setting this feature.
  amrex::Real T_fixed{-1.0};

  // all coupling types need the specific heats to transform the
  // reaction Jacobian elements from T to e
  amrex::Real cv{};

  // dx is useful for estimating timescales for equilibriation
  amrex::Real dx{};

  amrex::Real mu{};
  amrex::Real mu_e{};
  amrex::Real y_e{};
  amrex::Real eta{};
  amrex::Real abar{};
  amrex::Real zbar{};

#ifdef NSE_NET
  amrex::Real mu_p{};
  amrex::Real mu_n{};
#endif

#ifdef NSE
  bool nse{};
#endif

  // diagnostics
  int n_rhs{}, n_jac{}, n_step{};

  // Was the burn successful?
  bool success{};

  // integrator error code
  short error_code{};

};


inline
std::ostream& operator<< (std::ostream& o, burn_t const& burn_state)
{
  o << std::setprecision(12);
  o << "rho = " << burn_state.rho << std::endl;
  o << "T =   " << burn_state.T << std::endl;
  o << "xn = ";
  for (double X : burn_state.xn) {
    o << X << " ";
  }
  o << std::endl;
  o << "(i, j, k) = " << burn_state.i << " " << burn_state.j << " " << burn_state.k << std::endl;
#if NAUX_NET > 0
  o << "aux = ";
  for (double A : burn_state.aux) {
    o << A << " ";
  }
  o << std::endl;
#endif

  o << "y[SRHO] = " << burn_state.y[SRHO] << std::endl;
  o << "y[SEINT] = " << burn_state.y[SEINT] << std::endl;

  o << "y[SFS:] = ";
  for (int n = 0; n < NumSpec; n++) {
      o << burn_state.y[SFS+n] << " ";
  }
  o << std::endl;

#if NAUX_NET > 0
  o << "y[SFX:] = ";
  for (int n = 0; n < NumAux; n++) {
      o << burn_state.y[SFX+n] << " ";
  }
  o << std::endl;
#endif

  o << "ydot_a[SRHO] = " << burn_state.ydot_a[SRHO] << std::endl;
  o << "ydot_a[SEINT] = " << burn_state.ydot_a[SEINT] << std::endl;

  o << "ydot_a[SFS:] = ";
  for (int n = 0; n < NumSpec; n++) {
      o << burn_state.ydot_a[SFS+n] << " ";
  }
  o << std::endl;

#if NAUX_NET > 0
  o << "ydot_a[SFX:] = ";
  for (int n = 0; n < NumAux; n++) {
      o << burn_state.ydot_a[SFX+n] << " ";
  }
  o << std::endl;
#endif
  return o;
}


// Given an eos type, copy the data relevant to the burn type.

template <typename T, typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eos_to_burn (const T& eos_state, BurnT& burn_state)
{
    burn_state.rho  = eos_state.rho;
    burn_state.T    = eos_state.T;
    burn_state.e    = eos_state.e;
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = eos_state.xn[n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        burn_state.aux[n] = eos_state.aux[n];
    }
#endif

    burn_state.cv   = eos_state.cv;

    burn_state.y_e  = eos_state.y_e;
    burn_state.eta  = eos_state.eta;
    burn_state.abar = eos_state.abar;
    burn_state.zbar = eos_state.zbar;

#ifdef NSE_NET
    burn_state.mu_p = eos_state.mu_p;
    burn_state.mu_n = eos_state.mu_n;
#endif
}



// Given a burn type, copy the data relevant to the eos type.
// Note that when doing simplified SDC integration, we should
// avoid using this interface because the energy includes a
// contribution from the advection term. However this is useful
// for instantaneous RHS evaluations.

template <typename BurnT, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_eos (const BurnT& burn_state, T& eos_state)
{
    eos_state.rho  = burn_state.rho;
    eos_state.T    = burn_state.T;
    eos_state.e    = burn_state.e;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = burn_state.xn[n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = burn_state.aux[n];
    }
#endif

    eos_state.cv   = burn_state.cv;

    eos_state.y_e  = burn_state.y_e;
    eos_state.eta  = burn_state.eta;
    eos_state.abar = burn_state.abar;
    eos_state.zbar = burn_state.zbar;

#ifdef NSE_NET
    eos_state.mu_p = burn_state.mu_p;
    eos_state.mu_n = burn_state.mu_n;
#endif
}


template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void normalize_abundances_sdc_burn (BurnT& state)
{

    // Constrain the partial densities in burn_t state to sum to the
    // density.
    //
    // This is meant to be used upon exit, and we assume that
    // vode_to_burn was already called

    // we have state.y[SRHO] to constrain to

    amrex::Real rhoX_sum = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        state.y[SFS+n] = amrex::max(amrex::min(state.y[SFS+n], state.y[SRHO]), 0.0_rt);
        rhoX_sum += state.y[SFS+n];
    }
    rhoX_sum /= state.y[SRHO];

    for (int n = 0; n < NumSpec; n++) {
        state.y[SFS+n] /= rhoX_sum;
    }

#ifdef AUX_THERMO
    // make the aux data consistent with the vode_state X's

    eos_re_t eos_state;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state.y[SFS+n] / state.y[SRHO];
    }

    set_aux_comp_from_X(eos_state);

    // also store it in the burn_t state

    for (int n = 0; n < NumAux; n++) {
        state.y[SFX+n] = state.y[SRHO] * eos_state.aux[n];
    }
#endif

}

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void normalize_abundances_burn (BurnT& state)
{

    amrex::Real sum = 0.0_rt;
    for (double &x : state.xn) {
        x = amrex::max(network_rp::small_x, amrex::min(1.0_rt, x));
        sum += x;
    }
    for (double &x : state.xn) {
        x /= sum;
    }

}

#endif
