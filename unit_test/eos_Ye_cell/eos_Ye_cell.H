#ifndef EOS_YE_CELL_H
#define EOS_YE_CELL_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <fstream>
#include <react_util.H>

using namespace problem_rp;

AMREX_INLINE
void write_out(std::ofstream& out, eos_t& curr_state){
    int ic12 = network_spec_index("carbon-12");
    int io16 = network_spec_index("oxygen-16");
    int ine23 = network_spec_index("neon-23");
    int ina23 = network_spec_index("sodium-23");
    //Header rho,T,p,C12,O16,Ne23,Na23,Ye,S,e,eta,

    out << curr_state.rho       << ',' 
        << curr_state.T         << ','
        << curr_state.p         << ','
        << curr_state.xn[ic12]  << ','
        << curr_state.xn[io16]  << ','
        << curr_state.xn[ine23] << ','
        << curr_state.xn[ina23] << ','
        << curr_state.y_e       << ','
        << curr_state.s         << ','
        << curr_state.e         << ','
        << curr_state.eta       << ','
        << std::endl;
}

AMREX_INLINE
void eos_Ye_cell_c()
{

    eos_t state;

    // check on 2 of rho,T,P are set.
    if (density > 0. && temperature > 0. && pressure > 0.){
        amrex::Error("All of density, temperature and pressure are set. Make one less that zero to run so EOS solves for missing component.");
    }


    // Set mass fractions to sanitize inputs for them
    amrex::Real massfractions[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        massfractions[n] = 0.;
    }

    //initialize mass fractions
    int ic12 = network_spec_index("carbon-12");
    int io16 = network_spec_index("oxygen-16");
    int ine23 = network_spec_index("neon-23");
    int ina23 = network_spec_index("sodium-23");

    massfractions[ic12] = XC12;
    massfractions[io16] = XO16;
    massfractions[ine23] = 0.5_rt * XUrca23;
    massfractions[ina23] = 0.5_rt * XUrca23;

    // normalize mass fractions check 
    // Make sure user set all the mass fractions to values in the interval [0, 1]
    amrex::Real sumX = XC12 + XO16 + XUrca23;
    for (int n = 0; n < NumSpec; ++n) {
        if (massfractions[n-1] < 0 || massfractions[n-1] > 1) {
            amrex::Error("mass fraction for " + short_spec_names_cxx[n-1] + " not initialized in the interval [0,1]!");
        }
        massfractions[n] /= sumX;
    }

    // open outfile and write header
    std::ofstream out_stream(outfile);
    out_stream << ",rho,T,p,C12,O16,Ne23,Na23,Ye,S,e,eta," << std::endl;

    for (int i = 0; i < nUrca; i++){

        // set Urca mass frac
        amrex::Real frac = static_cast<amrex::Real>(nUrca-1 - i)/static_cast<amrex::Real>(nUrca-1);

        massfractions[ine23] = frac * XUrca23;
        massfractions[ina23] =  (1. - frac) * XUrca23;

        // fill state with composition
        for (int n = 0; n < NumSpec; ++n) {
            state.xn[n] = massfractions[n];
        }

        if (pressure < 0.){
            state.T   = temperature;
            state.rho = density;
            eos(eos_input_rt, state);
        }
        else if (density < 0.){
            state.T = temperature;
            state.p = pressure;
            eos(eos_input_tp, state);
        } else if (temperature < 0.){
            state.rho = density;
            state.p   = pressure;
            eos(eos_input_rp, state);
        }

    // write out curr state
    out_stream << i << ',';
    write_out(out_stream, state);
    }

}
#endif
