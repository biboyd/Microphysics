#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_N13_C13_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 34, 1, 151, 1, 6> j_N13_C13_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 151> j_N13_C13_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 34> j_N13_C13_temp;

    AMREX_GPU_MANAGED table_t j_F18_O18_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_F18_O18_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_F18_O18_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_F18_O18_temp;

    AMREX_GPU_MANAGED table_t j_Mg25_Na25_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_Mg25_Na25_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_Mg25_Na25_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_Mg25_Na25_temp;

    AMREX_GPU_MANAGED table_t j_Na23_Ne23_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_Na23_Ne23_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_Na23_Ne23_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_Na23_Ne23_temp;

    AMREX_GPU_MANAGED table_t j_Ne21_F21_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_Ne21_F21_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_Ne21_F21_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_Ne21_F21_temp;


}


void init_tabular()
{

    amrex::Print() << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_N13_C13_meta.ntemp = 34;
    j_N13_C13_meta.nrhoy = 151;
    j_N13_C13_meta.nvars = 6;
    j_N13_C13_meta.nheader = 6;

    init_tab_info(j_N13_C13_meta, "13n-13c_electroncapture.dat", j_N13_C13_rhoy, j_N13_C13_temp, j_N13_C13_data);


    j_F18_O18_meta.ntemp = 39;
    j_F18_O18_meta.nrhoy = 152;
    j_F18_O18_meta.nvars = 6;
    j_F18_O18_meta.nheader = 6;

    init_tab_info(j_F18_O18_meta, "18f-18o_electroncapture.dat", j_F18_O18_rhoy, j_F18_O18_temp, j_F18_O18_data);


    j_Mg25_Na25_meta.ntemp = 39;
    j_Mg25_Na25_meta.nrhoy = 152;
    j_Mg25_Na25_meta.nvars = 6;
    j_Mg25_Na25_meta.nheader = 7;

    init_tab_info(j_Mg25_Na25_meta, "25mg-25na_electroncapture.dat", j_Mg25_Na25_rhoy, j_Mg25_Na25_temp, j_Mg25_Na25_data);


    j_Na23_Ne23_meta.ntemp = 39;
    j_Na23_Ne23_meta.nrhoy = 152;
    j_Na23_Ne23_meta.nvars = 6;
    j_Na23_Ne23_meta.nheader = 7;

    init_tab_info(j_Na23_Ne23_meta, "23na-23ne_electroncapture.dat", j_Na23_Ne23_rhoy, j_Na23_Ne23_temp, j_Na23_Ne23_data);


    j_Ne21_F21_meta.ntemp = 39;
    j_Ne21_F21_meta.nrhoy = 152;
    j_Ne21_F21_meta.nvars = 6;
    j_Ne21_F21_meta.nheader = 6;

    init_tab_info(j_Ne21_F21_meta, "21ne-21f_electroncapture.dat", j_Ne21_F21_rhoy, j_Ne21_F21_temp, j_Ne21_F21_data);



}
