#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_Na23_Ne23_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_Na23_Ne23_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_Na23_Ne23_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_Na23_Ne23_temp;

    AMREX_GPU_MANAGED table_t j_Ne23_Na23_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_Ne23_Na23_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_Ne23_Na23_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_Ne23_Na23_temp;

    AMREX_GPU_MANAGED table_t j_Na25_Mg25_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_Na25_Mg25_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_Na25_Mg25_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_Na25_Mg25_temp;

    AMREX_GPU_MANAGED table_t j_Mg25_Na25_meta;
    AMREX_GPU_MANAGED amrex::Array3D<amrex::Real, 1, 39, 1, 152, 1, 6> j_Mg25_Na25_data;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 152> j_Mg25_Na25_rhoy;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, 39> j_Mg25_Na25_temp;


}


void init_tabular()
{

    amrex::Print() << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_Na23_Ne23_meta.ntemp = 39;
    j_Na23_Ne23_meta.nrhoy = 152;
    j_Na23_Ne23_meta.nvars = 6;
    j_Na23_Ne23_meta.nheader = 7;

    init_tab_info(j_Na23_Ne23_meta, "23na-23ne_electroncapture.dat", j_Na23_Ne23_rhoy, j_Na23_Ne23_temp, j_Na23_Ne23_data);


    j_Ne23_Na23_meta.ntemp = 39;
    j_Ne23_Na23_meta.nrhoy = 152;
    j_Ne23_Na23_meta.nvars = 6;
    j_Ne23_Na23_meta.nheader = 5;

    init_tab_info(j_Ne23_Na23_meta, "23ne-23na_betadecay.dat", j_Ne23_Na23_rhoy, j_Ne23_Na23_temp, j_Ne23_Na23_data);


    j_Na25_Mg25_meta.ntemp = 39;
    j_Na25_Mg25_meta.nrhoy = 152;
    j_Na25_Mg25_meta.nvars = 6;
    j_Na25_Mg25_meta.nheader = 6;

    init_tab_info(j_Na25_Mg25_meta, "25na-25mg_betadecay.dat", j_Na25_Mg25_rhoy, j_Na25_Mg25_temp, j_Na25_Mg25_data);


    j_Mg25_Na25_meta.ntemp = 39;
    j_Mg25_Na25_meta.nrhoy = 152;
    j_Mg25_Na25_meta.nvars = 6;
    j_Mg25_Na25_meta.nheader = 7;

    init_tab_info(j_Mg25_Na25_meta, "25mg-25na_electroncapture.dat", j_Mg25_Na25_rhoy, j_Mg25_Na25_temp, j_Mg25_Na25_data);



}
