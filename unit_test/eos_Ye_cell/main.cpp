#include <iostream>

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <eos_Ye_cell.H>

using namespace problem_rp;

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  std::cout << "calling the EOS on a series of zone states..." << std::endl;

  //init_unit_test()
  init_extern_parameters();

  // C++ EOS initialization (must be done after init_extern_parameters)
  eos_init(small_temp, small_dens);

  // C++ Network, RHS, screening, rates initialization
  network_init();

  eos_Ye_cell_c();

  amrex::Finalize();
}
