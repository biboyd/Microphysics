#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices {
        -1, -1, 5, -1, -1, 4, -1,
        -1, 0, 3, -1, -1, 4, -1,
        -1, 1, 3, -1, -1, 5, -1,
        -1, 2, 3, -1, -1, 6, -1,
        -1, 1, 9, -1, -1, 11, -1,
        -1, 0, 11, -1, -1, 12, -1,
        -1, 3, 3, -1, 1, 9, -1,
        -1, 3, 3, -1, 2, 7, -1,
        -1, 2, 4, -1, 0, 6, -1,
        -1, 0, 5, -1, 1, 4, -1,
        -1, -1, 9, -1, -1, 8, -1,
        -1, -1, 8, -1, -1, 9, 11,
        -1, -1, 10, -1, -1, 12, 14,
        -1, -1, 12, -1, -1, 10, -1
    };
}
#endif

void actual_network_init()
{

}
