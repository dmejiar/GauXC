#pragma once
#include <cstdint>

namespace GauXC  {
namespace integrator::host {

template <typename F>
void mmat_mgga_host( int32_t   npts,
                     int32_t   nbf,
                     const F*  vtau,
                     const F*  dbasis_x,
                     const F*  dbasis_y,
                     const F*  dbasis_z,
                     F*        mmat_x,
		     F*        mmat_y,
		     F*        mmat_z,
		     F*        nbe_scr); 

}
}
