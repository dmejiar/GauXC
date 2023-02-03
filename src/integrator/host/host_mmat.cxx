#include "host_zmat.hpp"
#include "blas.hpp"

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
		     F*        nbe_scr ) {

  GauXC::blas::lacpy( 'A', nbf, npts, dbasis_x, nbf, 
                      mmat_x, nbf );
  GauXC::blas::lacpy( 'A', nbf, npts, dbasis_y, nbf, 
                      mmat_y, nbf );
  GauXC::blas::lacpy( 'A', nbf, npts, dbasis_z, nbf, 
                      mmat_z, nbf );

  for( int32_t i = 0; i < npts; ++i ) {

    const int32_t ioff = i * nbf;
    auto* x_col    = mmat_x + ioff;
    auto* y_col    = mmat_y + ioff;
    auto* z_col    = mmat_z + ioff;

    const F mgga_fact = 0.25*vtau[i];
    GauXC::blas::scal( nbf, mgga_fact, x_col, 1 );
    GauXC::blas::scal( nbf, mgga_fact, y_col, 1 );
    GauXC::blas::scal( nbf, mgga_fact, z_col, 1 );
  }

  GauXC::blas::syr2k( 'L', 'N', nbf, npts, F(1.), dbasis_x, nbf,
                       nbf, mmat_x, nbf, F(1.), nbe_scr, nbf );
  GauXC::blas::syr2k( 'L', 'N', nbf, npts, F(1.), dbasis_y, nbf,
                       nbf, mmat_y, nbf, F(1.), nbe_scr, nbf );
  GauXC::blas::syr2k( 'L', 'N', nbf, npts, F(1.), dbasis_z, nbf,
                       nbf, mmat_z, nbf, F(1.), nbe_scr, nbf );
  }
}

template
void mmat_mgga_host( int32_t    npts,
                     int32_t    nbf,
                     const float*  vtau,
                     const float*  dbasis_x,
		     const float*  dbasis_y,
		     const float*  dbasis_z,
                     float*        mmat_x,
		     float*        mmat_y,
		     float*        mmat_z,
		     float*        nbe_scr );

template
void mmat_mgga_host( int32_t    npts,
                     int32_t    nbf,
                     const double*  vtau,
                     const double*  dbasis_x,
		     const double*  dbasis_y,
		     const double*  dbasis_z,
                     double*        mmat_x,
		     double*        mmat_y,
		     double*        mmat_z,
		     double*        nbe_scr );

