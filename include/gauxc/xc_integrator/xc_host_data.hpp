#pragma once
#include <vector>
#include <cstdint>

#include <gauxc/gauxc_config.hpp>

#ifdef GAUXC_ENABLE_HOST
namespace GauXC {

template <typename F>
struct XCHostData {

  std::vector<F> eps;
  std::vector<F> lapl;
  std::vector<F> tau;
  std::vector<F> gamma;
  std::vector<F> vrho;
  std::vector<F> vgamma;
  std::vector<F> vlapl;
  std::vector<F> vtau;
 
  std::vector<F> mmat;
  std::vector<F> zmat;
  std::vector<F> nbe_scr;
  std::vector<F> den_scr;
  std::vector<F> basis_eval;
   

  XCHostData( size_t n_deriv, 
              size_t nbf,
              size_t max_npts, 
              size_t max_npts_x_nbe ) :
    eps( max_npts ),
    lapl( (n_deriv > 0) * max_npts ),
    tau( (n_deriv > 0) * max_npts),
    gamma( (n_deriv > 0) * max_npts ),
    vrho( max_npts ),
    vgamma( (n_deriv > 0) * max_npts ),
    vlapl( (n_deriv > 0) * max_npts ),
    vtau( (n_deriv > 0) * max_npts ),
    mmat( 3*n_deriv * max_npts_x_nbe ),
    zmat( max_npts_x_nbe ),
    nbe_scr( nbf * nbf ),
    den_scr( (3*n_deriv + 1) * max_npts ),
    basis_eval( (3*n_deriv + 1) * max_npts_x_nbe ) { }
   

};

}
#endif
