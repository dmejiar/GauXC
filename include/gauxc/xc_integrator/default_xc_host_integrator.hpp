#pragma once
#include "xc_integrator_impl.hpp"

namespace GauXC  {
namespace detail {


template <typename F>
struct XCHostData {

  std::vector<F> eps;
  std::vector<F> gamma;
  std::vector<F> vrho;
  std::vector<F> vgamma;
 
  std::vector<F> zmat;
  std::vector<F> nbe_scr;
  std::vector<F> den_scr;
  std::vector<F> basis_eval;
   

  XCHostData( size_t n_deriv, 
              size_t nbf,
              size_t max_npts, 
              size_t max_npts_x_nbe ) :
    eps( max_npts ),
    gamma( (n_deriv > 0) * max_npts ),
    vrho( max_npts ),
    vgamma( (n_deriv > 0) * max_npts ),
    zmat( max_npts_x_nbe ),
    nbe_scr( nbf * nbf ),
    den_scr( (3*n_deriv + 1) * max_npts ),
    basis_eval( (3*n_deriv + 1) * max_npts_x_nbe ) { }
   

};

template <typename MatrixType>
class DefaultXCHostIntegrator : public XCIntegratorImpl<MatrixType> {

  using base_type     = XCIntegratorImpl<MatrixType>;
  using matrix_type   = typename base_type::matrix_type;
  using value_type    = typename base_type::value_type;
  using basisset_type = typename base_type::basisset_type;
  using exc_vxc_type  = typename base_type::exc_vxc_type;
    
  std::shared_ptr< XCHostData< value_type > > host_data_;

  exc_vxc_type eval_exc_vxc_( const MatrixType& ) override; 

public:

  template <typename... Args>
  DefaultXCHostIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  DefaultXCHostIntegrator( const DefaultXCHostIntegrator& ) = default;
  DefaultXCHostIntegrator( DefaultXCHostIntegrator&& ) noexcept = default;

  ~DefaultXCHostIntegrator() noexcept = default;

};




template <typename F, size_t n_deriv>
void process_batches_host_replicated_p(
  XCWeightAlg            weight_alg,
  const functional_type& func,
  const BasisSet<F>&     basis,
  const Molecule   &     mol,
  const MolMeta    &     meta,
  XCHostData<F>    &     host_data,
  std::vector< XCTask >& local_work,
  const F*               P,
  F*                     VXC,
  F*                     exc,
  F*                     n_el
);



template <typename F, typename... Args>
void process_batches_host_replicated_p( size_t n_deriv, Args&&... args ) {
  if( n_deriv == 0 )
    process_batches_host_replicated_p<F,0>( std::forward<Args>(args)... );
  else if( n_deriv == 1 )
    process_batches_host_replicated_p<F,1>( std::forward<Args>(args)... );
  else
    throw std::runtime_error("MGGA NYI");
}




template <typename MatrixType>
typename DefaultXCHostIntegrator<MatrixType>::exc_vxc_type 
  DefaultXCHostIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& P ) {

  size_t nbf = this->basis_->nbf();

  //// TODO: Check that P is sane


  auto& tasks = this->load_balancer_->get_tasks();

  size_t max_npts       = this->load_balancer_->max_npts();
  size_t max_nbe        = this->load_balancer_->max_nbe();
  size_t max_npts_x_nbe = this->load_balancer_->max_npts_x_nbe();

  size_t n_deriv = this->func_->is_gga() ? 1 : 0;

  host_data_ = std::make_shared<XCHostData<value_type>>( 
    n_deriv, nbf, max_npts, max_npts_x_nbe 
  );



  matrix_type VXC( nbf, nbf );
  value_type  EXC, N_EL;

  process_batches_host_replicated_p< value_type>(
    n_deriv,XCWeightAlg::SSF, *this->func_, *this->basis_,
    this->load_balancer_->molecule(), this->load_balancer_->molmeta(),
    *host_data_, tasks, P.data(), VXC.data(), &EXC, &N_EL 
  );

  return exc_vxc_type{EXC, std::move(VXC)};

} 

}
}