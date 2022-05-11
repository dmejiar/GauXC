
#include "incore_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include "host/reference_local_host_work_driver.hpp"
#include <stdexcept>
#include "device/xc_device_aos_data.hpp"
#include <fstream>
#include <gauxc/util/unused.hpp>

#include "integrator_util/exx_screening.hpp"
#include "integrator_util/integral_bounds.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exx_( int64_t m, int64_t n, const value_type* P,
             int64_t ldp, value_type* K, int64_t ldk, 
             const IntegratorSettingsEXX& /*settings*/ ) { 


  const auto& basis = this->load_balancer_->basis();

  // Check that P / K are sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION("P/K Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION("P/K Must Have Same Dimension as Basis");
  if( ldp < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
  if( ldk < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDK");

  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();

  // Allocate Device memory
  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  auto device_data_ptr = 
    this->timer_.time_op("XCIntegrator.DeviceAlloc",
      [=](){ return lwd->create_device_data(); });

  // Compute local contributions to K and retrieve
  // data from device 
  this->timer_.time_op("XCIntegrator.LocalWork_EXX", [&](){
    exx_local_work_( basis, P, ldp, K, ldk, 
      tasks.begin(), tasks.end(), *device_data_ptr);
  });

  this->timer_.time_op("XCIntegrator.ImbalanceWait",[&](){
    MPI_Barrier(this->load_balancer_->comm());
  });  

  // Reduce Results in host mem
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){
    this->reduction_driver_->allreduce_inplace( K, nbf*nbf, ReductionOp::Sum );
  });
}



template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exx_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                       value_type* K, int64_t ldk,
                       host_task_iterator task_begin, host_task_iterator task_end,
                       XCDeviceData& device_data ) {

  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );

  // Setup Aliases
  const auto& mol   = this->load_balancer_->molecule();
  const auto& meta  = this->load_balancer_->molmeta();

  const auto nbf     = basis.nbf();
  const auto nshells = basis.nshells();


  // Get basis map and shell pairs
  BasisSetMap basis_map(basis,mol);
  ShellPairCollection shell_pairs(basis);




  // Sort tasks 
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
  };
  std::sort( task_begin, task_end, task_comparator );



  // TODO: Refactor this into separate function
  auto& lb_state = this->load_balancer_->state();

  // Modify weights if need be
  if( not lb_state.modified_weights_are_stored ) {

  integrator_term_tracker enabled_terms;
  enabled_terms.weights = true;

  this->timer_.time_op("XCIntegrator.Weights", [&]() { 
    const auto natoms = mol.natoms();
    device_data.reset_allocations();
    device_data.allocate_static_data_weights( natoms );
    device_data.send_static_data_weights( mol, meta );

    // Processes batches in groups that saturadate available device memory
    auto task_it = task_begin;
    while( task_it != task_end ) {
      
      // Determine next task batch, send relevant data to device (weights only)
      auto task_batch_end = 
        device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

      // Apply partition weights 
      lwd->partition_weights( &device_data );
      
      // Copy back to host data
      device_data.copy_weights_to_tasks( task_it, task_batch_end );

      // Update iterator
      task_it = task_batch_end;

    } // End loop over batches

    // Signal that we don't need to do weights again
    lb_state.modified_weights_are_stored = true;
  });

  }

#if 0
  // TODO: This turns off EK screening 
  std::vector<int32_t> full_shell_list(nshells);
  std::iota(full_shell_list.begin(),full_shell_list.end(),0);
  for( auto it = task_begin; it != task_end; ++it ) {
    it->cou_screening.shell_list = full_shell_list;
    it->cou_screening.nbe        = nbf;
  }
#else
  // Compute EXX screening modifications
  this->timer_.time_op("XCIntegrator.EXX_Screening", [&]() { 
    // Compute base screening quantities
    const size_t nb2 = basis.nbf() * basis.nbf();
    std::vector<double> P_abs(nb2);
    for( auto i = 0; i < nb2; ++i ) P_abs[i] = std::abs(P[i]);

    const size_t ns2 = nshells * nshells;
    std::vector<double> V_max(ns2);
    for( auto i = 0; i < nshells; ++i )
    for( auto j = 0; j < nshells; ++j ) {
      V_max[i + j*nshells] = 
        util::max_coulomb( basis.at(i), basis.at(j) );
    }
    // Create LocalHostWorkDriver
    LocalHostWorkDriver host_lwd(
      std::make_unique<ReferenceLocalHostWorkDriver>()
    );
    exx_ek_screening( basis, basis_map, P_abs.data(), basis.nbf(),
      V_max.data(), nshells, 1e-10, 1e-10, &host_lwd, task_begin,
      task_end );

    // Remove tasks with no coulomb shells
    task_end = std::stable_partition( task_begin, task_end,
      []( const auto& t ) { return t.cou_screening.shell_list.size() > 0; } );
  });
#endif

  // Populate submat maps
  device_data.populate_submat_maps( basis.nbf(), task_begin, task_end, basis_map );



  // Do EXX integration in task batches
  device_data.reset_allocations();
  device_data.allocate_static_data_exx( nbf, nshells );
  device_data.send_static_data_density_basis( P, ldp, basis );
  device_data.send_static_data_shell_pairs( basis, shell_pairs );

  // Zero integrands
  device_data.zero_exx_integrands();

  // Processes batches in groups that saturadate available device memory
  integrator_term_tracker enabled_terms;
  enabled_terms.exx = true;

  //GAUXC_GENERIC_EXCEPTION("DIE DIE DIE");
  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Determine next task batch, send relevant data to device (EXX only)
    task_it = 
      device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

    /*** Process the batches ***/

    // Evaluate collocation
    lwd->eval_collocation( &device_data );

    // Evaluate F(mu,i) = P(mu,nu) * B(nu,i)
    // mu runs over significant ek shells
    // nu runs over the bfn shell list
    // i runs over all points
    lwd->eval_exx_fmat( &device_data );

    // Compute G(mu,i) = w(i) * A(mu,nu,i) * F(nu,i)
    // mu/nu run over significant ek shells
    // i runs over all points
    lwd->eval_exx_gmat( &device_data, basis_map );

    // Increment K(mu,nu) += B(mu,i) * G(nu,i)
    // mu runs over bfn shell list
    // nu runs over ek shells
    // i runs over all points
    lwd->inc_exx_k( &device_data );

  } // Loop over batches of batches 

#if 0
  // Symmetrize K in device memory
  lwd->symmetrize_exx_k( &device_data );
#endif

  device_data.retrieve_exx_integrands( K, ldk );
}

}
}
