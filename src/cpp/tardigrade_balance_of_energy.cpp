/**
  ******************************************************************************
  * \file tardigrade_balance_of_energy.cpp
  ******************************************************************************
  * The source file for the equations associated with the balance of energy
  ******************************************************************************
  */

#include "tardigrade_balance_of_energy.h"
#include "tardigrade_balance_of_mass.h"
#include<numeric>
#include<algorithm>
#include<functional>

namespace tardigradeBalanceEquations{

    namespace balanceOfEnergy{

        template<class floatVector_iter, class secondOrderTensor_iter>
        void computeBalanceOfEnergyNonDivergence( const floatType &density, const floatType &density_dot,
                                                  const floatVector_iter &density_gradient_begin, const floatVector_iter &density_gradient_end,
                                                  const floatType &internal_energy, const floatType &internal_energy_dot,
                                                  const floatVector_iter &internal_energy_gradient_begin, const floatVector_iter &internal_energy_gradient_end,
                                                  const floatVector_iter &velocity_begin, const floatVector_iter &velocity_end,
                                                  const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                                  const secondOrderTensor_iter &cauchy_stress_begin, const secondOrderTensor_iter &cauchy_stress_end,
                                                  const floatType &volume_fraction,
                                                  const floatType &internal_heat_generation,
                                                  const floatVector_iter &net_interphase_force_begin, const floatVector_iter &net_interphase_force_end,
                                                  floatType &result ){
            /*!
             * Compute the non-divergence parts of the balance of energy i.e.
             * 
             * \f$ \frac{\partial}{\partial t}\left( \rho^{\alpha} e^{\alpha} \right) + \left( \rho^{\alpha} v_i^{\alpha} e^{\alpha} \right)_{,i} - \frac{1}{2} c^{\alpha} v_i^{\alpha} v_i^{\alpha} + \sum_{\beta} \pi_i^{\alpha \beta} v_i^{\alpha} - \phi^{\alpha}\sigma_{ji}^{\alpha}v_{i,j}^{\alpha} - \rho^{\alpha} r^{\alpha} \f$
             * 
             * \param &density: The apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot: The partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy: The internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot: The partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction: The volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation: The internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &result: The result of the non-divergence part of the balance of energy
             */
            
            // Compute the mass change rate 
            floatType mass_change_rate;

            tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, density_gradient_begin, density_gradient_end, velocity_begin, velocity_end,
                                                                             velocity_gradient_begin, velocity_gradient_end, mass_change_rate );

            // Compute the trace of the velocity gradient
            floatType trace_velocity_gradient = 0.;
            for ( unsigned int i = 0; i < dim; i++ ){
                trace_velocity_gradient += *( velocity_gradient_begin + dim * i + i );
            }

            result = density_dot * internal_energy + density * internal_energy_dot
                   + internal_energy * std::inner_product( density_gradient_begin, density_gradient_end, velocity_begin, 0. )
                   + density * internal_energy * trace_velocity_gradient
                   + density * std::inner_product( velocity_begin, velocity_end, internal_energy_gradient_begin, 0. )
                   - 0.5 * mass_change_rate * std::inner_product( velocity_begin, velocity_end, velocity_begin, 0. )
                   + std::inner_product( net_interphase_force_begin, net_interphase_force_end, velocity_begin, 0. )
                   - density * internal_heat_generation;

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    result -= volume_fraction * ( *( cauchy_stress_begin + dim * j + i ) ) * ( *( velocity_gradient_begin + dim * i + j ) );

                }

            }

        }

        template<class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out>
        void computeBalanceOfEnergyNonDivergence( const floatType &density, const floatType &density_dot,
                                                  const floatVector_iter &density_gradient_begin, const floatVector_iter &density_gradient_end,
                                                  const floatType &internal_energy, const floatType &internal_energy_dot,
                                                  const floatVector_iter &internal_energy_gradient_begin, const floatVector_iter &internal_energy_gradient_end,
                                                  const floatVector_iter &velocity_begin, const floatVector_iter &velocity_end,
                                                  const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                                  const secondOrderTensor_iter &cauchy_stress_begin, const secondOrderTensor_iter &cauchy_stress_end,
                                                  const floatType &volume_fraction,
                                                  const floatType &internal_heat_generation,
                                                  const floatVector_iter &net_interphase_force_begin, const floatVector_iter &net_interphase_force_end,
                                                  floatType &result,
                                                  floatType &dRdRho, floatType &dRdRhoDot, floatVector_iter_out dRdGradRho_begin, floatVector_iter_out dRdGradRho_end,
                                                  floatType &dRdE, floatType &dRdEDot, floatVector_iter_out dRdGradE_begin, floatVector_iter_out dRdGradE_end,
                                                  floatVector_iter_out dRdV_begin, floatVector_iter_out dRdV_end,
                                                  secondOrderTensor_iter_out dRdGradV_begin, secondOrderTensor_iter_out dRdGradV_end,
                                                  secondOrderTensor_iter_out dRdCauchy_begin, secondOrderTensor_iter_out dRdCauchy_end,
                                                  floatType &dRdPhi, floatType &dRdr,
                                                  floatVector_iter_out dRdpi_begin, floatVector_iter_out dRdpi_end ){
            /*!
             * Compute the non-divergence parts of the balance of energy including the Jacobians i.e.
             * 
             * \f$ \frac{\partial}{\partial t}\left( \rho^{\alpha} e^{\alpha} \right) + \left( \rho^{\alpha} v_i^{\alpha} e^{\alpha} \right)_{,i} - \frac{1}{2} c^{\alpha} v_i^{\alpha} v_i^{\alpha} + \sum_{\beta} \pi_i^{\alpha \beta} v_i^{\alpha} - \phi^{\alpha}\sigma_{ji}^{\alpha}v_{i,j}^{\alpha} - \rho^{\alpha} r^{\alpha} \f$
             * 
             * \param &density: The apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot: The partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy: The internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot: The partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction: The volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation: The internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &result: The result of the non-divergence part of the balance of energy
             * \param &dRdRho: The Jacobian w.r.t. the apparent density
             * \param &dRdRhoDot: The Jacobian w.r.t. the partial temporal derivative of the apparent density
             * \param &dRdGradRho_begin: The starting iterator of the Jacobian w.r.t. the spatial gradient of the apparent density
             * \param &dRdGradRho_end: The stopping iterator of the Jacobian w.r.t. the spatial gradient of the apparent density
             * \param &dRdE: The Jacobian w.r.t. the internal energy
             * \param &dRdEDot: The Jacobian w.r.t. the partial temporal derivative of the internal energy
             * \param &dRdGradE_begin: The starting iterator of the Jacobian w.r.t. the spatial gradient of the internal energy
             * \param &dRdGradE_end: The stopping iterator of the Jacobian w.r.t. the spatial gradient of the internal energy
             * \param &dRdV_begin: The starting iterator of the Jacobian w.r.t. the velocity
             * \param &dRdV_end: The stopping iterator of the Jacobian w.r.t. the velocity
             * \param &dRdGradV_begin: The starting iterator of the Jacobian w.r.t. the velocity gradient
             * \param &dRdGradV_end: The stopping iterator of the Jacobian w.r.t. the velocity gradient
             * \param &dRdCauchy_begin: The starting iterator of the Jacobian w.r.t. the Cauchy stress
             * \param &dRdCauchy_end: The stopping iterator of the Jacobian w.r.t. the Cauchy stress
             * \param &dRdPhi: The Jacobian w.r.t. the volume fraction
             * \param &dRdr: The Jacobian w.r.t. the internal heat generation
             * \param &dRdpi_begin: The starting iterator of the Jacobian w.r.t. the net interphase force
             * \param &dRdpi_end: The stopping iterator of the Jacobian w.r.t. the net interphase force
             */
            
            // Compute the mass change rate 
            floatType mass_change_rate;

            floatType dCdRho, dCdRhoDot;

            floatVector dCdGradRho, dCdV;

            secondOrderTensor dCdGradV;

            tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass( density, density_dot, density_gradient_begin, density_gradient_end, velocity_begin, velocity_end,
                                                                             velocity_gradient_begin, velocity_gradient_end, mass_change_rate,
                                                                             dCdRho, dCdRhoDot, std::begin( dCdGradRho ), std::end( dCdGradRho ),
                                                                             std::begin( dCdV ), std::end( dCdV ), std::begin( dCdGradV ), std::end( dCdGradV ) );

            // Compute the trace of the velocity gradient
            floatType trace_velocity_gradient = 0.;
            for ( unsigned int i = 0; i < dim; i++ ){
                trace_velocity_gradient += *( velocity_gradient_begin + dim * i + i );
            }

            floatType v_dot_v = std::inner_product( velocity_begin, velocity_end, velocity_begin, 0. );
            floatType grad_rho_dot_v = std::inner_product( density_gradient_begin, density_gradient_end, velocity_begin, 0. );

            result = density_dot * internal_energy + density * internal_energy_dot
                   + internal_energy * grad_rho_dot_v
                   + density * internal_energy * trace_velocity_gradient
                   + density * std::inner_product( velocity_begin, velocity_end, internal_energy_gradient_begin, 0. )
                   - 0.5 * mass_change_rate * v_dot_v
                   + std::inner_product( net_interphase_force_begin, net_interphase_force_end, velocity_begin, 0. )
                   - density * internal_heat_generation;

            dRdRho = internal_energy_dot + internal_energy * trace_velocity_gradient + std::inner_product( velocity_begin, velocity_end, internal_energy_gradient_begin, 0. )
                   - internal_heat_generation - 0.5 * dCdRho * v_dot_v;

            dRdRhoDot = internal_energy - 0.5 * dCdRhoDot * v_dot_v;

            std::transform( velocity_begin, velocity_end, dRdGradRho_begin, std::bind( std::multiplies<floatType>(), std::placeholders::_1, internal_energy ) );

            dRdE = density_dot + std::inner_product( density_gradient_begin, density_gradient_end, velocity_begin, 0. )
                 + density * trace_velocity_gradient;

            dRdEDot = density;

            std::transform( velocity_begin, velocity_end, dRdGradE_begin, std::bind( std::multiplies<floatType>(), std::placeholders::_1, density ) );

            std::fill( dRdV_begin, dRdV_end, 0. );

            std::fill( dRdGradV_begin, dRdGradV_end, 0. );

            std::fill( dRdCauchy_begin, dRdCauchy_end, 0. );

            dRdPhi = 0;

            dRdr = -density;

            std::copy( velocity_begin, velocity_end, dRdpi_begin );

            for ( unsigned int i = 0; i < dim; i++ ){

                *( dRdGradRho_begin + i ) -= 0.5 * dCdGradRho[ i ] * v_dot_v;

                *( dRdV_begin + i ) += internal_energy * ( *( density_gradient_begin + i ) )
                                     + density * ( *( internal_energy_gradient_begin + i ) )
                                     - 0.5 * dCdV[ i ] * v_dot_v
                                     - mass_change_rate * ( *( velocity_begin + i ) )
                                     + *( net_interphase_force_begin + i );

                *( dRdGradV_begin + dim * i + i ) += density * internal_energy;

                for ( unsigned int j = 0; j < dim; j++ ){

                    dRdPhi -= ( *( cauchy_stress_begin + dim * j + i ) ) * ( *( velocity_gradient_begin + dim * i + j ) );

                    *( dRdGradV_begin + dim * i + j ) -= 0.5 * dCdGradV[ dim * i + j ] * v_dot_v + volume_fraction * ( *( cauchy_stress_begin + dim * j + i ) );

                    *( dRdCauchy_begin + dim * i + j ) -= volume_fraction * ( *( velocity_gradient_begin + dim * j + i ) );

                }

            }

            result += volume_fraction * dRdPhi;

        }

        template<class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class scalarArray_iter_out>
        void computeBalanceOfEnergyNonDivergence( const scalarArray_iter &density_begin, const scalarArray_iter &density_end,
                                                  const scalarArray_iter &density_dot_begin, const scalarArray_iter &density_dot_end,
                                                  const floatVector_iter &density_gradient_begin, const floatVector_iter &density_gradient_end,
                                                  const scalarArray_iter &internal_energy_begin, const scalarArray_iter &internal_energy_end,
                                                  const scalarArray_iter &internal_energy_dot_begin, const scalarArray_iter &internal_energy_dot_end,
                                                  const floatVector_iter &internal_energy_gradient_begin, const floatVector_iter &internal_energy_gradient_end,
                                                  const floatVector_iter &velocity_begin, const floatVector_iter &velocity_end,
                                                  const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                                  const secondOrderTensor_iter &cauchy_stress_begin, const secondOrderTensor_iter &cauchy_stress_end,
                                                  const scalarArray_iter &volume_fraction_begin, const scalarArray_iter &volume_fraction_end,
                                                  const scalarArray_iter &internal_heat_generation_begin, const scalarArray_iter &internal_heat_generation_end,
                                                  const floatVector_iter &net_interphase_force_begin, const floatVector_iter &net_interphase_force_end,
                                                  scalarArray_iter_out result_begin, scalarArray_iter_out result_end ){
            /*!
             * Compute the non-divergence parts of the balance of energy i.e.
             * 
             * \f$ \frac{\partial}{\partial t}\left( \rho^{\alpha} e^{\alpha} \right) + \left( \rho^{\alpha} v_i^{\alpha} e^{\alpha} \right)_{,i} - \frac{1}{2} c^{\alpha} v_i^{\alpha} v_i^{\alpha} + \sum_{\beta} \pi_i^{\alpha \beta} v_i^{\alpha} - \phi^{\alpha}\sigma_{ji}^{\alpha}v_{i,j}^{\alpha} - \rho^{\alpha} r^{\alpha} \f$
             * 
             * \param &density_begin: The starting iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_end: The stopping iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot_begin: The starting iterator of the partial temporal derivative of the apparent density (dm / dv) of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_dot_end: The partial temporal derivative of the apparent density (dm / dv) of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_begin: The starting iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_end: The stopping iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot_begin: The starting iterator of the partial temporal derivative of the internal energy of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_dot_end: The stopping iterator of the partial temporal derivative of the internal energy of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation_begin: The starting iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &internal_heat_generation_end: The stopping iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &result_begin: The starting iterator of the result of the non-divergence part of the balance of energy
             * \param &result_end: The stopping iterator of the result of the non-divergence part of the balance of energy
             */

            for ( auto rho = density_begin; rho != density_end; rho++ ){

                unsigned int phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfEnergyNonDivergence( *( density_begin + phase ), *( density_dot_begin + phase ), density_gradient_begin + dim * phase, density_gradient_begin + dim * ( phase + 1 ),
                                                     *( internal_energy_begin + phase ), *( internal_energy_dot_begin + phase ),
                                                     internal_energy_gradient_begin + dim * phase, internal_energy_gradient_begin + dim * ( phase + 1 ),
                                                     velocity_begin + dim * phase, velocity_begin + dim * ( phase + 1 ),
                                                     velocity_gradient_begin + sot_dim * phase, velocity_gradient_begin + sot_dim * ( phase + 1 ),
                                                     cauchy_stress_begin + sot_dim * phase, cauchy_stress_begin + sot_dim * ( phase + 1 ),
                                                     *( volume_fraction_begin + phase ), *( internal_heat_generation_begin + phase ),
                                                     net_interphase_force_begin + dim * phase, net_interphase_force_begin + dim * ( phase + 1 ),
                                                     *( result_begin + phase ) );

            }

        }

        template<class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class scalarArray_iter_out, class floatVector_iter_out, class secondOrderTensor_iter_out>
        void computeBalanceOfEnergyNonDivergence( const scalarArray_iter &density_begin,                  const scalarArray_iter &density_end,
                                                  const scalarArray_iter &density_dot_begin,              const scalarArray_iter &density_dot_end,
                                                  const floatVector_iter &density_gradient_begin,         const floatVector_iter &density_gradient_end,
                                                  const scalarArray_iter &internal_energy_begin,          const scalarArray_iter &internal_energy_end,
                                                  const scalarArray_iter &internal_energy_dot_begin,      const scalarArray_iter &internal_energy_dot_end,
                                                  const floatVector_iter &internal_energy_gradient_begin, const floatVector_iter &internal_energy_gradient_end,
                                                  const floatVector_iter &velocity_begin,                 const floatVector_iter &velocity_end,
                                                  const secondOrderTensor_iter &velocity_gradient_begin,  const secondOrderTensor_iter &velocity_gradient_end,
                                                  const secondOrderTensor_iter &cauchy_stress_begin,      const secondOrderTensor_iter &cauchy_stress_end,
                                                  const scalarArray_iter &volume_fraction_begin,          const scalarArray_iter &volume_fraction_end,
                                                  const scalarArray_iter &internal_heat_generation_begin, const scalarArray_iter &internal_heat_generation_end,
                                                  const floatVector_iter &net_interphase_force_begin,     const floatVector_iter &net_interphase_force_end,
                                                  scalarArray_iter_out result_begin,                      scalarArray_iter_out result_end,
                                                  scalarArray_iter_out dRdRho_begin,                      scalarArray_iter_out dRdRho_end,
                                                  scalarArray_iter_out dRdRhoDot_begin,                   scalarArray_iter_out dRdRhoDot_end,
                                                  floatVector_iter_out dRdGradRho_begin,                  floatVector_iter_out dRdGradRho_end,
                                                  scalarArray_iter_out dRdE_begin,                        scalarArray_iter_out dRdE_end,
                                                  scalarArray_iter_out dRdEDot_begin,                     scalarArray_iter_out dRdEDot_end,
                                                  floatVector_iter_out dRdGradE_begin,                    floatVector_iter_out dRdGradE_end,
                                                  floatVector_iter_out dRdV_begin,                        floatVector_iter_out dRdV_end,
                                                  secondOrderTensor_iter_out dRdGradV_begin,              secondOrderTensor_iter_out dRdGradV_end,
                                                  secondOrderTensor_iter_out dRdCauchy_begin,             secondOrderTensor_iter_out dRdCauchy_end,
                                                  scalarArray_iter_out dRdPhi_begin,                      scalarArray_iter_out dRdPhi_end,
                                                  scalarArray_iter_out dRdr_begin,                        scalarArray_iter_out dRdr_end,
                                                  floatVector_iter_out dRdpi_begin,                       floatVector_iter_out dRdpi_end ){
            /*!
             * Compute the non-divergence parts of the balance of energy i.e.
             * 
             * \f$ \frac{\partial}{\partial t}\left( \rho^{\alpha} e^{\alpha} \right) + \left( \rho^{\alpha} v_i^{\alpha} e^{\alpha} \right)_{,i} - \frac{1}{2} c^{\alpha} v_i^{\alpha} v_i^{\alpha} + \sum_{\beta} \pi_i^{\alpha \beta} v_i^{\alpha} - \phi^{\alpha}\sigma_{ji}^{\alpha}v_{i,j}^{\alpha} - \rho^{\alpha} r^{\alpha} \f$
             * 
             * \param &density_begin: The starting iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_end: The stopping iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot_begin: The starting iterator of the partial temporal derivative of the apparent density (dm / dv) of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_dot_end: The partial temporal derivative of the apparent density (dm / dv) of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_begin: The starting iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_end: The stopping iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot_begin: The starting iterator of the partial temporal derivative of the internal energy of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_dot_end: The stopping iterator of the partial temporal derivative of the internal energy of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation_begin: The starting iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &internal_heat_generation_end: The stopping iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &result_begin: The starting iterator of the result of the non-divergence part of the balance of energy
             * \param &result_end: The stopping iterator of the result of the non-divergence part of the balance of energy
             * \param &dRdRho_begin: The starting iterator of the Jacobian w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the Jacobian w.r.t. the apparent density
             * \param &dRdRhoDot_begin: The starting iterator of the Jacobian w.r.t. the partial temporal derivative of the apparent density
             * \param &dRdRhoDot_end: The stopping iterator of the Jacobian w.r.t. the partial temporal derivative of the apparent density
             * \param &dRdGradRho_begin: The starting iterator of the Jacobian w.r.t. the spatial gradient of the apparent density
             * \param &dRdGradRho_end: The stopping iterator of the Jacobian w.r.t. the spatial gradient of the apparent density
             * \param &dRdE_begin: The starting iterator of the Jacobian w.r.t. the internal energy
             * \param &dRdE_end: The stopping iterator of the Jacobian w.r.t. the internal energy
             * \param &dRdEDot_begin: The starting iterator of the Jacobian w.r.t. the partial temporal derivative of the internal energy
             * \param &dRdEDot_end: The stopping iterator of the Jacobian w.r.t. the partial temporal derivative of the internal energy
             * \param &dRdGradE_begin: The starting iterator of the Jacobian w.r.t. the spatial gradient of the internal energy
             * \param &dRdGradE_end: The stopping iterator of the Jacobian w.r.t. the spatial gradient of the internal energy
             * \param &dRdV_begin: The starting iterator of the Jacobian w.r.t. the velocity
             * \param &dRdV_end: The stopping iterator of the Jacobian w.r.t. the velocity
             * \param &dRdGradV_begin: The starting iterator of the Jacobian w.r.t. the velocity gradient
             * \param &dRdGradV_end: The stopping iterator of the Jacobian w.r.t. the velocity gradient
             * \param &dRdCauchy_begin: The starting iterator of the Jacobian w.r.t. the Cauchy stress
             * \param &dRdCauchy_end: The stopping iterator of the Jacobian w.r.t. the Cauchy stress
             * \param &dRdPhi_begin: The starting iterator of the Jacobian w.r.t. the volume fraction
             * \param &dRdPhi_end: The stopping iterator of the Jacobian w.r.t. the volume fraction
             * \param &dRdr_begin: The starting iterator of the Jacobian w.r.t. the internal heat generation
             * \param &dRdr_end: The stopping iterator of the Jacobian w.r.t. the internal heat generation
             * \param &dRdpi_begin: The starting iterator of the Jacobian w.r.t. the net interphase force
             * \param &dRdpi_end: The stopping iterator of the Jacobian w.r.t. the net interphase force
             */

            for ( auto rho = density_begin; rho != density_end; rho++ ){

                unsigned int phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfEnergyNonDivergence( *( density_begin + phase ), *( density_dot_begin + phase ), density_gradient_begin + dim * phase, density_gradient_begin + dim * ( phase + 1 ),
                                                     *( internal_energy_begin + phase ), *( internal_energy_dot_begin + phase ),
                                                     internal_energy_gradient_begin + dim * phase, internal_energy_gradient_begin + dim * ( phase + 1 ),
                                                     velocity_begin + dim * phase, velocity_begin + dim * ( phase + 1 ),
                                                     velocity_gradient_begin + sot_dim * phase, velocity_gradient_begin + sot_dim * ( phase + 1 ),
                                                     cauchy_stress_begin + sot_dim * phase, cauchy_stress_begin + sot_dim * ( phase + 1 ),
                                                     *( volume_fraction_begin + phase ), *( internal_heat_generation_begin + phase ),
                                                     net_interphase_force_begin + dim * phase, net_interphase_force_begin + dim * ( phase + 1 ),
                                                     *( result_begin + phase ),
                                                     *( dRdRho_begin + phase ), *( dRdRhoDot_begin + phase ), dRdGradRho_begin + dim * phase, dRdGradRho_begin + dim * ( phase + 1 ),
                                                     *( dRdE_begin + phase ), *( dRdEDot_begin + phase ), dRdGradE_begin + dim * phase, dRdGradE_begin + dim * ( phase + 1 ),
                                                     dRdV_begin + dim * phase, dRdV_begin + dim * ( phase + 1 ),
                                                     dRdGradV_begin + sot_dim * phase, dRdGradV_begin + sot_dim * ( phase + 1 ),
                                                     dRdCauchy_begin + sot_dim * phase, dRdCauchy_begin + sot_dim * ( phase + 1 ),
                                                     *( dRdPhi_begin + phase ), *( dRdr_begin + phase ),
                                                     dRdpi_begin + dim * phase, dRdpi_begin + dim * ( phase + 1 ) );

            }

        }

        template<class floatVector_iter>
        void computeBalanceOfEnergyDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                               const floatVector_iter &heat_flux_begin,              const floatVector_iter &heat_flux_end,
                                               floatType &result ){
            /*!
             * Compute the divergence part of the balance of energy for a variationally based method i.e.:
             * 
             * \f$ \psi q_{i,i} = \left( \psi q_i \right)_{,i} - \psi_{,i} q_i \f$
             *
             * where we will compute the second term as the first term is the heat flux boundary condition.
             * 
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &result: The resulting divergence term
             */

            result = -std::inner_product( test_function_gradient_begin, test_function_gradient_end, heat_flux_begin, 0. );

        }

        template<class floatVector_iter, class floatVector_iter_out>
        void computeBalanceOfEnergyDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                               const floatVector_iter &heat_flux_begin,              const floatVector_iter &heat_flux_end,
                                               floatType &result,
                                               floatVector_iter_out dRdGradPsi_begin, floatVector_iter_out dRdGradPsi_end,
                                               floatVector_iter_out dRdq_begin,       floatVector_iter_out dRdq_end ){
            /*!
             * Compute the divergence part of the balance of energy for a variationally based method i.e.:
             * 
             * \f$ \psi q_{i,i} = \left( \psi q_i \right)_{,i} - \psi_{,i} q_i \f$
             *
             * where we will compute the second term as the first term is the heat flux boundary condition.
             * 
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &result: The resulting divergence term
             * \param &dRdGradPsi_begin: The starting iterator of the Jacobian of the result w.r.t. the test function gradient
             * \param &dRdGradPsi_end: The stopping iterator of the Jacobian of the result w.r.t. the test function gradient
             * \param &dRdq_begin: The starting iterator of the Jacobian of the result w.r.t. the heat flux
             * \param &dRdq_end: The stopping iterator of the Jacobian of the result w.r.t. the heat flux
             */

            result = -std::inner_product( test_function_gradient_begin, test_function_gradient_end, heat_flux_begin, 0. );

            std::transform( heat_flux_begin, heat_flux_end, dRdGradPsi_begin, std::bind( std::multiplies<floatType>(), std::placeholders::_1, -1. ) );

            std::transform( test_function_gradient_begin, test_function_gradient_end, dRdq_begin, std::bind( std::multiplies<floatType>(), std::placeholders::_1, -1. ) );

        }

        template<class floatVector_iter, class scalarArray_iter_out>
        void computeBalanceOfEnergyDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                               const floatVector_iter &heat_flux_begin,              const floatVector_iter &heat_flux_end,
                                               scalarArray_iter_out result_begin, scalarArray_iter_out result_end ){
            /*!
             * Compute the multiphase divergence part of the balance of energy for a variationally based method i.e.:
             * 
             * \f$ \psi q_{i,i} = \left( \psi q_i \right)_{,i} - \psi_{,i} q_i \f$
             *
             * where we will compute the second term as the first term is the heat flux boundary condition.
             * 
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &result_begin: The starting iterator of the resulting divergence term
             * \param &result_end: The stopping iterator of the resulting divergence term
             */

            for ( auto r = result_begin; r != result_end; r++ ){

                const unsigned int phase = ( unsigned int )( r - result_begin );

                computeBalanceOfEnergyDivergence( test_function_gradient_begin,  test_function_gradient_end,
                                                  heat_flux_begin + dim * phase, heat_flux_begin + dim * ( phase + 1 ),
                                                  *r );

            }

        }

        template<class floatVector_iter, class scalarArray_iter_out, class floatVector_iter_out>
        void computeBalanceOfEnergyDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                               const floatVector_iter &heat_flux_begin,              const floatVector_iter &heat_flux_end,
                                               scalarArray_iter_out result_begin, scalarArray_iter_out result_end,
                                               floatVector_iter_out dRdGradPsi_begin, floatVector_iter_out dRdGradPsi_end,
                                               floatVector_iter_out dRdq_begin,       floatVector_iter_out dRdq_end ){
            /*!
             * Compute the multiphase divergence part of the balance of energy for a variationally based method i.e.:
             * 
             * \f$ \psi q_{i,i} = \left( \psi q_i \right)_{,i} - \psi_{,i} q_i \f$
             *
             * where we will compute the second term as the first term is the heat flux boundary condition.
             * 
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &result_begin: The starting iterator of the resulting divergence term
             * \param &result_end: The stopping iterator of the resulting divergence term
             * \param &dRdGradPsi_begin: The starting iterator of the Jacobian of the result w.r.t. the test function gradient
             * \param &dRdGradPsi_end: The stopping iterator of the Jacobian of the result w.r.t. the test function gradient
             * \param &dRdq_begin: The starting iterator of the Jacobian of the result w.r.t. the heat flux
             * \param &dRdq_end: The stopping iterator of the Jacobian of the result w.r.t. the heat flux
             */

            for ( auto r = result_begin; r != result_end; r++ ){

                const unsigned int phase = ( unsigned int )( r - result_begin );

                computeBalanceOfEnergyDivergence( test_function_gradient_begin,   test_function_gradient_end,
                                                  heat_flux_begin + dim * phase,  heat_flux_begin + dim * ( phase + 1 ),
                                                  *r,
                                                  dRdGradPsi_begin + dim * phase, dRdGradPsi_begin + dim * ( phase + 1 ),
                                                  dRdq_begin + dim * phase,       dRdq_begin + dim * ( phase + 1 ) );

            }

        }

    }

}
