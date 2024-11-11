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
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bm{\sigma}} \f$ where $\bm{\sigma} = \phi \bar{\bm{\sigma}} \f$
             * \param &cauchy_stress_end: The end point of the true Cauchy stress \f$ \bar{\bm{\sigma}} \f$ where $\bm{\sigma} = \phi \bar{\bm{\sigma}} \f$
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
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bm{\sigma}} \f$ where $\bm{\sigma} = \phi \bar{\bm{\sigma}} \f$
             * \param &cauchy_stress_end: The end point of the true Cauchy stress \f$ \bar{\bm{\sigma}} \f$ where $\bm{\sigma} = \phi \bar{\bm{\sigma}} \f$
             * \param &volume_fraction: The volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation: The internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &result: The result of the non-divergence part of the balance of energy
             * \param &dRdRho: The Jacobian w.r.t. the apparent density
             * \param &dRdRhoDot: The Jacobian w.r.t. the partial temporal derivative of the apparent density
             * \param &dRdGradRho_begin: The starting iterator of the Jacobian w.r.t. the spatial gradient of the apparent density
             * \param &dRdGradRho_end: The starting iterator of the Jacobian w.r.t. the spatial gradient of the apparent density
             * \param &dRdE: The Jacobian w.r.t. the internal energy
             * \param &dRdEDot: The Jacobian w.r.t. the partial temporal derivative of the internal energy
             * \param &dRdGradE_begin: The starting iterator of the Jacobian w.r.t. the spatial gradient of the internal energy
             * \param &dRdGradE_end: The starting iterator of the Jacobian w.r.t. the spatial gradient of the internal energy
             * \param &dRdV_begin: The starting iterator of the Jacobian w.r.t. the velocity
             * \param &dRdV_end: The starting iterator of the Jacobian w.r.t. the velocity
             * \param &dRdGradV_begin: The starting iterator of the Jacobian w.r.t. the velocity gradient
             * \param &dRdGradV_end: The starting iterator of the Jacobian w.r.t. the velocity gradient
             * \param &dRdCauchy_begin: The starting iterator of the Jacobian w.r.t. the Cauchy stress
             * \param &dRdCauchy_end: The starting iterator of the Jacobian w.r.t. the Cauchy stress
             * \param &dRdPhi: The Jacobian w.r.t. the volume fraction
             * \param &dRdr: The Jacobian w.r.t. the internal heat generation
             * \param &dRdpi_begin: The starting iterator of the Jacobian w.r.t. the net interphase force
             * \param &dRdpi_end: The starting iterator of the Jacobian w.r.t. the net interphase force
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

    }

}
