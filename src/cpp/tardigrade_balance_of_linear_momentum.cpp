/**
  ******************************************************************************
  * \file tardigrade_balance_of_linear_momentum.cpp
  ******************************************************************************
  * The source file for the equations associated with the balance of linear
  * momentum
  ******************************************************************************
  */

#include "tardigrade_balance_of_linear_momentum.h"
#include<numeric>

namespace tardigradeBalanceEquations{

    namespace balanceOfLinearMomentum{

        template<class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out>
        void computeBalanceOfLinearMomentumNonDivergence( const floatType &density, const floatType &density_dot,
                                                          const floatVector_iter &density_gradient_begin,        const floatVector_iter &density_gradient_end,
                                                          const floatVector_iter &velocity_begin,                const floatVector_iter &velocity_end,
                                                          const floatVector_iter &velocity_dot_begin,            const floatVector_iter &velocity_dot_end,
                                                          const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                                          const floatVector_iter &body_force_begin,              const floatVector_iter &body_force_end,
                                                          floatVector_iter_out result_begin,                     floatVector_iter_out result_end ){
            /*!
             * Compute the non-divergence part of the balance of linear momentum where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho v_i \right) = -f_i \f$
             * 
             * and the non-divergence part is
             * 
             * \f$ -\left( \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho v_i \right) + f_i = 0 \f$
             * 
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. The expression we implement is therefore
             * 
             * \f$ -\rho_{,j} v_j v_i - \rho v_{i,j} v_j - \rho v_i v_{j,j} + \rho b_i - \frac{\partial \rho}{\partial t} v_i - \rho \frac{\partial v_i}{\partial t} \f$
             * 
             * \param &density: The mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot: The partial time derivative of the density \f$ \left( \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_gradient_begin: The starting point of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping point of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting point of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping point of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting point of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping point of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting point of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping point of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &body_force_begin: The starting point of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &body_force_end: The stopping point of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &result_begin: The starting point of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping point of the non-divergence part of the balance of linear momentum
             */

            // Compute the inner product of the velocity and the density gradient
            floatType grad_rho_dot_v = std::inner_product( density_gradient_begin, density_gradient_end, velocity_begin, 0. );

            // Compute the trace of the velocity gradient
            floatType trace_velocity_gradient = 0;

            for ( unsigned int i = 0; i < dim; i++ ){

                trace_velocity_gradient += *( velocity_gradient_begin + dim * i + i );

            }

            // Compute the balance of linear momentum
            for ( unsigned int i = 0; i < dim; i++ ){

                *( result_begin + i ) = density * ( ( *( body_force_begin + i) ) - ( *( velocity_dot_begin + i ) ) - trace_velocity_gradient * ( *( velocity_begin + i ) ) ) - grad_rho_dot_v * ( *( velocity_begin + i ) )
                            - density_dot * ( *( velocity_begin + i ) );

                for ( unsigned int j = 0; j < dim; j++ ){

                    *( result_begin + i ) -= density * ( *( velocity_gradient_begin + dim * i + j ) ) * ( *( velocity_begin + j ) );

                }

            }

        }

    }

}
