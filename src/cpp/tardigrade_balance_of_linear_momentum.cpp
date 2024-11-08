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

        template<class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out>
        void computeBalanceOfLinearMomentumNonDivergence( const floatType &density, const floatType &density_dot,
                                                          const floatVector_iter &density_gradient_begin,        const floatVector_iter &density_gradient_end,
                                                          const floatVector_iter &velocity_begin,                const floatVector_iter &velocity_end,
                                                          const floatVector_iter &velocity_dot_begin,            const floatVector_iter &velocity_dot_end,
                                                          const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                                          const floatVector_iter &body_force_begin,              const floatVector_iter &body_force_end,
                                                          floatVector_iter_out result_begin,                     floatVector_iter_out result_end, 
                                                          floatVector_iter_out dRdRho_begin,                     floatVector_iter_out dRdRho_end,
                                                          floatVector_iter_out dRdRhoDot_begin,                  floatVector_iter_out dRdRhoDot_end,
                                                          secondOrderTensor_iter_out dRdGradRho_begin,           secondOrderTensor_iter_out dRdGradRho_end,
                                                          secondOrderTensor_iter_out dRdV_begin,                 secondOrderTensor_iter_out dRdV_end,
                                                          secondOrderTensor_iter_out dRdVDot_begin,              secondOrderTensor_iter_out dRdVDot_end,
                                                          thirdOrderTensor_iter_out dRdGradV_begin,              thirdOrderTensor_iter_out dRdGradV_end,
                                                          secondOrderTensor_iter_out dRdB_begin,                 secondOrderTensor_iter_out dRdB_end ){
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
             * \param &dRdRho_begin: The starting point of the derivative of the result w.r.t. the density
             * \param &dRdRho_end: The stopping point of the derivative of the result w.r.t. the density
             * \param &dRdRhoDot_begin: The starting point of the derivative of the result w.r.t. the partial temporal derivative of the density
             * \param &dRdRhoDot_end: The stopping point of the derivative of the result w.r.t. the partial temporal derivative of the density
             * \param &dRdGradRho_begin: The starting point of the derivative of the result w.r.t. the spatial derivative of the density
             * \param &dRdGradRho_end: The stopping point of the derivative of the result w.r.t. the spatial derivative of the density
             * \param &dRdV_begin: The starting point of the derivative of the result w.r.t. the velocity
             * \param &dRdV_end: The stopping point of the derivative of the result w.r.t. the velocity
             * \param &dRdVDot_begin: The starting point of the derivative of the result w.r.t. the partial temporal derivative of the velocity
             * \param &dRdVDot_end: The stopping point of the derivative of the result w.r.t. the partial temporal derivative of the velocity
             * \param &dRdGradV_begin: The starting point of the derivative of the result w.r.t. the spatial derivative of the velocity
             * \param &dRdGradV_end: The stopping point of the derivative of the result w.r.t. the spatial derivative of the velocity
             * \param &dRdB_begin: The starting point of the derivative of the result w.r.t. the body force
             * \param &dRdB_end: The stopping point of the derivative of the result w.r.t. the body force
             */

            // Compute the inner product of the velocity and the density gradient
            floatType grad_rho_dot_v = std::inner_product( density_gradient_begin, density_gradient_end, velocity_begin, 0. );

            // Compute the trace of the velocity gradient
            floatType trace_velocity_gradient = 0;

            for ( unsigned int i = 0; i < dim; i++ ){

                trace_velocity_gradient += *( velocity_gradient_begin + dim * i + i );

            }

            // Initialize the gradients
            std::fill( dRdGradRho_begin, dRdGradRho_end, 0. );
            std::fill( dRdV_begin, dRdV_end, 0. );
            std::fill( dRdVDot_begin, dRdVDot_end, 0. );
            std::fill( dRdGradV_begin, dRdGradV_end, 0. );
            std::fill( dRdB_begin, dRdB_end, 0. );

            // Compute the balance of linear momentum

            for ( unsigned int i = 0; i < dim; i++ ){

                *( result_begin + i ) = density * ( ( *( body_force_begin + i) ) - ( *( velocity_dot_begin + i ) ) - trace_velocity_gradient * ( *( velocity_begin + i ) ) ) - grad_rho_dot_v * ( *( velocity_begin + i ) )
                                      - density_dot * ( *( velocity_begin + i ) );

                *( dRdRho_begin + i ) = ( ( *( body_force_begin + i) ) - ( *( velocity_dot_begin + i ) ) - trace_velocity_gradient * ( *( velocity_begin + i ) ) );

                *( dRdRhoDot_begin + i ) = -( *( velocity_begin + i ) );

                *( dRdV_begin + dim * i + i ) = -density * trace_velocity_gradient - grad_rho_dot_v - density_dot;

                *( dRdVDot_begin + dim * i + i ) = -density;

                *( dRdB_begin + dim * i + i ) = density;

                for ( unsigned int j = 0; j < dim; j++ ){

                    *( result_begin + i ) -= density * ( *( velocity_gradient_begin + dim * i + j ) ) * ( *( velocity_begin + j ) );

                    *( dRdRho_begin + i ) -= ( *( velocity_gradient_begin + dim * i + j ) ) * ( *( velocity_begin + j ) );

                    *( dRdGradRho_begin + dim * i + j ) = -( *( velocity_begin + i ) ) * ( *( velocity_begin + j ) );

                    *( dRdV_begin + dim * i + j ) -= ( *( density_gradient_begin + j ) ) * ( *( velocity_begin + i ) ) + density * ( *( velocity_gradient_begin + dim * i + j ) );

                    *( dRdGradV_begin + dim * dim * i + dim * j + j ) -= density * ( *( velocity_begin + i ) );

                    *( dRdGradV_begin + dim * dim * i + dim * i + j ) -= density * ( *( velocity_begin + j ) );

                }

            }

        }

        template<class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out>
        void computeBalanceOfLinearMomentumNonDivergence( const scalarArray_iter &density_begin,                 const scalarArray_iter &density_end,
                                                          const scalarArray_iter &density_dot_begin,             const scalarArray_iter &density_dot_end,
                                                          const floatVector_iter &density_gradient_begin,        const floatVector_iter &density_gradient_end,
                                                          const floatVector_iter &velocity_begin,                const floatVector_iter &velocity_end,
                                                          const floatVector_iter &velocity_dot_begin,            const floatVector_iter &velocity_dot_end,
                                                          const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                                          const floatVector_iter &body_force_begin,              const floatVector_iter &body_force_end,
                                                          floatVector_iter_out result_begin,                     floatVector_iter_out result_end ){
            /*!
             * Compute the non-divergence part of the multiphase balance of linear momentum where the full equation is
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
             * \param &density_begin: The starting point of the mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_end: The stopping point of the mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot_begin: The starting point of the partial time derivative of the density \f$ \left( \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_dot_end: The stopping point of the partial time derivative of the density \f$ \left( \frac{\partial}{\partial t} \rho \right) \f$
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

            for ( auto rho = density_begin; rho != density_end; rho++ ){

                unsigned int phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfLinearMomentumNonDivergence( *( density_begin + phase ),
                                                             *( density_dot_begin + phase ),
                                                             density_gradient_begin + dim * phase,      density_gradient_begin + dim * ( phase + 1 ),
                                                             velocity_begin + dim * phase,              velocity_begin + dim * ( phase + 1 ),
                                                             velocity_dot_begin + dim * phase,          velocity_dot_begin + dim * ( phase + 1 ),
                                                             velocity_gradient_begin + sot_dim * phase, velocity_gradient_begin + sot_dim * ( phase + 1 ),
                                                             body_force_begin + dim * phase,            body_force_begin + sot_dim * ( phase + 1 ),
                                                             result_begin + dim * phase,                result_begin + dim * ( phase + 1 ) );

            }

        }

        template<class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out>
        void computeBalanceOfLinearMomentumNonDivergence( const scalarArray_iter &density_begin,                 const scalarArray_iter &density_end,
                                                          const scalarArray_iter &density_dot_begin,             const scalarArray_iter &density_dot_end,
                                                          const floatVector_iter &density_gradient_begin,        const floatVector_iter &density_gradient_end,
                                                          const floatVector_iter &velocity_begin,                const floatVector_iter &velocity_end,
                                                          const floatVector_iter &velocity_dot_begin,            const floatVector_iter &velocity_dot_end,
                                                          const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
                                                          const floatVector_iter &body_force_begin,              const floatVector_iter &body_force_end,
                                                          floatVector_iter_out result_begin,                     floatVector_iter_out result_end, 
                                                          floatVector_iter_out dRdRho_begin,                     floatVector_iter_out dRdRho_end,
                                                          floatVector_iter_out dRdRhoDot_begin,                  floatVector_iter_out dRdRhoDot_end,
                                                          secondOrderTensor_iter_out dRdGradRho_begin,           secondOrderTensor_iter_out dRdGradRho_end,
                                                          secondOrderTensor_iter_out dRdV_begin,                 secondOrderTensor_iter_out dRdV_end,
                                                          secondOrderTensor_iter_out dRdVDot_begin,              secondOrderTensor_iter_out dRdVDot_end,
                                                          thirdOrderTensor_iter_out dRdGradV_begin,              thirdOrderTensor_iter_out dRdGradV_end,
                                                          secondOrderTensor_iter_out dRdB_begin,                 secondOrderTensor_iter_out dRdB_end ){
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
             * \param &density_begin: The starting point of the mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_end: The stopping point of the mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot_begin: The starting point of the partial time derivative of the density \f$ \left( \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_dot_end: The stopping point of the partial time derivative of the density \f$ \left( \frac{\partial}{\partial t} \rho \right) \f$
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
             * \param &dRdRho_begin: The starting point of the derivative of the result w.r.t. the density
             * \param &dRdRho_end: The stopping point of the derivative of the result w.r.t. the density
             * \param &dRdRhoDot_begin: The starting point of the derivative of the result w.r.t. the partial temporal derivative of the density
             * \param &dRdRhoDot_end: The stopping point of the derivative of the result w.r.t. the partial temporal derivative of the density
             * \param &dRdGradRho_begin: The starting point of the derivative of the result w.r.t. the spatial derivative of the density
             * \param &dRdGradRho_end: The stopping point of the derivative of the result w.r.t. the spatial derivative of the density
             * \param &dRdV_begin: The starting point of the derivative of the result w.r.t. the velocity
             * \param &dRdV_end: The stopping point of the derivative of the result w.r.t. the velocity
             * \param &dRdVDot_begin: The starting point of the derivative of the result w.r.t. the partial temporal derivative of the velocity
             * \param &dRdVDot_end: The stopping point of the derivative of the result w.r.t. the partial temporal derivative of the velocity
             * \param &dRdGradV_begin: The starting point of the derivative of the result w.r.t. the spatial derivative of the velocity
             * \param &dRdGradV_end: The stopping point of the derivative of the result w.r.t. the spatial derivative of the velocity
             * \param &dRdB_begin: The starting point of the derivative of the result w.r.t. the body force
             * \param &dRdB_end: The stopping point of the derivative of the result w.r.t. the body force
             */

            for ( auto rho = density_begin; rho != density_end; rho++ ){

                unsigned int phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfLinearMomentumNonDivergence( *( density_begin + phase ),
                                                             *( density_dot_begin + phase ),
                                                             density_gradient_begin + dim * phase,      density_gradient_begin + dim * ( phase + 1 ),
                                                             velocity_begin + dim * phase,              velocity_begin + dim * ( phase + 1 ),
                                                             velocity_dot_begin + dim * phase,          velocity_dot_begin + dim * ( phase + 1 ),
                                                             velocity_gradient_begin + sot_dim * phase, velocity_gradient_begin + sot_dim * ( phase + 1 ),
                                                             body_force_begin + dim * phase,            body_force_begin + sot_dim * ( phase + 1 ),
                                                             result_begin + dim * phase,                result_begin + dim * ( phase + 1 ),
                                                             dRdRho_begin + dim * phase,                dRdRho_begin + dim * ( phase + 1 ),
                                                             dRdRhoDot_begin + dim * phase,             dRdRhoDot_begin + dim * ( phase + 1 ),
                                                             dRdGradRho_begin + sot_dim * phase,        dRdGradRho_begin + sot_dim * ( phase + 1 ),
                                                             dRdV_begin + sot_dim * phase,              dRdV_begin + sot_dim * ( phase + 1 ),
                                                             dRdVDot_begin + sot_dim * phase,           dRdVDot_begin + sot_dim * ( phase + 1 ),
                                                             dRdGradV_begin + dim * dim * dim * phase,  dRdGradV_begin + dim * dim * dim * ( phase + 1 ),
                                                             dRdB_begin + sot_dim * phase,              dRdB_begin + sot_dim * ( phase + 1 ) );

            }

        }

        template<class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out>
        void computeBalanceOfLinearMomentumDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                                       const secondOrderTensor_iter &cauchy_stress_begin,    const secondOrderTensor_iter &cauchy_stress_end,
                                                       const floatType &volume_fraction,
                                                       floatVector_iter_out result_begin,                    floatVector_iter_out result_end ){
            /*!
             * Compute the divergence part of the balance of linear momentum where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho v_i \right) = -f_i \f$
             * 
             * and the divergence part is
             * 
             * \f$ \left( \sigma_{ji} \right)_{,j} \f$
             * 
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. In the context of a variationally-based solution technique
             * we multiply this by a test function \f$ \psi \f$ and write
             * 
             * \f$ \left( \psi \sigma_{ji} \right)_{,j} - \psi_{,j} \sigma_{ji} \f$. The divergence term is defined as the second term. We also introduce a volume fraction
             * \f$ \phi \f$ to prepare for multi-phase continuum such that we will actually compute
             * 
             * \f$ -\psi_{,j} \phi \bar{\sigma}_{ji} \f$
             *
             * where we made the definition that the effective Cauchy stress \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$ where \f$\bar{\bf{\sigma}} \f$ is the true
             * Cauchy stress of the phase.
             * 
             * \param &test_function_gradient_begin: The starting point of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping point of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting point of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &cauchy_stress_end: The stopping point of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &volume_fraction: The volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &result_begin: The starting point of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping point of the divergence part of the balance of linear momentum
             */

            std::fill( result_begin, result_end, 0. );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    *( result_begin + i ) -= ( *( test_function_gradient_begin + j ) ) * volume_fraction * ( *( cauchy_stress_begin + dim * j + i ) );

                }

            }

        }

        template<class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out>
        void computeBalanceOfLinearMomentumDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                                       const secondOrderTensor_iter &cauchy_stress_begin,    const secondOrderTensor_iter &cauchy_stress_end,
                                                       const floatType &volume_fraction,
                                                       floatVector_iter_out result_begin,                    floatVector_iter_out result_end,
                                                       secondOrderTensor_iter_out dRdGradPsi_begin,          secondOrderTensor_iter_out dRdGradPsi_end,
                                                       thirdOrderTensor_iter_out dRdCauchy_begin,            thirdOrderTensor_iter_out dRdCauchy_end,
                                                       floatVector_iter_out dRdPhi_begin,                    floatVector_iter_out dRdPhi_end ){
            /*!
             * Compute the divergence part of the balance of linear momentum and its Jacobians where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho v_i \right) = -f_i \f$
             * 
             * and the divergence part is
             * 
             * \f$ \left( \sigma_{ji} \right)_{,j} \f$
             * 
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. In the context of a variationally-based solution technique
             * we multiply this by a test function \f$ \psi \f$ and write
             * 
             * \f$ \left( \psi \sigma_{ji} \right)_{,j} - \psi_{,j} \sigma_{ji} \f$. The divergence term is defined as the second term. We also introduce a volume fraction
             * \f$ \phi \f$ to prepare for multi-phase continuum such that we will actually compute
             * 
             * \f$ -\psi_{,j} \phi \bar{\sigma}_{ji} \f$
             *
             * where we made the definition that the effective Cauchy stress \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$ where \f$\bar{\bf{\sigma}} \f$ is the true
             * Cauchy stress of the phase.
             * 
             * \param &test_function_gradient_begin: The starting point of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping point of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting point of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &cauchy_stress_end: The stopping point of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &volume_fraction: The volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &result_begin: The starting point of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping point of the divergence part of the balance of linear momentum
             * \param &dRdGradPsi_begin: The starting point of the derivative of the result w.r.t. the test function spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdGradPsi_end: The stopping point of the derivative of the result w.r.t. the test function spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdCauchy_begin: The starting point of the derivative of the result w.r.t. the Cauchy stress \f$ \left( \bf{\sigma} \right) \f$
             * \param &dRdCauchy_end: The stopping point of the derivative of the result w.r.t. the Cauchy stress \f$ \left( \bf{\sigma} \right) \f$
             * \param &dRdPhi_begin: The starting point of the derivative of the result w.r.t. the volume fraction \f$ \left( \bf{\phi} \right) \f$
             * \param &dRdPhi_end: The stopping point of the derivative of the result w.r.t. the volume fraction \f$ \left( \bf{\phi} \right) \f$
             */

            std::fill( result_begin, result_end, 0. );
            std::fill( dRdGradPsi_begin, dRdGradPsi_end, 0. );
            std::fill( dRdCauchy_begin, dRdCauchy_end, 0. );
            std::fill( dRdPhi_begin, dRdPhi_end, 0. );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    *( result_begin + i ) -= ( *( test_function_gradient_begin + j ) ) * volume_fraction * ( *( cauchy_stress_begin + dim * j + i ) );

                    *( dRdGradPsi_begin + dim * i + j ) -= volume_fraction * ( *( cauchy_stress_begin + dim * j + i ) );

                    *( dRdCauchy_begin + dim * dim * i + dim * j + i ) -= ( *( test_function_gradient_begin + j ) ) * volume_fraction;

                    *( dRdPhi_begin + i ) -= ( *( test_function_gradient_begin + j ) ) * ( *( cauchy_stress_begin + dim * j + i ) );

                }

            }

        }

        template<class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out>
        void computeBalanceOfLinearMomentumDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                                       const secondOrderTensor_iter &cauchy_stress_begin,    const secondOrderTensor_iter &cauchy_stress_end,
                                                       const scalarArray_iter &volume_fraction_begin,        const scalarArray_iter &volume_fraction_end,
                                                       floatVector_iter_out result_begin,                    floatVector_iter_out result_end ){
            /*!
             * Compute the divergence part of the balance of linear momentum where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho v_i \right) = -f_i \f$
             * 
             * and the divergence part is
             * 
             * \f$ \left( \sigma_{ji} \right)_{,j} \f$
             * 
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. In the context of a variationally-based solution technique
             * we multiply this by a test function \f$ \psi \f$ and write
             * 
             * \f$ \left( \psi \sigma_{ji} \right)_{,j} - \psi_{,j} \sigma_{ji} \f$. The divergence term is defined as the second term. We also introduce a volume fraction
             * \f$ \phi \f$ to prepare for multi-phase continuum such that we will actually compute
             * 
             * \f$ -\psi_{,j} \phi \bar{\sigma}_{ji} \f$
             *
             * where we made the definition that the effective Cauchy stress \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$ where \f$\bar{\bf{\sigma}} \f$ is the true
             * Cauchy stress of the phase.
             * 
             * \param &test_function_gradient_begin: The starting point of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping point of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting point of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &cauchy_stress_end: The stopping point of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &volume_fraction_begin: The starting point of the volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &volume_fraction_end: The stopping point of the volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &result_begin: The starting point of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping point of the divergence part of the balance of linear momentum
             */

            for ( auto phi = volume_fraction_begin; phi != volume_fraction_end; phi++ ){

                unsigned int phase = ( unsigned int )( phi - volume_fraction_begin );

                computeBalanceOfLinearMomentumDivergence(  test_function_gradient_begin + dim * phase, test_function_gradient_begin + dim * ( phase + 1 ),
                                                           cauchy_stress_begin + sot_dim * phase,      cauchy_stress_begin + sot_dim * ( phase + 1 ),
                                                          *phi,
                                                           result_begin + dim * phase,                 result_begin + dim * ( phase + 1 ) );

            }

        }

        template<class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out>
        void computeBalanceOfLinearMomentumDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                                       const secondOrderTensor_iter &cauchy_stress_begin,    const secondOrderTensor_iter &cauchy_stress_end,
                                                       const scalarArray_iter &volume_fraction_begin,        const scalarArray_iter &volume_fraction_end,
                                                       floatVector_iter_out result_begin,                    floatVector_iter_out result_end,
                                                       secondOrderTensor_iter_out dRdGradPsi_begin,          secondOrderTensor_iter_out dRdGradPsi_end,
                                                       thirdOrderTensor_iter_out dRdCauchy_begin,            thirdOrderTensor_iter_out dRdCauchy_end,
                                                       floatVector_iter_out dRdPhi_begin,                    floatVector_iter_out dRdPhi_end ){
            /*!
             * Compute the divergence part of the balance of linear momentum and its Jacobians where the full equation is
             *
             * \f$ \left( \sigma_{ji} - \rho v_i v_j \right)_{,j} + \rho b_i - \frac{\partial}{\partial t} \left( \rho v_i \right) = -f_i \f$
             * 
             * and the divergence part is
             * 
             * \f$ \left( \sigma_{ji} \right)_{,j} \f$
             * 
             * and \f$ f_i \f$ is the vector of additional forces not accounted for in these equations. In the context of a variationally-based solution technique
             * we multiply this by a test function \f$ \psi \f$ and write
             * 
             * \f$ \left( \psi \sigma_{ji} \right)_{,j} - \psi_{,j} \sigma_{ji} \f$. The divergence term is defined as the second term. We also introduce a volume fraction
             * \f$ \phi \f$ to prepare for multi-phase continuum such that we will actually compute
             * 
             * \f$ -\psi_{,j} \phi \bar{\sigma}_{ji} \f$
             *
             * where we made the definition that the effective Cauchy stress \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$ where \f$\bar{\bf{\sigma}} \f$ is the true
             * Cauchy stress of the phase.
             * 
             * \param &test_function_gradient_begin: The starting point of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping point of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting point of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &cauchy_stress_end: The stopping point of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &volume_fraction_begin: The starting point of the volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &volume_fraction_end: The stopping point of the volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &result_begin: The starting point of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping point of the divergence part of the balance of linear momentum
             * \param &dRdGradPsi_begin: The starting point of the derivative of the result w.r.t. the test function spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdGradPsi_end: The stopping point of the derivative of the result w.r.t. the test function spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdCauchy_begin: The starting point of the derivative of the result w.r.t. the Cauchy stress \f$ \left( \bf{\sigma} \right) \f$
             * \param &dRdCauchy_end: The stopping point of the derivative of the result w.r.t. the Cauchy stress \f$ \left( \bf{\sigma} \right) \f$
             * \param &dRdPhi_begin: The starting point of the derivative of the result w.r.t. the volume fraction \f$ \left( \bf{\phi} \right) \f$
             * \param &dRdPhi_end: The stopping point of the derivative of the result w.r.t. the volume fraction \f$ \left( \bf{\phi} \right) \f$
             */

            for ( auto phi = volume_fraction_begin; phi != volume_fraction_end; phi++ ){

                unsigned int phase = ( unsigned int )( phi - volume_fraction_begin );

                computeBalanceOfLinearMomentumDivergence(  test_function_gradient_begin + dim * phase, test_function_gradient_begin + dim * ( phase + 1 ),
                                                           cauchy_stress_begin + sot_dim * phase,      cauchy_stress_begin + sot_dim * ( phase + 1 ),
                                                          *phi,
                                                           result_begin + dim * phase,                 result_begin + dim * ( phase + 1 ),
                                                           dRdGradPsi_begin + dim * dim * phase,       dRdGradPsi_begin + dim * dim * ( phase + 1 ),
                                                           dRdCauchy_begin + dim * sot_dim * phase,    dRdCauchy_begin + dim * sot_dim * ( phase + 1 ),
                                                           dRdPhi_begin + dim * phase,                 dRdPhi_begin + dim * ( phase + 1 ) );

            }

        }

    }

}
