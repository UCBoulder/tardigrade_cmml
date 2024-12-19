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
#include<algorithm>
#include<functional>

namespace tardigradeBalanceEquations{

    namespace balanceOfLinearMomentum{

        template<
            int dim, typename density_type, typename density_dot_type,
            class density_gradient_iter, class velocity_iter, class velocity_dot_iter,
            class velocity_gradient_iter, class body_force_iter, class result_iter
        >
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            result_iter result_begin,                              result_iter result_end
        ){
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
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             */

            // Compute the inner product of the velocity and the density gradient
            using density_gradient_type  = typename std::iterator_traits<density_gradient_iter>::value_type;
            using velocity_type          = typename std::iterator_traits<velocity_iter>::value_type;
            using grad_rho_dot_type      = decltype( std::declval<density_gradient_type&>( ) * std::declval<velocity_type>( ) );
            using velocity_gradient_type = typename std::iterator_traits<velocity_gradient_iter>::value_type;

            grad_rho_dot_type grad_rho_dot_v = std::inner_product( density_gradient_begin, density_gradient_end, velocity_begin, 0. );

            // Compute the trace of the velocity gradient
            velocity_gradient_type trace_velocity_gradient = 0;

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

        template<
            int dim, typename density_type, typename density_dot_type, typename testFunction_type,
            class density_gradient_iter, class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter, class result_iter
        >
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            const testFunction_type &psi,
            result_iter result_begin,                              result_iter result_end
        ){

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
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             */

            using result_type = typename std::iterator_traits<result_iter>::value_type;
            using final_type  = decltype( std::declval<result_type&>( ) * std::declval<testFunction_type&>( ) );

            computeBalanceOfLinearMomentumNonDivergence<dim>(
                density, density_dot, density_gradient_begin, density_gradient_end,
                velocity_begin, velocity_end, velocity_dot_begin, velocity_dot_end,
                velocity_gradient_begin, velocity_gradient_end,
                body_force_begin, body_force_end,
                result_begin, result_end
            );

            std::transform(
                result_begin, result_end, result_begin,
                std::bind(
                    std::multiplies< typename std::iterator_traits<final_type>::value_type >( ),
                    std::placeholders::_1,
                    psi
                )
            );

        }

        template<
            int dim, typename density_type, typename density_dot_type,
            class density_gradient_iter, class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter,
            class result_iter, class dRdRho_iter, class dRdRhoDot_iter, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out
        >
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            result_iter result_begin,                              result_iter result_end,
            dRdRho_iter dRdRho_begin,                              dRdRho_iter dRdRho_end,
            dRdRhoDot_iter dRdRhoDot_begin,                        dRdRhoDot_iter dRdRhoDot_end,
            secondOrderTensor_iter_out dRdGradRho_begin,           secondOrderTensor_iter_out dRdGradRho_end,
            secondOrderTensor_iter_out dRdV_begin,                 secondOrderTensor_iter_out dRdV_end,
            secondOrderTensor_iter_out dRdVDot_begin,              secondOrderTensor_iter_out dRdVDot_end,
            thirdOrderTensor_iter_out dRdGradV_begin,              thirdOrderTensor_iter_out dRdGradV_end,
            secondOrderTensor_iter_out dRdB_begin,                 secondOrderTensor_iter_out dRdB_end
        ){
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
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             * \param &dRdRho_begin: The starting iterator of the derivative of the result w.r.t. the density
             * \param &dRdRho_end: The stopping iterator of the derivative of the result w.r.t. the density
             * \param &dRdRhoDot_begin: The starting iterator of the derivative of the result w.r.t. the partial temporal derivative of the density
             * \param &dRdRhoDot_end: The stopping iterator of the derivative of the result w.r.t. the partial temporal derivative of the density
             * \param &dRdGradRho_begin: The starting iterator of the derivative of the result w.r.t. the spatial derivative of the density
             * \param &dRdGradRho_end: The stopping iterator of the derivative of the result w.r.t. the spatial derivative of the density
             * \param &dRdV_begin: The starting iterator of the derivative of the result w.r.t. the velocity
             * \param &dRdV_end: The stopping iterator of the derivative of the result w.r.t. the velocity
             * \param &dRdVDot_begin: The starting iterator of the derivative of the result w.r.t. the partial temporal derivative of the velocity
             * \param &dRdVDot_end: The stopping iterator of the derivative of the result w.r.t. the partial temporal derivative of the velocity
             * \param &dRdGradV_begin: The starting iterator of the derivative of the result w.r.t. the spatial derivative of the velocity
             * \param &dRdGradV_end: The stopping iterator of the derivative of the result w.r.t. the spatial derivative of the velocity
             * \param &dRdB_begin: The starting iterator of the derivative of the result w.r.t. the body force
             * \param &dRdB_end: The stopping iterator of the derivative of the result w.r.t. the body force
             */

            // Compute the inner product of the velocity and the density gradient
            using density_gradient_type  = typename std::iterator_traits<density_gradient_iter>::value_type;
            using velocity_type          = typename std::iterator_traits<velocity_iter>::value_type;
            using grad_rho_dot_type      = decltype( std::declval<density_gradient_type&>( ) * std::declval<velocity_type>( ) );
            using velocity_gradient_type = typename std::iterator_traits<velocity_gradient_iter>::value_type;

            grad_rho_dot_type grad_rho_dot_v = std::inner_product( density_gradient_begin, density_gradient_end, velocity_begin, 0. );

            // Compute the trace of the velocity gradient
            velocity_gradient_type trace_velocity_gradient = 0;

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

        template<
            int dim, class density_iter, class density_dot_iter, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out>
        void computeBalanceOfLinearMomentumNonDivergence( const density_iter &density_begin,                     const density_iter &density_end,
                                                          const density_dot_iter &density_dot_begin,             const density_dot_iter &density_dot_end,
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
             * \param &density_begin: The starting iterator of the mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_end: The stopping iterator of the mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \left( \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \left( \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density dot vectors are not the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_dot_end - velocity_dot_begin ), "The density and velocity dot vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( body_force_end - body_force_begin ), "The density and body force vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( result_end - result_begin ), "The density and result vectors are of inconsistent sizes" );

            for ( auto rho = density_begin; rho != density_end; rho++ ){

                unsigned int phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfLinearMomentumNonDivergence<dim>(
                    *( density_begin + phase ),
                    *( density_dot_begin + phase ),
                    density_gradient_begin + dim * phase,        density_gradient_begin + dim * ( phase + 1 ),
                    velocity_begin + dim * phase,                velocity_begin + dim * ( phase + 1 ),
                    velocity_dot_begin + dim * phase,            velocity_dot_begin + dim * ( phase + 1 ),
                    velocity_gradient_begin + dim * dim * phase, velocity_gradient_begin + dim * dim * ( phase + 1 ),
                    body_force_begin + dim * phase,              body_force_begin + dim * dim * ( phase + 1 ),
                    result_begin + dim * phase,                  result_begin + dim * ( phase + 1 )
                );

            }

        }

        template<int dim, class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out>
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
             * \param &density_begin: The starting iterator of the mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_end: The stopping iterator of the mass density per unit current volume \f$ \left( \rho \right) \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \left( \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \left( \frac{\partial}{\partial t} \rho \right) \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \left( \rho_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ \left( v_i \right) \f$
             * \param &velocity_dot_begin: The starting iterator of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_dot_end: The stopping iterator of the partial time derivative of the velocity \f$ \left( \frac{\partial}{\partial t} v_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ \left( v_{i,j} \right) \f$
             * \param &body_force_begin: The starting iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &body_force_end: The stopping iterator of the body force per unit mass \f$ \left( b_i \right) \f$
             * \param &result_begin: The starting iterator of the non-divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the non-divergence part of the balance of linear momentum
             * \param &dRdRho_begin: The starting iterator of the derivative of the result w.r.t. the density
             * \param &dRdRho_end: The stopping iterator of the derivative of the result w.r.t. the density
             * \param &dRdRhoDot_begin: The starting iterator of the derivative of the result w.r.t. the partial temporal derivative of the density
             * \param &dRdRhoDot_end: The stopping iterator of the derivative of the result w.r.t. the partial temporal derivative of the density
             * \param &dRdGradRho_begin: The starting iterator of the derivative of the result w.r.t. the spatial derivative of the density
             * \param &dRdGradRho_end: The stopping iterator of the derivative of the result w.r.t. the spatial derivative of the density
             * \param &dRdV_begin: The starting iterator of the derivative of the result w.r.t. the velocity
             * \param &dRdV_end: The stopping iterator of the derivative of the result w.r.t. the velocity
             * \param &dRdVDot_begin: The starting iterator of the derivative of the result w.r.t. the partial temporal derivative of the velocity
             * \param &dRdVDot_end: The stopping iterator of the derivative of the result w.r.t. the partial temporal derivative of the velocity
             * \param &dRdGradV_begin: The starting iterator of the derivative of the result w.r.t. the spatial derivative of the velocity
             * \param &dRdGradV_end: The stopping iterator of the derivative of the result w.r.t. the spatial derivative of the velocity
             * \param &dRdB_begin: The starting iterator of the derivative of the result w.r.t. the body force
             * \param &dRdB_end: The stopping iterator of the derivative of the result w.r.t. the body force
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(                   ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density dot vectors are not the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_dot_end - velocity_dot_begin ), "The density and velocity dot vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( body_force_end - body_force_begin ), "The density and body force vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( result_end - result_begin ), "The density and result vectors are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdRho_end - dRdRho_begin ), "The density and dRdRho are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdRhoDot_end - dRdRhoDot_begin ), "The density and dRdRhoDot are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdGradRho_end - dRdGradRho_begin ), "The density and dRdGradRho are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdV_end - dRdV_begin ), "The density and dRdV are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdGradV_end - dRdGradV_begin ), "The density and dRdGradV are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdB_end - dRdB_begin ), "The density and dRdB are of inconsistent sizes" );

            for ( auto rho = density_begin; rho != density_end; rho++ ){

                unsigned int phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfLinearMomentumNonDivergence<dim>(
                    *( density_begin + phase ),
                    *( density_dot_begin + phase ),
                    density_gradient_begin + dim * phase,        density_gradient_begin + dim * ( phase + 1 ),
                    velocity_begin + dim * phase,                velocity_begin + dim * ( phase + 1 ),
                    velocity_dot_begin + dim * phase,            velocity_dot_begin + dim * ( phase + 1 ),
                    velocity_gradient_begin + dim * dim * phase, velocity_gradient_begin + dim * dim * ( phase + 1 ),
                    body_force_begin + dim * phase,              body_force_begin + dim * dim * ( phase + 1 ),
                    result_begin + dim * phase,                  result_begin + dim * ( phase + 1 ),
                    dRdRho_begin + dim * phase,                  dRdRho_begin + dim * ( phase + 1 ),
                    dRdRhoDot_begin + dim * phase,               dRdRhoDot_begin + dim * ( phase + 1 ),
                    dRdGradRho_begin + dim * dim * phase,        dRdGradRho_begin + dim * dim * ( phase + 1 ),
                    dRdV_begin + dim * dim * phase,              dRdV_begin + dim * dim * ( phase + 1 ),
                    dRdVDot_begin + dim * dim * phase,           dRdVDot_begin + dim * dim * ( phase + 1 ),
                    dRdGradV_begin + dim * dim * dim * phase,    dRdGradV_begin + dim * dim * dim * ( phase + 1 ),
                    dRdB_begin + dim * dim * phase,              dRdB_begin + dim * dim * ( phase + 1 ) );

            }

        }

        template<int dim, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out>
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
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &volume_fraction: The volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &result_begin: The starting iterator of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the divergence part of the balance of linear momentum
             */

            std::fill( result_begin, result_end, 0. );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    *( result_begin + i ) -= ( *( test_function_gradient_begin + j ) ) * volume_fraction * ( *( cauchy_stress_begin + dim * j + i ) );

                }

            }

        }

        template<int dim, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out>
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
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &volume_fraction: The volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &result_begin: The starting iterator of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the divergence part of the balance of linear momentum
             * \param &dRdGradPsi_begin: The starting iterator of the derivative of the result w.r.t. the test function spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdGradPsi_end: The stopping iterator of the derivative of the result w.r.t. the test function spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdCauchy_begin: The starting iterator of the derivative of the result w.r.t. the Cauchy stress \f$ \left( \bf{\sigma} \right) \f$
             * \param &dRdCauchy_end: The stopping iterator of the derivative of the result w.r.t. the Cauchy stress \f$ \left( \bf{\sigma} \right) \f$
             * \param &dRdPhi_begin: The starting iterator of the derivative of the result w.r.t. the volume fraction \f$ \left( \bf{\phi} \right) \f$
             * \param &dRdPhi_end: The stopping iterator of the derivative of the result w.r.t. the volume fraction \f$ \left( \bf{\phi} \right) \f$
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

        template<int dim, class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out>
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
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &result_begin: The starting iterator of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the divergence part of the balance of linear momentum
             */

            for ( auto phi = volume_fraction_begin; phi != volume_fraction_end; phi++ ){

                unsigned int phase = ( unsigned int )( phi - volume_fraction_begin );

                computeBalanceOfLinearMomentumDivergence<dim>(
                    test_function_gradient_begin,            test_function_gradient_end,
                    cauchy_stress_begin + dim * dim * phase, cauchy_stress_begin + dim * dim * ( phase + 1 ),
                    *phi,
                    result_begin + dim * phase,              result_begin + dim * ( phase + 1 )
                );

            }

        }

        template<int dim, class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out>
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
             * \param &test_function_gradient_begin: The starting iterator of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \left( \bar{\bf{\sigma}} \right) \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of the current phase \f$ \left( \phi \right) \f$
             * \param &result_begin: The starting iterator of the divergence part of the balance of linear momentum
             * \param &result_end: The stopping iterator of the divergence part of the balance of linear momentum
             * \param &dRdGradPsi_begin: The starting iterator of the derivative of the result w.r.t. the test function spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdGradPsi_end: The stopping iterator of the derivative of the result w.r.t. the test function spatial gradient \f$ \left( \psi_{,i} \right) \f$
             * \param &dRdCauchy_begin: The starting iterator of the derivative of the result w.r.t. the Cauchy stress \f$ \left( \bf{\sigma} \right) \f$
             * \param &dRdCauchy_end: The stopping iterator of the derivative of the result w.r.t. the Cauchy stress \f$ \left( \bf{\sigma} \right) \f$
             * \param &dRdPhi_begin: The starting iterator of the derivative of the result w.r.t. the volume fraction \f$ \left( \bf{\phi} \right) \f$
             * \param &dRdPhi_end: The stopping iterator of the derivative of the result w.r.t. the volume fraction \f$ \left( \bf{\phi} \right) \f$
             */

            for ( auto phi = volume_fraction_begin; phi != volume_fraction_end; phi++ ){

                unsigned int phase = ( unsigned int )( phi - volume_fraction_begin );

                computeBalanceOfLinearMomentumDivergence<dim>(
                    test_function_gradient_begin,            test_function_gradient_end,
                    cauchy_stress_begin + dim * dim * phase,   cauchy_stress_begin + dim * dim * ( phase + 1 ),
                    *phi,
                    result_begin + dim * phase,                result_begin + dim * ( phase + 1 ),
                    dRdGradPsi_begin + dim * dim * phase,      dRdGradPsi_begin + dim * dim * ( phase + 1 ),
                    dRdCauchy_begin + dim * dim * dim * phase, dRdCauchy_begin + dim * dim * dim * ( phase + 1 ),
                    dRdPhi_begin + dim * phase,                dRdPhi_begin + dim * ( phase + 1 )
                );

            }

        }

    }

}
