/**
  ******************************************************************************
  * \file tardigrade_balance_of_mass.cpp
  ******************************************************************************
  * The source file for the equations associated with the balance of mass
  ******************************************************************************
  */

#include "tardigrade_balance_of_mass.h"
#include<numeric>
#include<algorithm>
#include<functional>

namespace tardigradeBalanceEquations{

    namespace balanceOfMass{

        void computeBalanceOfMass(
            const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
            const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &result
        ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient: The value of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity: The value of the velocity \f$ v_i \f$
             * \param &velocity_gradient: The value of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &result: The net mass change per unit volume \f$ c \f$
             */

            computeBalanceOfMass<global_dim>(
                density, density_dot,
                std::begin( density_gradient ),  std::end( density_gradient ),
                std::begin( velocity ),          std::end( velocity ),
                std::begin( velocity_gradient ), std::end( velocity_gradient ),
                result
            );

        }

        void computeBalanceOfMass(
            const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
            const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &result, 
            floatType   &dRdRho, floatType         &dRdRhoDot, floatVector &dRdGradRho,
            floatVector &dRdV,   secondOrderTensor &dRdGradV
        ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate including the Jacobians
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient: The value of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity: The value of the velocity \f$ v_i \f$
             * \param &velocity_gradient: The value of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &result: The net mass change per unit volume \f$ c \f$
             * \param &dRdRho: The derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dRdRhoDot: The derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dRdGradRho: The derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dRdV: The derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dRdGradV: The derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             */

            computeBalanceOfMass<global_dim>(
                density,  density_dot,
                std::begin( density_gradient ),  std::end( density_gradient ),
                std::begin( velocity ),          std::end( velocity ),
                std::begin( velocity_gradient ), std::end( velocity_gradient ),
                result, dRdRho, dRdRhoDot,
                std::begin( dRdGradRho ), std::end( dRdGradRho ),
                std::begin( dRdV ),       std::end( dRdV ),
                std::begin( dRdGradV ),   std::end( dRdGradV )
            );
        }

        template<
            int dim, typename density_type, typename densityDot_type, typename result_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            result_type &result
        ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &result: The net mass change per unit volume \f$ c \f$
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( ( unsigned int )( density_gradient_end - density_gradient_begin ) == ( unsigned int )( velocity_end - velocity_begin ), "The density gradient and the velocity must have the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ) == dim * dim, "The velocity gradient has a size of " + std::to_string( ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ) ) + " and must have a size of " + std::to_string( dim * dim ) );

            result = std::inner_product( density_gradient_begin, density_gradient_end,
                                                   velocity_begin, density_dot );

            for ( unsigned int i = 0; i < dim; ++i ){

                result += density * ( *( velocity_gradient_begin + dim * i + i ) );

            }

        }

        template<
            int dim, typename density_type, typename densityDot_type, typename result_type,
            typename testFunction_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &test_function,
            result_type &result
        ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &test_function: The value of the test function \f$ \psi \f$
             * \param &result: The net mass change per unit volume \f$ c \f$
             */

            computeBalanceOfMass<dim>(
                density, density_dot, density_gradient_begin, density_gradient_end,
                velocity_begin, velocity_end, velocity_gradient_begin, velocity_gradient_end,
                result
            );

            result *= test_function;

        }

        template<
            int dim, typename density_type, typename densityDot_type, typename result_type,
            typename dRdRho_type, typename dRdRhoDot_type, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class dRdGradRho_iter_out,
            class dRdV_iter_out, class dRdGradV_iter_out
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            result_type &result,
            dRdRho_type &dRdRho,                  dRdRhoDot_type &dRdRhoDot,
            dRdGradRho_iter_out dRdGradRho_begin, dRdGradRho_iter_out dRdGradRho_end,
            dRdV_iter_out dRdV_begin,             dRdV_iter_out dRdV_end,
            dRdGradV_iter_out dRdGradV_begin,     dRdGradV_iter_out dRdGradV_end
        ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &result: The net mass change per unit volume \f$ c \f$
             * \param &dRdRho: The derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dRdRhoDot: The derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dRdGradRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dRdGradRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dRdV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dRdV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dRdGradV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &dRdGradV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             */

            computeBalanceOfMass<dim>(
                density, density_dot,
                density_gradient_begin, density_gradient_end,
                velocity_begin, velocity_end,
                velocity_gradient_begin, velocity_gradient_end, result
            );

            dRdRho = 0;

            dRdRhoDot = 1;

            std::copy( velocity_begin, velocity_end, dRdGradRho_begin );

            std::copy( density_gradient_begin, density_gradient_end, dRdV_begin );

            std::fill( dRdGradV_begin, dRdGradV_end, 0 );

            for ( unsigned int i = 0; i < dim; ++i ){

                dRdRho += ( *( velocity_gradient_begin + dim * i + i ) );

                ( *( dRdGradV_begin + dim * i + i ) ) = density;

            }

        }

        template<
            int dim, typename density_type, typename densityDot_type, typename result_type,
            typename testFunction_type, typename interpolationFunction_type,
            typename dRdRho_type, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class interpolationFunctionGradient_iter,
            class dRdU_iter_out, class dRdUMesh_iter_out,
            typename dDensityDotdDensity_type, typename dUDotdU_type
        >
        void computeBalanceOfMass(
            const density_type &density,                                  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,           const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                          const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,         const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &test_function,                       const interpolationFunction_type &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const dDensityDotdDensity_type &dDensityDotdDensity,          const dUDotdU_type &dUDotdU,
            result_type &result,
            dRdRho_type &dRdRho,
            dRdU_iter_out dRdU_begin,             dRdU_iter_out dRdU_end,
            dRdUMesh_iter_out dRdUMesh_begin,     dRdUMesh_iter_out dRdUMesh_end
        ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * By assuming that the function is being used for a Galerkin integration method we can reduce the
             * dimensionality.
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             * 
             * We note that dRdUMesh also includes the term associated with the volumetric change ( \f$ c \phi_{,i} \f$ )
             * where \f$ \phi_{,i} \f$ is the interpolation function
             *
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &test_function: The test function \f$ \psi \f$
             * \param &interpolation_function: The interpolation function \f$ \phi \f$
             * \param &interpolation_function_gradient_begin: The starting iterator of the spatial gradient of the interpolation function \f$ \phi_{,i} \f$
             * \param &interpolation_function_gradient_end: The stopping iterator of the spatial gradient of the interpolation function \f$ \phi_{,i} \f$
             * \param &dDensityDotdDensity: The derivative of the time-derivative of the density w.r.t. the density (based on timestep and integration scheme)
             * \param &dUDotdU: The derivative of the time-derivative of the displacement w.r.t. the displacement (may not be mesh displacement)
             * \param &result: The net mass change per unit volume \f$ c \f$
             * \param &dRdRho: The derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dRdU_begin: The starting iterator of the derivative of the mass change rate w.r.t. the displacement (may not be mesh displacement)
             * \param &dRdU_end: The stopping iterator of the derivative of the mass change rate w.r.t. the displacement (may not be mesh displacement)
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the mass change rate w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the mass change rate w.r.t. the mesh displacement
             */

            density_type dRdRhoDot;

            std::array< density_type, dim > dRdGradRho;
            std::array< typename std::iterator_traits<dRdU_iter_out>::value_type, dim > dRdV;
            std::array< typename std::iterator_traits<dRdU_iter_out>::value_type, dim * dim > dRdGradV;

            computeBalanceOfMass<dim>(
                density, density_dot,
                density_gradient_begin,   density_gradient_end,
                velocity_begin,           velocity_end,
                velocity_gradient_begin,  velocity_gradient_end,
                result,
                dRdRho,                   dRdRhoDot,
                std::begin( dRdGradRho ), std::end( dRdGradRho ),
                std::begin( dRdV ),       std::end( dRdV ),
                std::begin( dRdGradV ),   std::end( dRdGradV )
            );

            // Set the mass change rate
            result *= test_function;

            // Assemble the derivatives w.r.t. the density
            dRdRho += dRdRhoDot * dDensityDotdDensity;
            dRdRho *= test_function * interpolation_function;

            dRdRho += test_function * std::inner_product( std::begin( dRdGradRho ), std::end( dRdGradRho ), interpolation_function_gradient_begin, dRdRho_type( ) );

            // Assemble the derivatives w.r.t. the mesh displacement
            std::fill( dRdU_begin, dRdU_end, 0 );
            std::fill( dRdUMesh_begin, dRdUMesh_end, 0 );

            for ( unsigned int a = 0; a < dim; ++a ){

                *( dRdU_begin + a ) += test_function * dRdV[ a ] * dUDotdU * interpolation_function;

                *( dRdUMesh_begin + a ) += result * ( *( interpolation_function_gradient_begin + a ) );

                for ( unsigned int i = 0; i < dim; ++i ){

                    *( dRdU_begin + a ) += test_function * dRdGradV[ dim * a + i ] * ( *( interpolation_function_gradient_begin + i ) ) * dUDotdU;

                    *( dRdUMesh_begin + a ) -= test_function * dRdGradRho[ i ] * ( *( density_gradient_begin + a ) ) * ( *( interpolation_function_gradient_begin + i ) );

                    for ( unsigned int j = 0; j < dim; ++j ){

                        *( dRdUMesh_begin + a ) -= test_function * dRdGradV[ dim * j + i ] * ( *( velocity_gradient_begin + dim * j + a ) ) * ( *( interpolation_function_gradient_begin + i ) );

                    }

                }

            }

        }

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class result_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            result_iter_out result_begin,                               result_iter_out result_end
        ){
            /*!
             * Compute the balance of mass for a multi-phase continuum returning the values of the mass-change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density_begin: The starting iterator of the density \f$ \rho \f$
             * \param &density_end: The stopping iterator of the density \f$ \rho \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &result_begin: The starting iterator of the net mass change per unit volume \f$ c \f$
             * \param &result_end: The stopping iterator of the net mass change per unit volume \f$ c \f$
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density dot arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( result_end - result_begin ), "The density and result arrays must be the same length" );

            unsigned int phase;

            for ( auto rho = density_begin; rho != density_end; ++rho ){

                phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfMass<dim>( *( density_begin + phase ), *( density_dot_begin + phase ),
                                           density_gradient_begin + dim * phase,        density_gradient_begin + dim * ( phase + 1 ),
                                           velocity_begin + dim * phase,                velocity_begin + dim * ( phase + 1 ),
                                           velocity_gradient_begin + dim * dim * phase, velocity_gradient_begin + dim * dim * ( phase + 1 ),
                                           *( result_begin + phase ) );

            }

        }

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            typename testFunction_type,
            class velocity_iter, class velocityGradient_iter, class result_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &test_function,
            result_iter_out result_begin,                               result_iter_out result_end
        ){
            /*!
             * Compute the balance of mass for a multi-phase continuum returning the values of the mass-change rate
             * 
             * Variational formulation with a test function
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density_begin: The starting iterator of the density \f$ \rho \f$
             * \param &density_end: The stopping iterator of the density \f$ \rho \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &test_function: The value of the test function \f$ \psi \f$
             * \param &result_begin: The starting iterator of the net mass change per unit volume \f$ c \f$
             * \param &result_end: The stopping iterator of the net mass change per unit volume \f$ c \f$
             */

            computeBalanceOfMass<dim>(
                density_begin, density_end, density_dot_begin, density_dot_end,
                density_gradient_begin, density_gradient_end,
                velocity_begin, velocity_end,
                velocity_gradient_begin, velocity_gradient_end,
                result_begin, result_end
            );

            std::transform(
                result_begin, result_end, result_begin,
                std::bind(
                    std::multiplies< typename std::iterator_traits<result_iter_out>::value_type >( ),
                    std::placeholders::_1,
                    test_function
                )
            );

        }

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class result_iter_out,
            class dRdRho_iter_out, class dRdRhoDot_iter_out, class dRdGradRho_iter_out,
            class dRdV_iter_out, class dRdGradV_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                    const density_iter &density_end,
            const densityDot_iter &density_dot_begin,             const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            result_iter_out result_begin,                         result_iter_out result_end,
            dRdRho_iter_out dRdRho_begin,                         dRdRho_iter_out dRdRho_end,
            dRdRhoDot_iter_out dRdRhoDot_begin,                   dRdRhoDot_iter_out dRdRhoDot_end,
            dRdGradRho_iter_out dRdGradRho_begin,                 dRdGradRho_iter_out dRdGradRho_end,
            dRdV_iter_out dRdV_begin,                             dRdV_iter_out dRdV_end,
            dRdGradV_iter_out dRdGradV_begin,                     dRdGradV_iter_out dRdGradV_end
        ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density_begin: The starting iterator of the value of the density \f$ \rho \f$
             * \param &density_end: The stopping iterator of the value of the density \f$ \rho \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &result_begin: The starting iterator of the net mass change per unit volume \f$ c \f$
             * \param &result_end: The stopping iterator of the net mass change per unit volume \f$ c \f$
             * \param &dRdRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dRdRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dRdRhoDot_begin: The starting iterator of the derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dRdRhoDot_end: The stopping iterator of the derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dRdGradRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dRdGradRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dRdV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dRdV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dRdGradV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &dRdGradV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density dot arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( result_end - result_begin ), "The density and result arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdRho_end - dRdRho_begin ), "The density and dRdRho arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdRhoDot_end - dRdRhoDot_begin ), "The density and dRdRhoDot arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdGradRho_end - dRdGradRho_begin ), "The density and dRdGradRho arrays must have consistent lengths" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdV_end - dRdV_begin ), "The density and dRdV arrays must have consistent lengths" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdGradV_end - dRdGradV_begin ), "The density and dRdGradV arrays must have consistent lengths" );

            unsigned int phase;

            for ( auto rho = density_begin; rho != density_end; ++rho ){

                phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfMass<dim>( *( density_begin + phase ),                *( density_dot_begin + phase ),
                                           density_gradient_begin + dim * phase,      density_gradient_begin + dim * ( phase + 1 ),
                                           velocity_begin + dim * phase,              velocity_begin + dim * ( phase + 1 ),
                                           velocity_gradient_begin + dim * dim * phase, velocity_gradient_begin + dim * dim * ( phase + 1 ),
                                           *( result_begin + phase ),       *( dRdRho_begin + phase ), *( dRdRhoDot_begin + phase ),
                                           dRdGradRho_begin + dim * phase,            dRdGradRho_begin + dim * ( phase + 1 ),
                                           dRdV_begin + dim * phase,                  dRdV_begin + dim * ( phase + 1 ),
                                           dRdGradV_begin + dim * dim * phase,         dRdGradV_begin + dim * dim * ( phase + 1 ) );

            }

        }

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            typename testFunction_type, typename interpolationFunction_type,
            class velocity_iter, class velocityGradient_iter, class interpolationFunctionGradient_iter,
            class result_iter_out,
            class dRdRho_iter_out, class dRdU_iter_out, class dRdUMesh_iter_out,
            class dDensityDotdDensity_iter, class dUDotdU_iter
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                            const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                     const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,           const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                          const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,         const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &test_function,                       const interpolationFunction_type &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const dDensityDotdDensity_iter &dDensityDotdDensity_begin,    const dDensityDotdDensity_iter &dDensityDotdDensity_end,      
            const dUDotdU_iter &dUDotdU_begin,                            const dUDotdU_iter &dUDotdU_end,
            result_iter_out result_begin,                                 result_iter_out result_end,
            dRdRho_iter_out dRdRho_begin,         dRdRho_iter_out dRdRho_end,
            dRdU_iter_out dRdU_begin,             dRdU_iter_out dRdU_end,
            dRdUMesh_iter_out dRdUMesh_begin,     dRdUMesh_iter_out dRdUMesh_end
        ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             *
             * \param &density_begin: The starting iterator of the value of the density \f$ \rho \f$
             * \param &density_end: The stopping iterator of the value of the density \f$ \rho \f$
             * \param &density_dot_begin: The starting iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_dot_end: The stopping iterator of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &test_function: The value of the test function \f$ \psi \f$
             * \param &interpolation_function: The value of the interpolation function \f$ \phi \f$
             * \param &interpolation_function_gradient_begin: The starting iterator of the spatial gradient of the interpolation function \f$ \phi_{,i} \f$
             * \param &interpolation_function_gradient_end: The stopping iterator of the spatial gradient of the interpolation function \f$ \phi_{,i} \f$
             * \param &dDensityDotdDensity_begin: The starting iterator of the derivative of the time rate of change of the density w.r.t. the density
             * \param &dDensityDotdDensity_end: The stopping iterator of the derivative of the time rate of change of the density w.r.t. the density
             * \param &dUDotdU_begin: The starting iterator of the derivative of the time rate of change of the phase DOF w.r.t. the phase DOF
             * \param &dUDotdU_end: The stopping iterator of the derivative of the time rate of change of the phase DOF w.r.t. the phase DOF
             * \param &result_begin: The starting iterator of the net mass change per unit volume \f$ c \f$
             * \param &result_end: The stopping iterator of the net mass change per unit volume \f$ c \f$
             * \param &dRdRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$.
             *     The derivatives are stored as the derivative of each phase w.r.t. the phase's density in that order.
             * \param &dRdRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             *     The derivatives are stored as the derivative of each phase w.r.t. the phase's density in that order.
             * \param &dRdU_begin: The starting iterator of the derivative of the mass change rate w.r.t. the phase displacement \f$ u_{i} \f$
             *     The derivatives are stored as the derivative of each phase w.r.t. the phase's displacement in that order
             * \param &dRdU_end: The stopping iterator of the derivative of the mass change rate w.r.t. the phase displacement \f$ u_{i} \f$
             *     The derivatives are stored as the derivative of each phase w.r.t. the phase's displacement in that order
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the mass change rate w.r.t. the mesh displacement \f$ u_{i} \f$
             *     The derivatives are stored as the derivative of each phase w.r.t. the mesh displacement in that order
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the mass change rate w.r.t. the mesh displacement \f$ u_{i} \f$
             *     The derivatives are stored as the derivative of each phase w.r.t. the mesh displacement in that order
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density dot arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( result_end - result_begin ), "The density and result arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dDensityDotdDensity_end - dDensityDotdDensity_begin ), "The density and the derivative of the density time derivative w.r.t. the density must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dUDotdU_end - dUDotdU_begin ), "The density and the derivative of the dof time derivative w.r.t. the dof must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdRho_end - dRdRho_begin ), "The density and dRdRho arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdU_end - dRdU_begin ), "The density and dRdU arrays must have consistent lengths" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dRdUMesh_end - dRdUMesh_begin ), "The density and dRdUMesh arrays must have consistent lengths" );

            unsigned int phase;

            for ( auto rho = density_begin; rho != density_end; ++rho ){

                phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfMass<dim>(
                    *( density_begin + phase ),                  *( density_dot_begin + phase ),
                    density_gradient_begin + dim * phase,        density_gradient_begin + dim * ( phase + 1 ),
                    velocity_begin + dim * phase,                velocity_begin + dim * ( phase + 1 ),
                    velocity_gradient_begin + dim * dim * phase, velocity_gradient_begin + dim * dim * ( phase + 1 ),
                    test_function,                               interpolation_function,
                    interpolation_function_gradient_begin,       interpolation_function_gradient_end,
                    *( dDensityDotdDensity_begin + phase ),      *( dUDotdU_begin + phase ),
                    *( result_begin + phase ),
                    *( dRdRho_begin + phase ),
                    dRdU_begin + dim * phase,                    dRdU_begin + dim * ( phase + 1 ),
                    dRdUMesh_begin + dim * phase,                dRdUMesh_begin + dim * ( phase + 1 )
                );

            }

        }

        template<
            int dim, int matresponse_dim, typename density_type, typename densityDot_type, typename result_type,
            typename testFunction_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter,
            class material_response_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin, const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,  const velocityGradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const testFunction_type &test_function,
            result_type &result
        ){

            /*!
             * A balance of mass function for a general material response problem where the change in mass may be
             * a function of many different variables.
             * 
             * The mass rate of change value \f$ c \f$ is assumed to be located after the stress measure
             * which has dimensions of dim * dim and the model-calculated internal energy which has a dimension of 1
             * 
             * The material response vector is assumed to be organized as
             * 
             * \f$ \sigma, e, c, b, \pi, q, r, \Pi \f$
             *
             * where \f$ \sigma \f$ is the Cauchy stress, \f$ e \f$ is the predicted internal energy per unit mass,
             * \f$ c \f$ is the mass change rate per unit volume, \f$ b \f$ is the body force vector, \f$ \pi \f$ is
             * the net inter-phase force, \f$ q \f$ is the heat flux vector, \f$ r \f$ is the volumetric heat generation
             * per unit mass, and \f$ \Pi \f$ is the net interphase heat transfer. We note that additional outputs can be
             * added to the material response vector, but they must be done after the above quantities.
             * 
             * 
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &material_response_begin: The starting iterator of the material response
             * \param &material_response_end: The stopping iterator of the material response
             * \param &test_function: The value of the test function \f$ \psi \f$
             * \param &result: The net mass change per unit volume \f$ c \f$
             */

            // Compute the non-mass change parts of the balance of mass
            computeBalanceOfMass<dim>(
                density, density_dot, density_gradient_begin, density_gradient_end,
                velocity_begin, velocity_end, velocity_gradient_begin, velocity_gradient_end,
                test_function,
                result
            );

            // Add in the contributions from the change in mass
            result -= test_function * ( *( material_response_begin + matresponse_dim * matresponse_dim + 1 ) );

        }

        template<
            int dim, int matresponse_dim, int matresponse_values, class density_iter, class densityDot_iter, class result_iter,
            typename testFunction_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter,
            class material_response_iter
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                     const density_iter &density_end,
            const densityDot_iter &density_dot_begin,              const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,    const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,  const velocityGradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const testFunction_type &test_function,
            result_iter result_begin, result_iter result_end
        ){

            /*!
             * A balance of mass function for a general material response problem where the change in mass may be
             * a function of many different variables. Evaluates for all phases.
             * 
             * The mass rate of change value \f$ c \f$ is assumed to be located after the stress measure
             * which has dimensions of dim * dim and the model-calculated internal energy which has a dimension of 1
             * 
             * The material response vector is assumed to be organized as
             * 
             * \f$ \sigma, e, c, b, \pi, q, r, \Pi \f$
             *
             * where \f$ \sigma \f$ is the Cauchy stress, \f$ e \f$ is the predicted internal energy per unit mass,
             * \f$ c \f$ is the mass change rate per unit volume, \f$ b \f$ is the body force vector, \f$ \pi \f$ is
             * the net inter-phase force, \f$ q \f$ is the heat flux vector, \f$ r \f$ is the volumetric heat generation
             * per unit mass, and \f$ \Pi \f$ is the net interphase heat transfer. We note that additional outputs can be
             * added to the material response vector, but they must be done after the above quantities.
             * 
             * 
             * \param &density_begin: The starting iterator of the value of the density \f$ \rho \f$
             * \param &density_end: The stopping iterator of the value of the density \f$ \rho \f$
             * \param &density_dot_begin: The starting iterator of the value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_dot_end: The stopping iterator of the value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &material_response_begin: The starting iterator of the material response the spatial dimension is matresponse_dim
             *     and the number of values defined for each phase is matresponse_values
             * \param &material_response_end: The stopping iterator of the material response the spatial dimension is matresponse_dim
             *     and the number of values defined for each phase is matresponse_values
             * \param &test_function: The value of the test function \f$ \psi \f$
             * \param &result_begin: The starting iterator of the net mass change per unit volume \f$ c \f$
             * \param &result_end: The stopping iterator of the net mass change per unit volume \f$ c \f$
             */

            TARDIGRADE_ERROR_TOOLS_EVAL(
                const unsigned int num_phases = ( unsigned int )( density_end - density_begin );
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                num_phases == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density dot arrays must have the same length"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * num_phases == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient arrays must have consistent lengths"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * num_phases == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity arrays must have consistent lengths"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * dim * num_phases == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient arrays must have consistent lengths"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                matresponse_values * num_phases == ( unsigned int )( material_response_end - material_response_begin ), "The density and material response arrays must have consistent lengths"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                num_phases == ( unsigned int )( result_end - result_begin ), "The density and result arrays must have the same length"
            )

            for ( auto v = std::pair< unsigned int, density_iter >( 0, density_begin ); v.second != density_end; ++v.first, ++v.second ){

                computeBalanceOfMass<dim, matresponse_dim>(
                    *v.second, *( density_dot_begin + v.first ), density_gradient_begin + dim * v.first, density_gradient_begin + dim * ( v.first + 1 ),
                    velocity_begin + dim * v.first, velocity_begin + dim * ( v.first + 1 ), velocity_gradient_begin + dim * dim * v.first, velocity_gradient_begin + dim * dim * ( v.first + 1 ),
                    material_response_begin + matresponse_values * v.first, material_response_begin + matresponse_values * ( v.first + 1 ),
                    test_function,
                    *( result_begin + v.first )
                );

            }

        }

    }

}
