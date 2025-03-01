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
            class velocity_iter, class velocityGradient_iter, class dRdGradRho_iter,
            class dRdV_iter, class dRdGradV_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            result_type &result,
            dRdRho_type &dRdRho,                  dRdRhoDot_type &dRdRhoDot,
            dRdGradRho_iter dRdGradRho_begin, dRdGradRho_iter dRdGradRho_end,
            dRdV_iter dRdV_begin,             dRdV_iter dRdV_end,
            dRdGradV_iter dRdGradV_begin,     dRdGradV_iter dRdGradV_end
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
            class dRdU_iter, class dRdUMesh_iter,
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
            dRdU_iter dRdU_begin,             dRdU_iter dRdU_end,
            dRdUMesh_iter dRdUMesh_begin,     dRdUMesh_iter dRdUMesh_end
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
            std::array< typename std::iterator_traits<dRdU_iter>::value_type, dim > dRdV;
            std::array< typename std::iterator_traits<dRdU_iter>::value_type, dim * dim > dRdGradV;

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
            class velocity_iter, class velocityGradient_iter, class result_iter
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            result_iter result_begin,                               result_iter result_end
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
            class velocity_iter, class velocityGradient_iter, class result_iter
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &test_function,
            result_iter result_begin,                               result_iter result_end
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
                    std::multiplies< typename std::iterator_traits<result_iter>::value_type >( ),
                    std::placeholders::_1,
                    test_function
                )
            );

        }

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class result_iter,
            class dRdRho_iter, class dRdRhoDot_iter, class dRdGradRho_iter,
            class dRdV_iter, class dRdGradV_iter
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                    const density_iter &density_end,
            const densityDot_iter &density_dot_begin,             const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            result_iter result_begin,                         result_iter result_end,
            dRdRho_iter dRdRho_begin,                         dRdRho_iter dRdRho_end,
            dRdRhoDot_iter dRdRhoDot_begin,                   dRdRhoDot_iter dRdRhoDot_end,
            dRdGradRho_iter dRdGradRho_begin,                 dRdGradRho_iter dRdGradRho_end,
            dRdV_iter dRdV_begin,                             dRdV_iter dRdV_end,
            dRdGradV_iter dRdGradV_begin,                     dRdGradV_iter dRdGradV_end
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
            class result_iter,
            class dRdRho_iter, class dRdU_iter, class dRdUMesh_iter,
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
            result_iter result_begin,                                 result_iter result_end,
            dRdRho_iter dRdRho_begin,         dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin,             dRdU_iter dRdU_end,
            dRdUMesh_iter dRdUMesh_begin,     dRdUMesh_iter dRdUMesh_end
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
            int dim, int mass_change_index, typename density_type, typename densityDot_type, typename result_type,
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
            result -= test_function * ( *( material_response_begin + mass_change_index ) );

        }

        template<
            int dim, int mass_change_index, class density_iter, class densityDot_iter, class result_iter,
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
             * mass_change_index is the index of the material response vector that represents the mass change rate
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
             * \param &material_response_begin: The starting iterator of the material response the spatial dimension is material_response_dim
             *     and the number of values defined for each phase is material_response_size
             * \param &material_response_end: The stopping iterator of the material response the spatial dimension is material_response_dim
             *     and the number of values defined for each phase is material_response_size
             * \param &test_function: The value of the test function \f$ \psi \f$
             * \param &result_begin: The starting iterator of the net mass change per unit volume \f$ c \f$
             * \param &result_end: The stopping iterator of the net mass change per unit volume \f$ c \f$
             */

            const unsigned int num_phases = ( unsigned int )( density_end - density_begin );
            const unsigned int material_response_size = ( unsigned int )( material_response_end - material_response_begin ) / num_phases;

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
                material_response_size * num_phases == ( unsigned int )( material_response_end - material_response_begin ), "The density and material response arrays must have consistent lengths"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                num_phases == ( unsigned int )( result_end - result_begin ), "The density and result arrays must have the same length"
            )

            for ( auto v = std::pair< unsigned int, density_iter >( 0, density_begin ); v.second != density_end; ++v.first, ++v.second ){

                computeBalanceOfMass<dim, mass_change_index>(
                    *v.second, *( density_dot_begin + v.first ), density_gradient_begin + dim * v.first, density_gradient_begin + dim * ( v.first + 1 ),
                    velocity_begin + dim * v.first, velocity_begin + dim * ( v.first + 1 ), velocity_gradient_begin + dim * dim * v.first, velocity_gradient_begin + dim * dim * ( v.first + 1 ),
                    material_response_begin + material_response_size * v.first, material_response_begin + material_response_size * ( v.first + 1 ),
                    test_function,
                    *( result_begin + v.first )
                );

            }

        }

        template<
            int dim, int mass_change_index, int material_response_dim, int material_response_num_dof,
            typename density_type, typename densityDot_type, typename result_type,
            typename testFunction_type, typename interpolationFunction_type,
            class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter,
            class material_response_iter, class material_response_jacobian_iter,
            class interpolationFunctionGradient_iter,
            class full_material_response_dof_gradient_iter,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter, class dRdTheta_iter, class dRdE_iter, class dRdZ_iter,
            class dRdUMesh_iter,
            typename dDensityDotdDensity_type, typename dUDotdU_type,
            int density_index,
            int displacement_index,
            int velocity_index,
            int temperature_index,
            int internal_energy_index,
            int additional_dof_index
        >
        void computeBalanceOfMass(
            const density_type &density,                                  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,           const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                          const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,         const velocityGradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin,        const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const testFunction_type &test_function,                       const interpolationFunction_type &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dDensityDotdDensity_type &dDensityDotdDensity,          const dUDotdU_type &dUDotdU,
            const unsigned int phase,
            result_type &result,
            dRdRho_iter dRdRho_begin,         dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin,             dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin,             dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin,     dRdTheta_iter dRdTheta_end,
            dRdE_iter dRdE_begin,             dRdE_iter dRdE_end,
            dRdZ_iter dRdZ_begin,             dRdZ_iter dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin,     dRdUMesh_iter dRdUMesh_end
        ){

            /*!
             * A balance of mass function for a general material response problem where the change in mass may be
             * a function of many different variables.
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
             * We assume that the material_response_jacobian is organized such that there are nphases times the number of
             * phase degrees of freedom plus the additional dof followed by the gradients of these degrees of freedom
             * for the columns of the Jacobian. If the degrees of freedom are \f$ rho, w, \f$ and \f$ v \f$ then
             * the columns of the jacobian would be 
             * 
             * 
             * \f$ \rho, w, v, \nabla \rho, \nabla w, \nabla v \f$
             * 
             * If the degrees of freedom are vector quantities, then they, and their gradients, must be stored in
             * row-major form.
             * 
             * In the case of a multi-phase problem we expect the following structure:
             * 
             * \f$ \rho^1, \rho^2, \ldots w^1, w^2, \ldots, v^1, v^2, \ldots, \nabla \rho^1, \nabla \rho^2, ... \f$
             *
             * where the superscript indicates the phase. In this function we assume a default order of
             * 
             * \f$ \rho, w, v, \theta, e, z \f$
             * 
             * Where \f$ v \f$ is the velocity which may be equal to \f$ u \f$ or it's time derivative.
             * 
             * mass_change_index is the index of the material response vector that represents the mass change rate
             * material_response_dim is the spatial dimension of the material-response Jacobian
             * material_response_num_dof are the number of degrees of freedom in the material-response Jacobian for each phase
             * We note that there are 3 + 2 * material_response_dim expected degrees of freedom.
             * 
             * \param &density: The value of the density \f$ \rho \f$
             * \param &density_dot: The value of the partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &velocity_begin: The starting iterator of the velocity \f$ v_i \f$
             * \param &velocity_end: The stopping iterator of the velocity \f$ v_i \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &test_function: The test function \f$ \psi \f$
             * \param &interpolation_function: The interpolation function \f$ \phi \f$
             * \param &interpolation_function_gradient_begin: The starting iterator of the spatial gradient of the interpolation function \f$ \phi_{,i} \f$
             * \param &interpolation_function_gradient_end: The stopping iterator of the spatial gradient of the interpolation function \f$ \phi_{,i} \f$
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of all of the degrees of freedom used by the material response
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of all of the degrees of freedom used by the material response
             * \param &dDensityDotdDensity: The derivative of the time-derivative of the density w.r.t. the density (based on timestep and integration scheme)
             * \param &dUDotdU: The derivative of the time-derivative of the displacement w.r.t. the displacement (may not be mesh displacement)
             * \param &phase: The phase the balance equation applies to
             * \param &result: The net mass change per unit volume \f$ c \f$
             * \param &dRdRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dRdRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dRdU_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial dof (may not be displacement)
             * \param &dRdU_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial dof (may not be displacement)
             * \param &dRdW_begin: The starting iterator of the derivative of the mass change rate w.r.t. the displacement (may not be mesh displacement)
             * \param &dRdW_end: The stopping iterator of the derivative of the mass change rate w.r.t. the displacement (may not be mesh displacement)
             * \param &dRdTheta_begin: The starting iterator of the derivative of the mass change rate w.r.t. the temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the mass change rate w.r.t. the temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the mass change rate w.r.t. the internal energy
             * \param &dRdE_end: The stopping iterator of the derivative of the mass change rate w.r.t. the internal energy
             * \param &dRdZ_begin: The starting iterator of the derivative of the mass change rate w.r.t. the additional dof
             * \param &dRdZ_end: The stopping iterator of the derivative of the mass change rate w.r.t. the additional dof
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the mass change rate w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the mass change rate w.r.t. the mesh displacement
             */

            using dRdRho_type = typename std::iterator_traits<dRdRho_iter>::value_type;

            using dRdU_type = typename std::iterator_traits<dRdU_iter>::value_type;

            dRdRho_type phase_dRdRho;

            std::array< dRdU_type, dim > phase_dRdU;

            // Compute the non-mass change parts of the balance of mass
            computeBalanceOfMass<dim>(
                density, density_dot, density_gradient_begin, density_gradient_end,
                velocity_begin, velocity_end, velocity_gradient_begin, velocity_gradient_end,
                test_function, interpolation_function,
                interpolation_function_gradient_begin, interpolation_function_gradient_end,
                dDensityDotdDensity, dUDotdU,
                result,
                phase_dRdRho,
                std::begin( phase_dRdU ),   std::end( phase_dRdU ),
                dRdUMesh_begin,             dRdUMesh_end
            );

            // Add in the contributions from the change in mass
            result -= test_function * ( *( material_response_begin + mass_change_index ) );

            // Zero out the Jacobians
            std::fill( dRdRho_begin,   dRdRho_end,   0 );
            std::fill( dRdU_begin,     dRdU_end,     0 );
            std::fill( dRdW_begin,     dRdW_end,     0 );
            std::fill( dRdTheta_begin, dRdTheta_end, 0 );
            std::fill( dRdE_begin,     dRdE_end,     0 );
            std::fill( dRdZ_begin,     dRdZ_end,     0 );

            // Set the number of phases
            const unsigned int nphases = ( unsigned int )( dRdRho_end - dRdRho_begin );
            constexpr unsigned int num_phase_dof = 3 + 2 * material_response_dim;
            const unsigned int num_additional_dof = ( material_response_num_dof - num_phase_dof );

            // Add the material response contributions to the density Jacobian
            for ( auto p = std::pair< unsigned int, dRdRho_iter >( 0, dRdRho_begin ); p.second != dRdRho_end; ++p.first, ++p.second ){

                // DOF value contributions
                *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + nphases * density_index + p.first ) ) * interpolation_function;

                // DOF spatial gradient contributions
                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * density_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            *( dRdRho_begin + phase ) += phase_dRdRho;

            // Add the material response contributions to the U vector
            for ( auto p = std::pair< unsigned int, dRdU_iter >( 0, dRdU_begin ); p.second != dRdU_end; ++p.first, ++p.second ){

                // DOF value contributions
                *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + nphases * velocity_index + p.first ) ) * dUDotdU * interpolation_function;

                // DOF spatial gradient contributions
                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * velocity_index + p.first ) + a ) ) * dUDotdU * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            for ( unsigned int a = 0; a < dim; ++a ){

                *( dRdU_begin + dim * phase + a ) += phase_dRdU[ a ];

            }

            // Add the material response contributions to the displacement Jacobian
            for ( auto p = std::pair< unsigned int, dRdW_iter >( 0, dRdW_begin ); p.second != dRdW_end; ++p.first, ++p.second ){

                // DOF value contributions
                *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + nphases * displacement_index + p.first ) ) * interpolation_function;

                // DOF spatial gradient contributions
                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * displacement_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // Add the material response contributions to the temperature
            for ( auto p = std::pair< unsigned int, dRdTheta_iter >( 0, dRdTheta_begin ); p.second != dRdTheta_end; ++p.first, ++p.second ){

                // DOF value contributions
                *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + nphases * temperature_index + p.first ) ) * interpolation_function;

                // DOF spatial gradient contributions
                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * temperature_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // Add the material response contributions to the internal energy
            for ( auto p = std::pair< unsigned int, dRdE_iter >( 0, dRdE_begin ); p.second != dRdE_end; ++p.first, ++p.second ){

                // DOF value contributions
                *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + nphases * internal_energy_index + p.first ) ) * interpolation_function;

                // DOF spatial gradient contributions
                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * internal_energy_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // Add the material response contributions to the additional degrees of freedom
            for ( auto p = std::pair< unsigned int, dRdZ_iter >( 0, dRdZ_begin ); p.second != dRdZ_end; ++p.first, ++p.second ){

                // DOF value contributions
                *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + nphases * additional_dof_index + p.first ) ) * interpolation_function;

                // DOF spatial gradient contributions
                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second -= test_function * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * additional_dof_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // Add the material response contributions to the mesh displacement Jacobian
            for ( unsigned int i = 0; i < ( nphases * num_phase_dof + num_additional_dof ); ++i ){

                for ( unsigned int j = 0; j < material_response_dim; ++j ){

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *( dRdUMesh_begin + a ) += 
                            test_function
                            * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * mass_change_index + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * i + j ) )
                            * ( *( full_material_response_dof_gradient_begin + material_response_dim * i + a ) )
                            * ( *( interpolation_function_gradient_begin + j ) );

                    }

                }

            }

            for ( unsigned int a = 0; a < dim; ++a ){

                *( dRdUMesh_begin + a ) -= test_function * ( *( material_response_begin + mass_change_index ) ) * ( *( interpolation_function_gradient_begin + a ) );

            }

        }

        template<
            int dim, int mass_change_index, int material_response_dim, int material_response_num_dof,
            class density_iter, class densityDot_iter, class result_iter,
            typename testFunction_type, typename interpolationFunction_type,
            class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter,
            class material_response_iter, class material_response_jacobian_iter,
            class interpolationFunctionGradient_iter,
            class full_material_response_dof_gradient_iter,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter, class dRdTheta_iter, class dRdE_iter, class dRdZ_iter,
            class dRdUMesh_iter,
            typename dDensityDotdDensity_type, typename dUDotdU_type,
            int density_index,
            int displacement_index,
            int velocity_index,
            int temperature_index,
            int internal_energy_index,
            int additional_dof_index
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                            const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                     const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,           const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                          const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,         const velocityGradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin,        const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const testFunction_type &test_function,                       const interpolationFunction_type &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dDensityDotdDensity_type &dDensityDotdDensity,          const dUDotdU_type &dUDotdU,
            result_iter result_begin,         result_iter result_end,
            dRdRho_iter dRdRho_begin,         dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin,             dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin,             dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin,     dRdTheta_iter dRdTheta_end,
            dRdE_iter dRdE_begin,             dRdE_iter dRdE_end,
            dRdZ_iter dRdZ_begin,             dRdZ_iter dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin,     dRdUMesh_iter dRdUMesh_end
        ){

            /*!
             * A balance of mass function for a general material response problem where the change in mass may be
             * a function of many different variables. Evaluates for all phases.
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
             * We assume that the material_response_jacobian is organized such that there are ( nphases * num_phase_dof + num_additional_dof )
             * degrees of freedom followed by the gradients of these degrees of freedom for the columns of the Jacobian.
             * If the degrees of freedom are \f$ rho, w, \f$ and \f$ v \f$ then the columns of the jacobian would be 
             * 
             * 
             * \f$ \rho, w, v, \nabla \rho, \nabla w, \nabla v \f$
             * 
             * If the degrees of freedom are vector quantities, then they, and their gradients, must be stored in
             * row-major form.
             * 
             * In the case of a multi-phase problem we expect the following structure:
             * 
             * \f$ \rho^1, \rho^2, \ldots w^1, w^2, \ldots, v^1, v^2, \ldots, \nabla \rho^1, \nabla \rho^2, ... \f$
             *
             * where the superscript indicates the phase. In this function we assume a default order of
             * 
             * \f$ \rho, w, v, \theta, e, z \f$
             * 
             * Where \f$ v \f$ is the velocity which may be equal to \f$ u \f$ or it's time derivative.
             * 
             * mass_change_index is the index of the material response vector that represents the mass change rate
             * material_response_dim is the spatial dimension of the material-response Jacobian
             * ( nphases * num_phase_dof + num_additional_dof ) are the number of degrees of freedom in the material-response Jacobian
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
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &test_function: The test function \f$ \psi \f$
             * \param &interpolation_function: The interpolation function \f$ \phi \f$
             * \param &interpolation_function_gradient_begin: The starting iterator of the spatial gradient of the interpolation function \f$ \phi_{,i} \f$
             * \param &interpolation_function_gradient_end: The stopping iterator of the spatial gradient of the interpolation function \f$ \phi_{,i} \f$
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of all of the degrees of freedom used by the material response
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of all of the degrees of freedom used by the material response
             * \param &dDensityDotdDensity: The derivative of the time-derivative of the density w.r.t. the density (based on timestep and integration scheme)
             * \param &dUDotdU: The derivative of the time-derivative of the displacement w.r.t. the displacement (may not be mesh displacement)
             * \param &result_begin: The starting iterator of the net mass change per unit volume \f$ c \f$
             * \param &result_end: The stopping iterator of the net mass change per unit volume \f$ c \f$
             * \param &dRdRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dRdRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dRdU_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial dof (may not be displacement)
             * \param &dRdU_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial dof (may not be displacement)
             * \param &dRdW_begin: The starting iterator of the derivative of the mass change rate w.r.t. the displacement (may not be mesh displacement)
             * \param &dRdW_end: The stopping iterator of the derivative of the mass change rate w.r.t. the displacement (may not be mesh displacement)
             * \param &dRdTheta_begin: The starting iterator of the derivative of the mass change rate w.r.t. the temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the mass change rate w.r.t. the temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the mass change rate w.r.t. the internal energy
             * \param &dRdE_end: The stopping iterator of the derivative of the mass change rate w.r.t. the internal energy
             * \param &dRdZ_begin: The starting iterator of the derivative of the mass change rate w.r.t. the additional dof
             * \param &dRdZ_end: The stopping iterator of the derivative of the mass change rate w.r.t. the additional dof
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the mass change rate w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the mass change rate w.r.t. the mesh displacement
             */

            using density_type     = typename std::iterator_traits<density_iter>::value_type;

            using density_dot_type = typename std::iterator_traits<densityDot_iter>::value_type;

            using result_type      = typename std::iterator_traits<result_iter>::value_type;

            const unsigned int nphases = ( unsigned int )( density_end - density_begin );

            const unsigned int material_response_size = ( unsigned int )( material_response_end - material_response_begin ) / nphases;

            constexpr unsigned int num_phase_dof = 3 + 2 * material_response_dim;

            constexpr unsigned int num_additional_dof = material_response_num_dof - num_phase_dof;

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( density_dot_end - density_dot_begin ), "The length of density dot and density must be the same"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The length of the density gradient and the density must be consistent"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim == ( unsigned int )( velocity_end - velocity_begin ), "The length of the velocity and the density must be consistent"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim * dim == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The length of the velocity gradient and the density must be consistent"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                mass_change_index < material_response_size, "The material response vector must be larger than the mass-change index times the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size * ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) == ( unsigned int )( material_response_jacobian_end - material_response_jacobian_begin ),
                "The material response jacobian but have a consistent size with the material response vector and the material_response_num_dof"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim == ( unsigned int )( interpolation_function_gradient_end - interpolation_function_gradient_begin ),
                "The interpolation function gradient must have a size of dim"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( result_end - result_begin ),
                "The result vector must be the same size as the density vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * 1 == ( unsigned int )( dRdRho_end - dRdRho_begin ),
                "dRdRho must have a consistent size with the density vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * material_response_dim == ( unsigned int )( dRdU_end - dRdU_begin ),
                "dRdU must have a consistent size with the density vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * material_response_dim == ( unsigned int )( dRdW_end - dRdW_begin ),
                "dRdW must have a consistent size with the density vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * 1 == ( unsigned int )( dRdTheta_end - dRdTheta_begin ),
                "dRdTheta must have a consistent size with the density vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * 1 == ( unsigned int )( dRdE_end - dRdE_begin ),
                "dRdE must have a consistent size with the density vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * num_additional_dof == ( unsigned int )( dRdZ_end - dRdZ_begin ),
                "dRdZ must have a consistent size with the density vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim == ( unsigned int )( dRdUMesh_end - dRdUMesh_begin ),
                "dRdUMesh must have a consistent size with the density vector"
            )

            for(
                auto v = std::pair< unsigned int, density_iter >( 0, density_begin );
                v.second != density_end;
                ++v.first, ++v.second
            ){

                computeBalanceOfMass
                <
                    dim, mass_change_index, material_response_dim, material_response_num_dof,
                    density_type, density_dot_type, result_type,
                    testFunction_type, interpolationFunction_type,
                    densityGradient_iter,
                    velocity_iter, velocityGradient_iter,
                    material_response_iter, material_response_jacobian_iter,
                    interpolationFunctionGradient_iter,
                    full_material_response_dof_gradient_iter,
                    dRdRho_iter, dRdU_iter, dRdW_iter, dRdTheta_iter, dRdE_iter, dRdZ_iter,
                    dRdUMesh_iter,
                    dDensityDotdDensity_type, dUDotdU_type,
                    density_index,
                    displacement_index,
                    velocity_index,
                    temperature_index,
                    internal_energy_index,
                    additional_dof_index
                >
                (
                    *( density_begin + v.first ),                                 *( density_dot_begin + v.first ),
                    density_gradient_begin  + dim * v.first,                      density_gradient_begin  + dim * ( v.first + 1 ),
                    velocity_begin          + dim * v.first,                      velocity_begin          + dim * ( v.first + 1 ),
                    velocity_gradient_begin + dim * dim * v.first,                velocity_gradient_begin + dim * dim * ( v.first + 1 ),
                    material_response_begin + material_response_size * v.first,   material_response_begin + material_response_size * ( v.first + 1 ),
                    material_response_jacobian_begin + material_response_size * ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * v.first,
                    material_response_jacobian_begin + material_response_size * ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( v.first + 1 ),
                    test_function, interpolation_function,
                    interpolation_function_gradient_begin, interpolation_function_gradient_end,
                    full_material_response_dof_gradient_begin,
                    full_material_response_dof_gradient_end,
                    dDensityDotdDensity, dUDotdU,
                    v.first,
                    *( result_begin + v.first ),
                    dRdRho_begin   + nphases *   1 * v.first, dRdRho_begin   + nphases * 1   * ( v.first + 1 ),
                    dRdU_begin     + nphases * dim * v.first, dRdU_begin     + nphases * dim * ( v.first + 1 ),
                    dRdW_begin     + nphases * dim * v.first, dRdW_begin     + nphases * dim * ( v.first + 1 ),
                    dRdTheta_begin + nphases *   1 * v.first, dRdTheta_begin + nphases *   1 * ( v.first + 1 ),
                    dRdE_begin     + nphases *   1 * v.first, dRdE_begin     + nphases *   1 * ( v.first + 1 ),
                    dRdZ_begin     + num_additional_dof * v.first,
                    dRdZ_begin     + num_additional_dof * ( v.first + 1 ),
                    dRdUMesh_begin +       1 * dim * v.first, dRdUMesh_begin +       1 * dim * ( v.first + 1 )
                );

            }

        }

    }

}
