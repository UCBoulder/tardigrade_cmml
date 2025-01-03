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
            const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &mass_change_rate
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
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             */

            computeBalanceOfMass<global_dim>(
                density, density_dot,
                std::begin( density_gradient ),  std::end( density_gradient ),
                std::begin( velocity ),          std::end( velocity ),
                std::begin( velocity_gradient ), std::end( velocity_gradient ),
                mass_change_rate
            );

        }

        void computeBalanceOfMass(
            const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
            const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &mass_change_rate, 
            floatType   &dCdRho, floatType         &dCdRhoDot, floatVector &dCdGradRho,
            floatVector &dCdV,   secondOrderTensor &dCdGradV
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
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             * \param &dCdRho: The derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdRhoDot: The derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dCdGradRho: The derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dCdV: The derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dCdGradV: The derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             */

            computeBalanceOfMass<global_dim>(
                density,  density_dot,
                std::begin( density_gradient ),  std::end( density_gradient ),
                std::begin( velocity ),          std::end( velocity ),
                std::begin( velocity_gradient ), std::end( velocity_gradient ),
                mass_change_rate, dCdRho, dCdRhoDot,
                std::begin( dCdGradRho ), std::end( dCdGradRho ),
                std::begin( dCdV ),       std::end( dCdV ),
                std::begin( dCdGradV ),   std::end( dCdGradV )
            );
        }

        template<
            int dim, typename density_type, typename densityDot_type, typename c_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            c_type &mass_change_rate
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
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( ( unsigned int )( density_gradient_end - density_gradient_begin ) == ( unsigned int )( velocity_end - velocity_begin ), "The density gradient and the velocity must have the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ) == dim * dim, "The velocity gradient has a size of " + std::to_string( ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ) ) + " and must have a size of " + std::to_string( dim * dim ) );

            mass_change_rate = std::inner_product( density_gradient_begin, density_gradient_end,
                                                   velocity_begin, density_dot );

            for ( unsigned int i = 0; i < dim; ++i ){

                mass_change_rate += density * ( *( velocity_gradient_begin + dim * i + i ) );

            }

        }

        template<
            int dim, typename density_type, typename densityDot_type, typename c_type,
            typename testFunction_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &test_function,
            c_type &mass_change_rate
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
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             */

            computeBalanceOfMass<dim>(
                density, density_dot, density_gradient_begin, density_gradient_end,
                velocity_begin, velocity_end, velocity_gradient_begin, velocity_gradient_end,
                mass_change_rate
            );

            mass_change_rate *= test_function;

        }

        template<
            int dim, typename density_type, typename densityDot_type, typename c_type,
            typename dCdRho_type, typename dCdRhoDot_type, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class dCdGradRho_iter_out,
            class dCdV_iter_out, class dCdGradV_iter_out
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            c_type &mass_change_rate,
            dCdRho_type &dCdRho,                  dCdRhoDot_type &dCdRhoDot,
            dCdGradRho_iter_out dCdGradRho_begin, dCdGradRho_iter_out dCdGradRho_end,
            dCdV_iter_out dCdV_begin,             dCdV_iter_out dCdV_end,
            dCdGradV_iter_out dCdGradV_begin,     dCdGradV_iter_out dCdGradV_end
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
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             * \param &dCdRho: The derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdRhoDot: The derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dCdGradRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dCdGradRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dCdV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dCdV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dCdGradV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &dCdGradV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             */

            computeBalanceOfMass<dim>(
                density, density_dot,
                density_gradient_begin, density_gradient_end,
                velocity_begin, velocity_end,
                velocity_gradient_begin, velocity_gradient_end, mass_change_rate
            );

            dCdRho = 0;

            dCdRhoDot = 1;

            std::copy( velocity_begin, velocity_end, dCdGradRho_begin );

            std::copy( density_gradient_begin, density_gradient_end, dCdV_begin );

            std::fill( dCdGradV_begin, dCdGradV_end, 0 );

            for ( unsigned int i = 0; i < dim; ++i ){

                dCdRho += ( *( velocity_gradient_begin + dim * i + i ) );

                ( *( dCdGradV_begin + dim * i + i ) ) = density;

            }

        }

        template<
            int dim, typename density_type, typename densityDot_type, typename c_type,
            typename testFunction_type, typename interpolationFunction_type,
            typename dCdRho_type, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class interpolationFunctionGradient_iter,
            class dCdU_iter_out, class dCdUMesh_iter_out,
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
            c_type &mass_change_rate,
            dCdRho_type &dCdRho,
            dCdU_iter_out dCdU_begin,             dCdU_iter_out dCdU_end,
            dCdUMesh_iter_out dCdUMesh_begin,     dCdUMesh_iter_out dCdUMesh_end
        ){
            /*!
             * Compute the value of the balance of mass returning the value of the mass change rate
             * 
             * By assuming that the function is being used for a Galerkin integration method we can reduce the
             * dimensionality.
             * 
             * \f$ \frac{\partial \rho}{\partial t} + \left( \rho v_i \right)_{,i} = c \f$
             * 
             * We note that dCdUMesh also includes the term associated with the volumetric change ( \f$ c \phi_{,i} \f$ )
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
             * \param &mass_change_rate: The rate of change of the mass \f$ c \f$
             * \param &dCdRho: The derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdU_begin: The starting iterator of the derivative of the mass change rate w.r.t. the displacement (may not be mesh displacement)
             * \param &dCdU_end: The stopping iterator of the derivative of the mass change rate w.r.t. the displacement (may not be mesh displacement)
             * \param &dCdUMesh_begin: The starting iterator of the derivative of the mass change rate w.r.t. the mesh displacement
             * \param &dCdUMesh_end: The stopping iterator of the derivative of the mass change rate w.r.t. the mesh displacement
             */

            density_type dCdRhoDot;

            std::array< density_type, dim > dCdGradRho;
            std::array< typename std::iterator_traits<dCdU_iter_out>::value_type, dim > dCdV;
            std::array< typename std::iterator_traits<dCdU_iter_out>::value_type, dim * dim > dCdGradV;

            computeBalanceOfMass<dim>(
                density, density_dot,
                density_gradient_begin,   density_gradient_end,
                velocity_begin,           velocity_end,
                velocity_gradient_begin,  velocity_gradient_end,
                mass_change_rate,
                dCdRho,                   dCdRhoDot,
                std::begin( dCdGradRho ), std::end( dCdGradRho ),
                std::begin( dCdV ),       std::end( dCdV ),
                std::begin( dCdGradV ),   std::end( dCdGradV )
            );

            // Set the mass change rate
            mass_change_rate *= test_function;

            // Assemble the derivatives w.r.t. the density
            dCdRho += dCdRhoDot * dDensityDotdDensity;
            dCdRho *= test_function * interpolation_function;

            dCdRho += test_function * std::inner_product( std::begin( dCdGradRho ), std::end( dCdGradRho ), interpolation_function_gradient_begin, dCdRho_type( ) );

            // Assemble the derivatives w.r.t. the mesh displacement
            std::fill( dCdU_begin, dCdU_end, 0 );
            std::fill( dCdUMesh_begin, dCdUMesh_end, 0 );

            for ( unsigned int a = 0; a < dim; ++a ){

                *( dCdU_begin + a ) += test_function * dCdV[ a ] * dUDotdU * interpolation_function;

                *( dCdUMesh_begin + a ) += mass_change_rate * ( *( interpolation_function_gradient_begin + a ) );

                for ( unsigned int i = 0; i < dim; ++i ){

                    *( dCdU_begin + a ) += test_function * dCdGradV[ dim * a + i ] * ( *( interpolation_function_gradient_begin + i ) ) * dUDotdU;

                    *( dCdUMesh_begin + a ) -= test_function * dCdGradRho[ i ] * ( *( density_gradient_begin + a ) ) * ( *( interpolation_function_gradient_begin + i ) );

                    for ( unsigned int j = 0; j < dim; ++j ){

                        *( dCdUMesh_begin + a ) -= test_function * dCdGradV[ dim * j + i ] * ( *( velocity_gradient_begin + dim * j + a ) ) * ( *( interpolation_function_gradient_begin + i ) );

                    }

                }

            }

        }

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class mass_change_rate_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            mass_change_rate_iter_out mass_change_rate_begin,           mass_change_rate_iter_out mass_change_rate_end
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
             * \param &mass_change_rate_begin: The starting iterator of the rate of change of the mass \f$ c \f$
             * \param &mass_change_rate_end: The stopping iterator of the rate of change of the mass \f$ c \f$
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density dot arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( mass_change_rate_end - mass_change_rate_begin ), "The density and mass change rate arrays must be the same length" );

            unsigned int phase;

            for ( auto rho = density_begin; rho != density_end; ++rho ){

                phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfMass<dim>( *( density_begin + phase ), *( density_dot_begin + phase ),
                                           density_gradient_begin + dim * phase,        density_gradient_begin + dim * ( phase + 1 ),
                                           velocity_begin + dim * phase,                velocity_begin + dim * ( phase + 1 ),
                                           velocity_gradient_begin + dim * dim * phase, velocity_gradient_begin + dim * dim * ( phase + 1 ),
                                           *( mass_change_rate_begin + phase ) );

            }

        }

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            typename testFunction_type,
            class velocity_iter, class velocityGradient_iter, class mass_change_rate_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &test_function,
            mass_change_rate_iter_out mass_change_rate_begin,           mass_change_rate_iter_out mass_change_rate_end
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
             * \param &mass_change_rate_begin: The starting iterator of the rate of change of the mass \f$ c \f$
             * \param &mass_change_rate_end: The stopping iterator of the rate of change of the mass \f$ c \f$
             */

            computeBalanceOfMass<dim>(
                density_begin, density_end, density_dot_begin, density_dot_end,
                density_gradient_begin, density_gradient_end,
                velocity_begin, velocity_end,
                velocity_gradient_begin, velocity_gradient_end,
                mass_change_rate_begin, mass_change_rate_end
            );

            std::transform(
                mass_change_rate_begin, mass_change_rate_end, mass_change_rate_begin,
                std::bind(
                    std::multiplies< typename std::iterator_traits<mass_change_rate_iter_out>::value_type >( ),
                    std::placeholders::_1,
                    test_function
                )
            );

        }

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class mass_change_rate_iter_out,
            class dCdRho_iter_out, class dCdRhoDot_iter_out, class dCdGradRho_iter_out,
            class dCdV_iter_out, class dCdGradV_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                    const density_iter &density_end,
            const densityDot_iter &density_dot_begin,             const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            mass_change_rate_iter_out mass_change_rate_begin,     mass_change_rate_iter_out mass_change_rate_end,
            dCdRho_iter_out dCdRho_begin,                         dCdRho_iter_out dCdRho_end,
            dCdRhoDot_iter_out dCdRhoDot_begin,                   dCdRhoDot_iter_out dCdRhoDot_end,
            dCdGradRho_iter_out dCdGradRho_begin,                 dCdGradRho_iter_out dCdGradRho_end,
            dCdV_iter_out dCdV_begin,                             dCdV_iter_out dCdV_end,
            dCdGradV_iter_out dCdGradV_begin,                     dCdGradV_iter_out dCdGradV_end
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
             * \param &mass_change_rate_begin: The starting iterator of the rate of change of the mass \f$ c \f$
             * \param &mass_change_rate_end: The stopping iterator of the rate of change of the mass \f$ c \f$
             * \param &dCdRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdRhoDot_begin: The starting iterator of the derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dCdRhoDot_end: The stopping iterator of the derivative of the mass change rate w.r.t. partial time derivative of the density \f$ \frac{\partial \rho}{\partial t} \f$
             * \param &dCdGradRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dCdGradRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the density \f$ \rho_{,i} \f$
             * \param &dCdV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dCdV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the velocity \f$ v_{i} \f$
             * \param &dCdGradV_begin: The starting iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             * \param &dCdGradV_end: The stopping iterator of the derivative of the mass change rate w.r.t. the spatial gradient of the velocity \f$ v_{i,j} \f$
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density dot arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( mass_change_rate_end - mass_change_rate_begin ), "The density and mass change rate arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dCdRho_end - dCdRho_begin ), "The density and dCdRho arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dCdRhoDot_end - dCdRhoDot_begin ), "The density and dCdRhoDot arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dCdGradRho_end - dCdGradRho_begin ), "The density and dCdGradRho arrays must have consistent lengths" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dCdV_end - dCdV_begin ), "The density and dCdV arrays must have consistent lengths" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dCdGradV_end - dCdGradV_begin ), "The density and dCdGradV arrays must have consistent lengths" );

            unsigned int phase;

            for ( auto rho = density_begin; rho != density_end; ++rho ){

                phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfMass<dim>( *( density_begin + phase ),                *( density_dot_begin + phase ),
                                           density_gradient_begin + dim * phase,      density_gradient_begin + dim * ( phase + 1 ),
                                           velocity_begin + dim * phase,              velocity_begin + dim * ( phase + 1 ),
                                           velocity_gradient_begin + dim * dim * phase, velocity_gradient_begin + dim * dim * ( phase + 1 ),
                                           *( mass_change_rate_begin + phase ),       *( dCdRho_begin + phase ), *( dCdRhoDot_begin + phase ),
                                           dCdGradRho_begin + dim * phase,            dCdGradRho_begin + dim * ( phase + 1 ),
                                           dCdV_begin + dim * phase,                  dCdV_begin + dim * ( phase + 1 ),
                                           dCdGradV_begin + dim * dim * phase,         dCdGradV_begin + dim * dim * ( phase + 1 ) );

            }

        }

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            typename testFunction_type, typename interpolationFunction_type,
            class velocity_iter, class velocityGradient_iter, class interpolationFunctionGradient_iter,
            class mass_change_rate_iter_out,
            class dCdRho_iter_out, class dCdU_iter_out, class dCdUMesh_iter_out,
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
            mass_change_rate_iter_out mass_change_rate_begin,             mass_change_rate_iter_out mass_change_rate_end,
            dCdRho_iter_out dCdRho_begin,         dCdRho_iter_out dCdRho_end,
            dCdU_iter_out dCdU_begin,             dCdU_iter_out dCdU_end,
            dCdUMesh_iter_out dCdUMesh_begin,     dCdUMesh_iter_out dCdUMesh_end
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
             * \param &mass_change_rate_begin: The starting iterator of the rate of change of the mass \f$ c \f$
             * \param &mass_change_rate_end: The stopping iterator of the rate of change of the mass \f$ c \f$
             * \param &dCdRho_begin: The starting iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdRho_end: The stopping iterator of the derivative of the mass change rate w.r.t. density \f$ \rho \f$
             * \param &dCdU_begin: The starting iterator of the derivative of the mass change rate w.r.t. the phase displacement \f$ u_{i} \f$
             * \param &dCdU_end: The stopping iterator of the derivative of the mass change rate w.r.t. the phase displacement \f$ u_{i} \f$
             * \param &dCdUMesh_begin: The starting iterator of the derivative of the mass change rate w.r.t. the mesh displacement \f$ u_{i} \f$
             * \param &dCdUMesh_end: The stopping iterator of the derivative of the mass change rate w.r.t. the mesh displacement \f$ u_{i} \f$
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density dot arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient arrays are of inconsistent length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( mass_change_rate_end - mass_change_rate_begin ), "The density and mass change rate arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dDensityDotdDensity_end - dDensityDotdDensity_begin ), "The density and the derivative of the density time derivative w.r.t. the density must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dUDotdU_end - dUDotdU_begin ), "The density and the derivative of the dof time derivative w.r.t. the dof must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(             ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dCdRho_end - dCdRho_begin ), "The density and dCdRho arrays must be the same length" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dCdU_end - dCdU_begin ), "The density and dCdU arrays must have consistent lengths" );

            TARDIGRADE_ERROR_TOOLS_CHECK(       dim * ( unsigned int )( density_end - density_begin ) == ( unsigned int )( dCdUMesh_end - dCdUMesh_begin ), "The density and dCdUMesh arrays must have consistent lengths" );

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
                    *( mass_change_rate_begin + phase ),
                    *( dCdRho_begin + phase ),
                    dCdU_begin + dim * phase,                    dCdU_begin + dim * ( phase + 1 ),
                    dCdUMesh_begin + dim * phase,                dCdUMesh_begin + dim * ( phase + 1 )
                );

            }

        }

    }

}
