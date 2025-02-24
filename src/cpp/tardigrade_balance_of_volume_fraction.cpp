/**
  ******************************************************************************
  * \file tardigrade_balance_of_volume_fraction.cpp
  ******************************************************************************
  * The source file for the equations associated with the balance of the
  * volume fraction
  ******************************************************************************
  */

#include "tardigrade_balance_of_volume_fraction.h"
#include<numeric>
#include<algorithm>
#include<functional>

namespace tardigradeBalanceEquations{

    namespace balanceOfVolumeFraction{

        template<
            int dim,
            typename density_type, class velocity_iter,
            typename volume_fraction_type, typename volume_fraction_dot_type,
            class volume_fraction_gradient_iter,
            typename mass_change_rate_type,
            typename rest_density_type,
            typename trace_mass_change_velocity_gradient_type,
            typename test_function_type,
            typename result_type
        >
        void computeBalanceOfVolumeFraction(
            const density_type                  &density,
            const velocity_iter                 &velocity_begin,                 const velocity_iter                 &velocity_end,
            const volume_fraction_type          &volume_fraction,
            const volume_fraction_dot_type      &volume_fraction_dot,
            const volume_fraction_gradient_iter &volume_fraction_gradient_begin, const volume_fraction_gradient_iter &volume_fraction_gradient_end,
            const mass_change_rate_type         &mass_change_rate,
            const rest_density_type             &rest_density,
            const trace_mass_change_velocity_gradient_type &trace_mass_change_velocity_gradient,
            const test_function_type &test_function,
            result_type &result,
            const double volume_fraction_tolerance
        ){
            /*!
             * Compute the balance of the volume fraction for a reacting continuum
             * 
             * \param &density: The current apparent density \f$ \left( \rho^{\alpha} \right) \f$
             * \param &velocity_begin: The starting iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &velocity_end: The stopping iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &volume_fraction: The current volume fraction \f$ \left( \phi^{\alpha} \right) \f$
             * \param &volume_fraction_dot: The partial time derivative of the current volume fraction \f$ \left( \frac{\partial \phi^{\alpha}}{\partial t} \right) \f$
             * \param &volume_fraction_gradient_begin: The starting iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &volume_fraction_gradient_end: The stopping iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &mass_change_rate: The current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &rest_density: The rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &trace_mass_change_velocity_gradient: The current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &test_function: The current value of the test function
             * \param &result: The value of the balance of volume fraction
             * \param volume_fraction_tolerance: The tolerance of the volume fraction where if it is less than the tolerance, the
             *     true density is assumed to be the rest density.
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(
                ( unsigned int )( velocity_end - velocity_begin ) == dim,
                "The velocity must be the same size as the spatial dimension"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                ( unsigned int )( volume_fraction_gradient_end - volume_fraction_gradient_begin ) == dim,
                "The volume fraction gradient must be the same size as the spatial dimension"
            );

            using true_density_type = decltype( std::declval< density_type& >( ) / std::declval< volume_fraction_type& >( ) );

            result = volume_fraction_dot + volume_fraction * trace_mass_change_velocity_gradient;

            for ( unsigned int i = 0; i < dim; ++i ){

                result += ( *( volume_fraction_gradient_begin + i ) ) * ( *( velocity_begin + i ) );

            }

            // Compute the true density
            true_density_type true_density = rest_density;

            if ( volume_fraction >= volume_fraction_tolerance ){

                true_density = ( density / volume_fraction );

            }

            result -= mass_change_rate / true_density;

            result *= test_function;

        }

        template<
            int dim, int mass_change_rate_index, int trace_mass_change_velocity_gradient_index,
            typename density_type, class velocity_iter,
            typename volume_fraction_type, typename volume_fraction_dot_type,
            class volume_fraction_gradient_iter,
            class material_response_iter,
            typename rest_density_type,
            typename test_function_type,
            typename result_type
        >
        void computeBalanceOfVolumeFraction(
            const density_type                  &density,
            const velocity_iter                 &velocity_begin,                 const velocity_iter                 &velocity_end,
            const volume_fraction_type          &volume_fraction,
            const volume_fraction_dot_type      &volume_fraction_dot,
            const volume_fraction_gradient_iter &volume_fraction_gradient_begin, const volume_fraction_gradient_iter &volume_fraction_gradient_end,
            const material_response_iter        &material_response_begin,
            const material_response_iter        &material_response_end,
            const rest_density_type             &rest_density,
            const test_function_type &test_function,
            result_type &result,
            const double volume_fraction_tolerance
        ){
            /*!
             * Compute the balance of the volume fraction for a reacting continuum
             * 
             * \param &density: The current apparent density \f$ \left( \rho^{\alpha} \right) \f$
             * \param &velocity_begin: The starting iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &velocity_end: The stopping iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &volume_fraction: The current volume fraction \f$ \left( \phi^{\alpha} \right) \f$
             * \param &volume_fraction_dot: The partial time derivative of the current volume fraction \f$ \left( \frac{\partial \phi^{\alpha}}{\partial t} \right) \f$
             * \param &volume_fraction_gradient_begin: The starting iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &volume_fraction_gradient_end: The stopping iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &material_response_begin: The starting iterator for the material response vector. This vector defines the mass change rate and the
             *     trace of the mass-change velocity gradient which are located at their respective indices (see template parameters)
             * \param &material_response_end: The stopping iterator for the material response vector. This vector defines the mass change rate and the
             *     trace of the mass-change velocity gradient which are located at their respective indices (see template parameters)
             * \param &rest_density: The rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &test_function: The current value of the test function
             * \param &result: The value of the balance of volume fraction
             * \param volume_fraction_tolerance: The tolerance of the volume fraction where if it is less than the tolerance, the
             *     true density is assumed to be the rest density.
             */

//            std::cout << "density                  : " << density << "\n";
//            std::cout << "velocity                 : "; for ( auto _v = velocity_begin; _v != velocity_end; ++_v ){ std::cout << *_v << " "; } std::cout << "\n";
//            std::cout << "volume_fraction          : " << volume_fraction << "\n";
//            std::cout << "volume_fraction_dot      : " << volume_fraction_dot << "\n";
//            std::cout << "volume_fraction_gradient : "; for ( auto _v = volume_fraction_gradient_begin; _v != volume_fraction_gradient_end; ++_v ){ std::cout << *_v << " "; } std::cout << "\n";
//            std::cout << "MR                       : "; for ( auto _v = material_response_begin; _v != material_response_end; ++_v ){ std::cout << *_v << " "; } std::cout << "\n";

            computeBalanceOfVolumeFraction<dim>(
                density,
                velocity_begin, velocity_end,
                volume_fraction, volume_fraction_dot, volume_fraction_gradient_begin, volume_fraction_gradient_end,
                *( material_response_begin + mass_change_rate_index ),
                rest_density,
                *( material_response_begin + trace_mass_change_velocity_gradient_index ),
                test_function, result, volume_fraction_tolerance
            );

        }

        template<
            int dim,
            typename density_type, class velocity_iter,
            typename volume_fraction_type, typename volume_fraction_dot_type,
            class volume_fraction_gradient_iter,
            typename mass_change_rate_type,
            typename rest_density_type,
            typename trace_mass_change_velocity_gradient_type,
            typename test_function_type,
            typename interpolation_function_type,
            class interpolation_function_gradient_iter,
            typename dUDotdU_type,
            typename dVolumeFractionDotdVolumeFraction_type,
            typename result_type,
            typename dRdRho_type,
            class dRdU_iter,
            typename dRdVolumeFraction_type,
            typename dRdC_type,
            typename dRdTraceVA_type,
            class dRdUMesh_iter
        >
        void computeBalanceOfVolumeFraction(
            const density_type                  &density,
            const velocity_iter                 &velocity_begin,                 const velocity_iter                 &velocity_end,
            const volume_fraction_type          &volume_fraction,
            const volume_fraction_dot_type      &volume_fraction_dot,
            const volume_fraction_gradient_iter &volume_fraction_gradient_begin, const volume_fraction_gradient_iter &volume_fraction_gradient_end,
            const mass_change_rate_type         &mass_change_rate,
            const rest_density_type             &rest_density,
            const trace_mass_change_velocity_gradient_type &trace_mass_change_velocity_gradient,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const dUDotdU_type dUDotdU, const dVolumeFractionDotdVolumeFraction_type dVolumeFractionDotdVolumeFraction,
            result_type &result,
            dRdRho_type &dRdRho,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end,
            dRdVolumeFraction_type &dRdVolumeFraction,
            dRdC_type &dRdC,
            dRdTraceVA_type &dRdTraceVA,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end,
            const double volume_fraction_tolerance
        ){
            /*!
             * Compute the balance of the volume fraction for a reacting continuum
             * 
             * \param &density: The current apparent density \f$ \left( \rho^{\alpha} \right) \f$
             * \param &velocity_begin: The starting iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &velocity_end: The stopping iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &volume_fraction: The current volume fraction \f$ \left( \phi^{\alpha} \right) \f$
             * \param &volume_fraction_dot: The partial time derivative of the current volume fraction \f$ \left( \frac{\partial \phi^{\alpha}}{\partial t} \right) \f$
             * \param &volume_fraction_gradient_begin: The starting iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &volume_fraction_gradient_end: The stopping iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &mass_change_rate: The current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &rest_density: The rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &trace_mass_change_velocity_gradient: The current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &test_function: The current value of the test function
             * \param &interpolation_function: The current value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the current value of the spatial gradient of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the current value of the spatial gradient of the interpolation function
             * \param &dUDotdU: The derivative of the phase velocity w.r.t. the phase displacement dof
             * \param &dVolumeFractionDotdVolumeFraction: The derivative of the partial time derivative of the volume fraction phase w.r.t. the volume fraction
             * \param &result: The value of the balance of volume fraction
             * \param &dRdRho: The derivative of the residual w.r.t. the apparent density
             * \param &dRdU_begin: The starting iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdU_end: The stopping iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdVolumeFraction: The derivative of the residual w.r.t. the volume fraction
             * \param &dRdC: The derivative of the residual w.r.t. the mass change rate per unit volume
             * \param &dRdTraceVA: The derivative of the residual w.r.t. the trace of the velocity gradient associated with the change in mass
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param volume_fraction_tolerance: The tolerance of the volume fraction where if it is less than the tolerance, the
             *     true density is assumed to be the rest density.
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(
                ( unsigned int )( velocity_end - velocity_begin ) == dim,
                "The velocity must be the same size as the spatial dimension"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                ( unsigned int )( volume_fraction_gradient_end - volume_fraction_gradient_begin ) == dim,
                "The volume fraction gradient must be the same size as the spatial dimension"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                ( unsigned int )( dRdU_end - dRdU_begin ) == dim,
                "The derivative of the residual w.r.t. the phase displacement dof must be the same size as the spatial dimension"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                ( unsigned int )( dRdUMesh_end - dRdUMesh_begin ) == dim,
                "The derivative of the residual w.r.t. the mesh displacement must be the same size as the spatial dimension"
            );

            using true_density_type = decltype( std::declval< density_type& >( ) / std::declval< volume_fraction_type& >( ) );

            result = volume_fraction_dot + volume_fraction * trace_mass_change_velocity_gradient;

            dRdVolumeFraction = interpolation_function * ( dVolumeFractionDotdVolumeFraction + trace_mass_change_velocity_gradient );

            std::fill( dRdU_begin, dRdU_end, 0 );

            std::fill( dRdUMesh_begin, dRdUMesh_end, 0 );

            for ( unsigned int i = 0; i < dim; ++i ){

                result += ( *( volume_fraction_gradient_begin + i ) ) * ( *( velocity_begin + i ) );

                *( dRdU_begin + i ) += test_function * ( *( volume_fraction_gradient_begin + i ) ) * interpolation_function * dUDotdU;

                dRdVolumeFraction += ( *( interpolation_function_gradient_begin + i ) ) * ( *( velocity_begin + i ) );

                for ( unsigned int a = 0; a < dim; ++a ){

                    *( dRdUMesh_begin + a ) -= test_function * ( *( volume_fraction_gradient_begin + a ) ) * ( *( interpolation_function_gradient_begin + i ) ) * ( *( velocity_begin + i ) );

                }

            }

            // Compute the true density
            true_density_type true_density = rest_density;

            density_type dGammadVolumeFraction = 0;
            volume_fraction_type dGammadRho    = 0;

            if ( volume_fraction >= volume_fraction_tolerance ){

                true_density = ( density / volume_fraction );

                dGammadVolumeFraction = -density / ( volume_fraction * volume_fraction ) * interpolation_function;

                dGammadRho = interpolation_function / volume_fraction;

            }

            result -= mass_change_rate / true_density;

            result *= test_function;

            dRdRho = test_function * mass_change_rate / ( true_density * true_density ) * dGammadRho;

            dRdVolumeFraction += mass_change_rate / ( true_density * true_density ) * dGammadVolumeFraction;

            dRdVolumeFraction *= test_function;

            for ( unsigned int a = 0; a < dim; ++a ){

                *( dRdUMesh_begin + a ) += result * ( *( interpolation_function_gradient_begin + a ) );

            }

            dRdC = -test_function / true_density;

            dRdTraceVA = test_function * volume_fraction;

        }

        template<
            int dim,
            class density_iter, class velocity_iter,
            class volume_fraction_iter, class volume_fraction_dot_iter,
            class volume_fraction_gradient_iter,
            class mass_change_rate_iter,
            class rest_density_iter,
            class trace_mass_change_velocity_gradient_iter,
            typename test_function_type,
            class result_iter
        >
        void computeBalanceOfVolumeFraction(
            const density_iter                  &density_begin,                  const density_iter                  &density_end,
            const velocity_iter                 &velocity_begin,                 const velocity_iter                 &velocity_end,
            const volume_fraction_iter          &volume_fraction_begin,          const volume_fraction_iter          &volume_fraction_end,
            const volume_fraction_dot_iter      &volume_fraction_dot_begin,      const volume_fraction_dot_iter      &volume_fraction_dot_end,
            const volume_fraction_gradient_iter &volume_fraction_gradient_begin, const volume_fraction_gradient_iter &volume_fraction_gradient_end,
            const mass_change_rate_iter         &mass_change_rate_begin,         const mass_change_rate_iter         &mass_change_rate_end,
            const rest_density_iter             &rest_density_begin,             const rest_density_iter             &rest_density_end,
            const trace_mass_change_velocity_gradient_iter &trace_mass_change_velocity_gradient_begin,
            const trace_mass_change_velocity_gradient_iter &trace_mass_change_velocity_gradient_end,
            const test_function_type &test_function,
            result_iter result_begin, result_iter result_end,
            const double volume_fraction_tolerance
        ){
            /*!
             * Compute the balance of the volume fraction for a multiphase reacting continuum
             * 
             * \param &density_begin: The starting iterator of the current apparent density \f$ \left( \rho^{\alpha} \right) \f$
             * \param &density_end: The stopping iterator of the current apparent density \f$ \left( \rho^{\alpha} \right) \f$
             * \param &velocity_begin: The starting iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &velocity_end: The stopping iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &volume_fraction_begin: The starting iterator of the current volume fraction \f$ \left( \phi^{\alpha} \right) \f$
             * \param &volume_fraction_end: The stopping iterator of the current volume fraction \f$ \left( \phi^{\alpha} \right) \f$
             * \param &volume_fraction_dot_begin: The starting iterator of the partial time derivative of the current volume fraction \f$ \left( \frac{\partial \phi^{\alpha}}{\partial t} \right) \f$
             * \param &volume_fraction_dot_end: The stopping iterator of the partial time derivative of the current volume fraction \f$ \left( \frac{\partial \phi^{\alpha}}{\partial t} \right) \f$
             * \param &volume_fraction_gradient_begin: The starting iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &volume_fraction_gradient_end: The stopping iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &mass_change_rate_begin: The starting iterator of the current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &mass_change_rate_end: The stopping iterator of the current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &rest_density_begin: The starting iterator of the rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &rest_density_end: The stopping iterator of the rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &trace_mass_change_velocity_gradient_begin: The starting iterator of the current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &trace_mass_change_velocity_gradient_end: The stopping iterator of the current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &test_function: The current value of the test function
             * \param &result_begin: The starting iterator of the value of the balance of volume fraction
             * \param &result_end: The stopping iterator of the value of the balance of volume fraction
             * \param volume_fraction_tolerance: The tolerance of the volume fraction where if it is less than the tolerance, the
             *     true density is assumed to be the rest density.
             */

            TARDIGRADE_ERROR_TOOLS_EVAL(
                unsigned int nphases = ( unsigned int )( density_end - density_begin );
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( velocity_end - velocity_begin ),
                "The velocity size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( volume_fraction_end - volume_fraction_begin ),
                "The volume fraction size is not equal to the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( volume_fraction_dot_end - volume_fraction_dot_begin ),
                "The partial derivative of the volume fraction w.r.t. time size is not equal to the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( volume_fraction_gradient_end - volume_fraction_gradient_begin ),
                "The volume fraction gradient size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( mass_change_rate_end - mass_change_rate_begin ),
                "The mass change rate size is not equal to the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( trace_mass_change_velocity_gradient_end - trace_mass_change_velocity_gradient_begin ),
                "The trace of the mass change velocity gradient size is not equal to the number of phases"
            );

            for ( auto rho = std::pair< unsigned int, density_iter >( 0, density_begin ); rho.second != density_end; ++rho.first, ++rho.second ){

                computeBalanceOfVolumeFraction<dim>(
                    *rho.second,
                    velocity_begin + dim * rho.first,                 velocity_begin + dim * ( rho.first + 1 ),
                    *( volume_fraction_begin + rho.first ),           *( volume_fraction_dot_begin + rho.first ),
                    volume_fraction_gradient_begin + dim * rho.first, volume_fraction_gradient_begin + dim * ( rho.first + 1 ),
                    *( mass_change_rate_begin + rho.first ),
                    *( rest_density_begin + rho.first ),
                    *( trace_mass_change_velocity_gradient_begin + rho.first ),
                    test_function,
                    *( result_begin + rho.first ),
                    volume_fraction_tolerance
                );

            }

        }

        template<
            int dim, int mass_change_rate_index, int trace_mass_change_velocity_gradient_index,
            class density_iter, class velocity_iter,
            class volume_fraction_iter, class volume_fraction_dot_iter,
            class volume_fraction_gradient_iter,
            class material_response_iter,
            class rest_density_iter,
            typename test_function_type,
            class result_iter
        >
        void computeBalanceOfVolumeFraction(
            const density_iter                  &density_begin,                  const density_iter                  &density_end,
            const velocity_iter                 &velocity_begin,                 const velocity_iter                 &velocity_end,
            const volume_fraction_iter          &volume_fraction_begin,          const volume_fraction_iter          &volume_fraction_end,
            const volume_fraction_dot_iter      &volume_fraction_dot_begin,      const volume_fraction_dot_iter      &volume_fraction_dot_end,
            const volume_fraction_gradient_iter &volume_fraction_gradient_begin, const volume_fraction_gradient_iter &volume_fraction_gradient_end,
            const material_response_iter        &material_response_begin,        const material_response_iter        &material_response_end,
            const rest_density_iter             &rest_density_begin,             const rest_density_iter             &rest_density_end,
            const test_function_type &test_function,
            result_iter result_begin, result_iter result_end,
            const double volume_fraction_tolerance
        ){
            /*!
             * Compute the balance of the volume fraction for a multiphase reacting continuum
             * 
             * \param &density_begin: The starting iterator of the current apparent density \f$ \left( \rho^{\alpha} \right) \f$
             * \param &density_end: The stopping iterator of the current apparent density \f$ \left( \rho^{\alpha} \right) \f$
             * \param &velocity_begin: The starting iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &velocity_end: The stopping iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &volume_fraction_begin: The starting iterator of the current volume fraction \f$ \left( \phi^{\alpha} \right) \f$
             * \param &volume_fraction_end: The stopping iterator of the current volume fraction \f$ \left( \phi^{\alpha} \right) \f$
             * \param &volume_fraction_dot_begin: The starting iterator of the partial time derivative of the current volume fraction \f$ \left( \frac{\partial \phi^{\alpha}}{\partial t} \right) \f$
             * \param &volume_fraction_dot_end: The stopping iterator of the partial time derivative of the current volume fraction \f$ \left( \frac{\partial \phi^{\alpha}}{\partial t} \right) \f$
             * \param &volume_fraction_gradient_begin: The starting iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &volume_fraction_gradient_end: The stopping iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &material_response_begin: The starting iterator for the material response vector. This vector defines the mass change rate and the
             *     trace of the mass-change velocity gradient which are located at their respective indices (see template parameters)
             * \param &material_response_end: The stopping iterator for the material response vector. This vector defines the mass change rate and the
             *     trace of the mass-change velocity gradient which are located at their respective indices (see template parameters)
             * \param &rest_density_begin: The starting iterator of the rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &rest_density_end: The stopping iterator of the rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &test_function: The current value of the test function
             * \param &result_begin: The starting iterator of the value of the balance of volume fraction
             * \param &result_end: The stopping iterator of the value of the balance of volume fraction
             * \param volume_fraction_tolerance: The tolerance of the volume fraction where if it is less than the tolerance, the
             *     true density is assumed to be the rest density.
             */

            const unsigned int nphases = ( unsigned int )( density_end - density_begin );
            const unsigned int material_response_size = ( unsigned int )( material_response_end - material_response_begin ) / nphases;

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( velocity_end - velocity_begin ),
                "The velocity size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( volume_fraction_end - volume_fraction_begin ),
                "The volume fraction size is not equal to the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( volume_fraction_dot_end - volume_fraction_dot_begin ),
                "The partial derivative of the volume fraction w.r.t. time size is not equal to the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( volume_fraction_gradient_end - volume_fraction_gradient_begin ),
                "The volume fraction gradient size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size == ( unsigned int )( material_response_end - material_response_begin ),
                "The material response vector size must be a integer multiple of the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                ( material_response_size > mass_change_rate_index ),
                "The material response vector size must larger than the index for the mass change rate"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                ( material_response_size > trace_mass_change_velocity_gradient_index ),
                "The material response vector size must larger than the index for the trace of the mass change rate velocity gradient"
            )

            for ( auto rho = std::pair< unsigned int, density_iter >( 0, density_begin ); rho.second != density_end; ++rho.first, ++rho.second ){

                computeBalanceOfVolumeFraction<dim, mass_change_rate_index, trace_mass_change_velocity_gradient_index>(
                    *rho.second,
                    velocity_begin + dim * rho.first,                 velocity_begin + dim * ( rho.first + 1 ),
                    *( volume_fraction_begin + rho.first ),           *( volume_fraction_dot_begin + rho.first ),
                    volume_fraction_gradient_begin + dim * rho.first, volume_fraction_gradient_begin + dim * ( rho.first + 1 ),
                    material_response_begin + material_response_size * rho.first,
                    material_response_begin + material_response_size * ( rho.first + 1 ),
                    *( rest_density_begin + rho.first ),
                    test_function,
                    *( result_begin + rho.first ),
                    volume_fraction_tolerance
                );

            }

        }

        template<
            int dim,
            class density_iter, class velocity_iter,
            class volume_fraction_iter, class volume_fraction_dot_iter,
            class volume_fraction_gradient_iter,
            class mass_change_rate_iter,
            class rest_density_iter,
            class trace_mass_change_velocity_gradient_iter,
            typename test_function_type,
            typename interpolation_function_type,
            class interpolation_function_gradient_iter,
            typename dUDotdU_type,
            typename dVolumeFractionDotdVolumeFraction_type,
            class result_iter,
            class dRdRho_iter,
            class dRdU_iter,
            class dRdVolumeFraction_iter,
            class dRdC_iter,
            class dRdTraceVA_iter,
            class dRdUMesh_iter
        >
        void computeBalanceOfVolumeFraction(
            const density_iter                  &density_begin,                  const density_iter                  &density_end,
            const velocity_iter                 &velocity_begin,                 const velocity_iter                 &velocity_end,
            const volume_fraction_iter          &volume_fraction_begin,          const volume_fraction_iter          &volume_fraction_end,
            const volume_fraction_dot_iter      &volume_fraction_dot_begin,      const volume_fraction_dot_iter      &volume_fraction_dot_end,
            const volume_fraction_gradient_iter &volume_fraction_gradient_begin, const volume_fraction_gradient_iter &volume_fraction_gradient_end,
            const mass_change_rate_iter         &mass_change_rate_begin,         const mass_change_rate_iter         &mass_change_rate_end,
            const rest_density_iter             &rest_density_begin,             const rest_density_iter             &rest_density_end,
            const trace_mass_change_velocity_gradient_iter &trace_mass_change_velocity_gradient_begin,
            const trace_mass_change_velocity_gradient_iter &trace_mass_change_velocity_gradient_end,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const dUDotdU_type dUDotdU,
            const dVolumeFractionDotdVolumeFraction_type dVolumeFractionDotdVolumeFraction,
            result_iter result_begin, result_iter result_end,
            dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin, dRdVolumeFraction_iter dRdVolumeFraction_end,
            dRdC_iter dRdC_begin, dRdC_iter dRdC_end,
            dRdTraceVA_iter dRdTraceVA_begin, dRdTraceVA_iter dRdTraceVA_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end,
            const double volume_fraction_tolerance
        ){
            /*!
             * Compute the balance of the volume fraction for a multiphase reacting continuum
             * 
             * \param &density_begin: The starting iterator of the current apparent density \f$ \left( \rho^{\alpha} \right) \f$
             * \param &density_end: The stopping iterator of the current apparent density \f$ \left( \rho^{\alpha} \right) \f$
             * \param &velocity_begin: The starting iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &velocity_end: The stopping iterator of the phase velocity \f$ \left( v_i^{\alpha} \right) \f$
             * \param &volume_fraction_begin: The starting iterator of the current volume fraction \f$ \left( \phi^{\alpha} \right) \f$
             * \param &volume_fraction_end: The stopping iterator of the current volume fraction \f$ \left( \phi^{\alpha} \right) \f$
             * \param &volume_fraction_dot_begin: The starting iterator of the partial time derivative of the current volume fraction \f$ \left( \frac{\partial \phi^{\alpha}}{\partial t} \right) \f$
             * \param &volume_fraction_dot_end: The stopping iterator of the partial time derivative of the current volume fraction \f$ \left( \frac{\partial \phi^{\alpha}}{\partial t} \right) \f$
             * \param &volume_fraction_gradient_begin: The starting iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &volume_fraction_gradient_end: The stopping iterator of the spatial gradient of the volume fraction \f$ \left( \phi_{,i}^{\alpha} \right) \f$
             * \param &mass_change_rate_begin: The starting iterator of the current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &mass_change_rate_end: The stopping iterator of the current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &rest_density_begin: The starting iterator of the rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &rest_density_end: The stopping iterator of the rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &trace_mass_change_velocity_gradient_begin: The starting iterator of the current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &trace_mass_change_velocity_gradient_end: The stopping iterator of the current mass change rate per current volume \f$ \left( c^{\alpha} \right) \f$
             * \param &test_function: The current value of the test function
             * \param &interpolation_function: The current value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the current value of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the current value of the interpolation function
             * \param &dUDotdU: The derivative of the phase velocity w.r.t. the phase displacement dof
             * \param &dVolumeFractionDotdVolumeFraction: The derivative of the partial time derivative of the volume fraction phase w.r.t. the volume fraction
             * \param &result_begin: The starting iterator of the value of the balance of volume fraction
             * \param &result_end: The stopping iterator of the value of the balance of volume fraction
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdU_begin: The starting iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdU_end: The stopping iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdVolumeFraction_begin: The starting iterator of the derivative of the residual w.r.t. the volume fraction
             * \param &dRdVolumeFraction_end: The stopping iterator of the derivative of the residual w.r.t. the volume fraction
             * \param &dRdC_begin: The starting iterator of the derivative of the residual w.r.t. the mass change rate per unit volume
             * \param &dRdC_end: The stopping iterator of the derivative of the residual w.r.t. the mass change rate per unit volume
             * \param &dRdTraceVA_begin: The starting iterator of the derivative of the residual w.r.t. the trace of the velocity gradient associated with the change in mass
             * \param &dRdTraceVA_end: The stopping iterator of the derivative of the residual w.r.t. the trace of the velocity gradient associated with the change in mass
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param volume_fraction_tolerance: The tolerance of the volume fraction where if it is less than the tolerance, the
             *     true density is assumed to be the rest density.
             */

            TARDIGRADE_ERROR_TOOLS_EVAL(
                unsigned int nphases = ( unsigned int )( density_end - density_begin );
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( velocity_end - velocity_begin ),
                "The velocity size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( volume_fraction_end - volume_fraction_begin ),
                "The volume fraction size is not equal to the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( volume_fraction_dot_end - volume_fraction_dot_begin ),
                "The partial derivative of the volume fraction w.r.t. time size is not equal to the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( volume_fraction_gradient_end - volume_fraction_gradient_begin ),
                "The volume fraction gradient size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( mass_change_rate_end - mass_change_rate_begin ),
                "The mass change rate size is not equal to the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( trace_mass_change_velocity_gradient_end - trace_mass_change_velocity_gradient_begin ),
                "The trace of the mass change velocity gradient size is not equal to the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( dRdRho_end - dRdRho_begin ),
                "The dRdRho size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( dRdU_end - dRdU_begin ),
                "The dRdU size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( dRdVolumeFraction_end - dRdVolumeFraction_begin ),
                "The dRdVolumeFraction size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( dRdC_end - dRdC_begin ),
                "The dRdC size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( dRdTraceVA_end - dRdTraceVA_begin ),
                "The dRdTraceVA size is inconsistent with the number of phases"
            );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( dRdUMesh_end - dRdUMesh_begin ),
                "The dRdUMesh size is inconsistent with the number of phases"
            );

            for ( auto rho = std::pair< unsigned int, density_iter >( 0, density_begin ); rho.second != density_end; ++rho.first, ++rho.second ){

                computeBalanceOfVolumeFraction<dim>(
                    *rho.second,
                    velocity_begin + dim * rho.first,                 velocity_begin + dim * ( rho.first + 1 ),
                    *( volume_fraction_begin + rho.first ),           *( volume_fraction_dot_begin + rho.first ),
                    volume_fraction_gradient_begin + dim * rho.first, volume_fraction_gradient_begin + dim * ( rho.first + 1 ),
                    *( mass_change_rate_begin + rho.first ),
                    *( rest_density_begin + rho.first ),
                    *( trace_mass_change_velocity_gradient_begin + rho.first ),
                    test_function,
                    interpolation_function,
                    interpolation_function_gradient_begin, interpolation_function_gradient_end,
                    dUDotdU, dVolumeFractionDotdVolumeFraction,
                    *( result_begin + rho.first ),
                    *( dRdRho_begin + rho.first ),
                    dRdU_begin + dim * rho.first, dRdU_begin + dim * ( rho.first + 1 ),
                    *( dRdVolumeFraction_begin + rho.first ),
                    *( dRdC_begin + rho.first ),
                    *( dRdTraceVA_begin + rho.first ),
                    dRdUMesh_begin + dim * rho.first, dRdUMesh_begin + dim * ( rho.first + 1 ),
                    volume_fraction_tolerance
                );

            }

        }

    }

}
