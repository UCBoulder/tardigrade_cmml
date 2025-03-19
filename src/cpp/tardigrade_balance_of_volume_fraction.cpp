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
            int dim, int material_response_dim, int mass_change_rate_index, int trace_mass_change_velocity_gradient_index,
            int material_response_num_dof,
            typename density_type, class velocity_iter,
            typename volume_fraction_type, typename volume_fraction_dot_type,
            class volume_fraction_gradient_iter,
            class material_response_iter,
            class material_response_jacobian_iter,
            typename rest_density_type,
            typename test_function_type,
            typename interpolation_function_type, class interpolation_function_gradient_iter,
            class full_material_response_dof_gradient_iter,
            typename dUDotdU_type, typename dVolumeFractionDotdVolumeFraction_type,
            typename result_type,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter, class dRdTheta_iter,
            class dRdE_iter, class dRdZ_iter, class dRdVolumeFraction_iter, class dRdUMesh_iter,
            int density_index        ,
            int displacement_index   ,
            int velocity_index       ,
            int temperature_index    ,
            int internal_energy_index,
            int additional_dof_index 
        >
        void computeBalanceOfVolumeFraction(
            const density_type                  &density,
            const velocity_iter                 &velocity_begin,                 const velocity_iter                 &velocity_end,
            const volume_fraction_type          &volume_fraction,
            const volume_fraction_dot_type      &volume_fraction_dot,
            const volume_fraction_gradient_iter &volume_fraction_gradient_begin, const volume_fraction_gradient_iter &volume_fraction_gradient_end,
            const material_response_iter          &material_response_begin,
            const material_response_iter          &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const rest_density_type               &rest_density,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU, const dVolumeFractionDotdVolumeFraction_type dVolumeFractionDotdVolumeFraction,
            const unsigned int phase,
            result_type &result,
            dRdRho_iter   dRdRho_begin,                     dRdRho_iter            dRdRho_end,
            dRdU_iter     dRdU_begin,                       dRdU_iter              dRdU_end,
            dRdW_iter     dRdW_begin,                       dRdW_iter              dRdW_end,
            dRdTheta_iter dRdTheta_begin,                   dRdTheta_iter          dRdTheta_end,
            dRdE_iter     dRdE_begin,                       dRdE_iter              dRdE_end,
            dRdZ_iter     dRdZ_begin,                       dRdZ_iter              dRdZ_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin, dRdVolumeFraction_iter dRdVolumeFraction_end,
            dRdUMesh_iter          dRdUMesh_begin,          dRdUMesh_iter          dRdUMesh_end,
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
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &rest_density: The rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &test_function: The current value of the test function
             * \param &interpolation_function: The current value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the current value of the spatial gradient of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the current value of the spatial gradient of the interpolation function
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the material response dof vector
             * \param &dUDotdU: The derivative of the phase velocity w.r.t. the phase displacement dof
             * \param &dVolumeFractionDotdVolumeFraction: The derivative of the partial time derivative of the volume fraction phase w.r.t. the volume fraction
             * \param &phase: The current active phase
             * \param &result: The value of the balance of volume fraction
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdU_begin: The starting iterator of the derivative of the residual w.r.t. the phase spatial DOF associated with velocity
             * \param &dRdU_end: The stopping iterator of the derivative of the residual w.r.t. the phase spatial DOF associated with velocity
             * \param &dRdW_begin: The starting iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdW_end: The stopping iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdTheta_begin: The starting iterator of the derivative of the residual w.r.t. the phase temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the residual w.r.t. the phase temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the residual w.r.t. the phase internal energy
             * \param &dRdE_end: The stopping iterator of the derivative of the residual w.r.t. the phase internal energy
             * \param &dRdZ_begin: The starting iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdZ_end: The stopping iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdVolumeFraction_begin: The starting iterator of the derivative of the residual w.r.t. the volume fraction
             * \param &dRdVolumeFraction_end: The stopping iterator of the derivative of the residual w.r.t. the volume fraction
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param volume_fraction_tolerance: The tolerance of the volume fraction where if it is less than the tolerance, the
             *     true density is assumed to be the rest density.
             */

            const unsigned int nphases = ( unsigned int )( dRdRho_end - dRdRho_begin );
            constexpr unsigned int num_phase_dof = 3 + 2 * material_response_dim;
            const unsigned int num_additional_dof = material_response_num_dof - num_phase_dof;

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( dRdU_end - dRdU_begin ), "The dRdU must be a consistent size with the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( dRdW_end - dRdW_begin ), "The dRdW must be a consistent size with the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( dRdTheta_end - dRdTheta_begin ), "The dRdTheta must be the same size as the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( dRdE_end - dRdE_begin ), "The dRdE must be the same size as the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( dRdVolumeFraction_end - dRdVolumeFraction_begin ), "The dRdVolumeFraction must be the same size as the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim == ( unsigned int )( dRdUMesh_end - dRdUMesh_begin ), "The dRdUMesh must be the same size as the dimension"
            )

            result_type dRdC_phase, dRdTraceVA_phase;

            std::fill( dRdRho_begin,            dRdRho_end,            0 );
            std::fill( dRdU_begin,              dRdU_end,              0 );
            std::fill( dRdW_begin,              dRdW_end,              0 );
            std::fill( dRdTheta_begin,          dRdTheta_end,          0 );
            std::fill( dRdE_begin,              dRdE_end,              0 );
            std::fill( dRdZ_begin,              dRdZ_end,              0 );
            std::fill( dRdVolumeFraction_begin, dRdVolumeFraction_end, 0 );
            std::fill( dRdUMesh_begin,          dRdUMesh_end,          0 );

            computeBalanceOfVolumeFraction<dim>(
                density,
                velocity_begin, velocity_end,
                volume_fraction, volume_fraction_dot, volume_fraction_gradient_begin, volume_fraction_gradient_end,
                *( material_response_begin + mass_change_rate_index ),
                rest_density,
                *( material_response_begin + trace_mass_change_velocity_gradient_index ),
                test_function, 
                interpolation_function, interpolation_function_gradient_begin, interpolation_function_gradient_end,
                dUDotdU, dVolumeFractionDotdVolumeFraction,
                result,
                *( dRdRho_begin + phase ),
                dRdU_begin + dim * phase,  dRdU_begin + dim * ( phase + 1 ),
                *( dRdVolumeFraction_begin + phase ),
                dRdC_phase,
                dRdTraceVA_phase,
                dRdUMesh_begin, dRdUMesh_end,
                volume_fraction_tolerance
            );

            // MASS CHANGE RATE
            // density
            for ( auto p = std::pair< unsigned int, dRdRho_iter >( 0, dRdRho_begin ); p.second != dRdRho_end; ++p.first, ++p.second ){

                *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + nphases * density_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * density_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // velocity
            for ( auto p = std::pair< unsigned int, dRdU_iter >( 0, dRdU_begin ); p.second != dRdU_end; ++p.first, ++p.second ){

                *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + nphases * velocity_index + p.first ) ) * interpolation_function * dUDotdU;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * velocity_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) ) * dUDotdU;

                }

            }

            // displacement
            for ( auto p = std::pair< unsigned int, dRdW_iter >( 0, dRdW_begin ); p.second != dRdW_end; ++p.first, ++p.second ){

                *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + nphases * displacement_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * displacement_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // temperature
            for ( auto p = std::pair< unsigned int, dRdTheta_iter >( 0, dRdTheta_begin ); p.second != dRdTheta_end; ++p.first, ++p.second ){

                *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + nphases * temperature_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * temperature_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // internal energy
            for ( auto p = std::pair< unsigned int, dRdE_iter >( 0, dRdE_begin ); p.second != dRdE_end; ++p.first, ++p.second ){

                *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + nphases * internal_energy_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * internal_energy_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // additional dof
            for ( auto p = std::pair< unsigned int, dRdZ_iter >( 0, dRdZ_begin ); p.second != dRdZ_end; ++p.first, ++p.second ){

                *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + nphases * additional_dof_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdC_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * additional_dof_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // mesh displacement
            for ( unsigned int I = 0; I < ( nphases * num_phase_dof + num_additional_dof ); ++I ){

                for ( unsigned int k = 0; k < material_response_dim; ++k ){

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *( dRdUMesh_begin + a ) -=
                            dRdC_phase
                            * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( mass_change_rate_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * I + k ) )
                            * ( *( full_material_response_dof_gradient_begin + material_response_dim * I + a ) )
                            * ( *( interpolation_function_gradient_begin + k ) );

                    }

                }

            }

            // TRACE MASS CHANGE VELOCITY GRADIENT
            // density
            for ( auto p = std::pair< unsigned int, dRdRho_iter >( 0, dRdRho_begin ); p.second != dRdRho_end; ++p.first, ++p.second ){

                *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + nphases * density_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * density_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // velocity
            for ( auto p = std::pair< unsigned int, dRdU_iter >( 0, dRdU_begin ); p.second != dRdU_end; ++p.first, ++p.second ){

                *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + nphases * velocity_index + p.first ) ) * interpolation_function * dUDotdU;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * velocity_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) ) * dUDotdU;

                }

            }

            // displacement
            for ( auto p = std::pair< unsigned int, dRdW_iter >( 0, dRdW_begin ); p.second != dRdW_end; ++p.first, ++p.second ){

                *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + nphases * displacement_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * displacement_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // temperature
            for ( auto p = std::pair< unsigned int, dRdTheta_iter >( 0, dRdTheta_begin ); p.second != dRdTheta_end; ++p.first, ++p.second ){

                *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + nphases * temperature_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * temperature_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // internal energy
            for ( auto p = std::pair< unsigned int, dRdE_iter >( 0, dRdE_begin ); p.second != dRdE_end; ++p.first, ++p.second ){

                *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + nphases * internal_energy_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * internal_energy_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // additional dof
            for ( auto p = std::pair< unsigned int, dRdZ_iter >( 0, dRdZ_begin ); p.second != dRdZ_end; ++p.first, ++p.second ){

                *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + nphases * additional_dof_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdTraceVA_phase * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * ( nphases * additional_dof_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // mesh displacement
            for ( unsigned int I = 0; I < ( nphases * num_phase_dof + num_additional_dof ); ++I ){

                for ( unsigned int k = 0; k < material_response_dim; ++k ){

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *( dRdUMesh_begin + a ) -=
                            dRdTraceVA_phase
                            * ( *( material_response_jacobian_begin + ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( trace_mass_change_velocity_gradient_index ) + ( nphases * num_phase_dof + num_additional_dof ) + material_response_dim * I + k ) )
                            * ( *( full_material_response_dof_gradient_begin + material_response_dim * I + a ) )
                            * ( *( interpolation_function_gradient_begin + k ) );

                    }

                }

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

        template<
            int dim, int material_response_dim, int mass_change_rate_index, int trace_mass_change_velocity_gradient_index,
            int material_response_num_dof,
            class density_iter, class velocity_iter,
            class volume_fraction_iter, class volume_fraction_dot_iter,
            class volume_fraction_gradient_iter,
            class material_response_iter,
            class material_response_jacobian_iter,
            class rest_density_iter,
            typename test_function_type,
            typename interpolation_function_type, class interpolation_function_gradient_iter,
            class full_material_response_dof_gradient_iter,
            typename dUDotdU_type, typename dVolumeFractionDotdVolumeFraction_type,
            class result_iter,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter, class dRdTheta_iter,
            class dRdE_iter, class dRdZ_iter, class dRdVolumeFraction_iter, class dRdUMesh_iter,
            int density_index        ,
            int displacement_index   ,
            int velocity_index       ,
            int temperature_index    ,
            int internal_energy_index,
            int additional_dof_index 
        >
        void computeBalanceOfVolumeFraction(
            const density_iter                  &density_begin,                  const density_iter                  &density_end,
            const velocity_iter                 &velocity_begin,                 const velocity_iter                 &velocity_end,
            const volume_fraction_iter          &volume_fraction_begin,          const volume_fraction_iter          &volume_fraction_end,
            const volume_fraction_dot_iter      &volume_fraction_dot_begin,      const volume_fraction_dot_iter      &volume_fraction_dot_end,
            const volume_fraction_gradient_iter &volume_fraction_gradient_begin, const volume_fraction_gradient_iter &volume_fraction_gradient_end,
            const material_response_iter          &material_response_begin,
            const material_response_iter          &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const rest_density_iter               &rest_density_begin,           const rest_density_iter             &rest_density_end,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU, const dVolumeFractionDotdVolumeFraction_type dVolumeFractionDotdVolumeFraction,
            result_iter   result_begin,                     result_iter            result_end,
            dRdRho_iter   dRdRho_begin,                     dRdRho_iter            dRdRho_end,
            dRdU_iter     dRdU_begin,                       dRdU_iter              dRdU_end,
            dRdW_iter     dRdW_begin,                       dRdW_iter              dRdW_end,
            dRdTheta_iter dRdTheta_begin,                   dRdTheta_iter          dRdTheta_end,
            dRdE_iter     dRdE_begin,                       dRdE_iter              dRdE_end,
            dRdZ_iter     dRdZ_begin,                       dRdZ_iter              dRdZ_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin, dRdVolumeFraction_iter dRdVolumeFraction_end,
            dRdUMesh_iter          dRdUMesh_begin,          dRdUMesh_iter          dRdUMesh_end,
            const double volume_fraction_tolerance
        ){
            /*!
             * Compute the balance of the volume fraction for a reacting continuum
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
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &rest_density_begin: The starting iterator of the rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &rest_density_end: The stopping iterator of the rest density of the material \f$ \left( \bar{\gamma}^{\alpha} \right) \f$
             * \param &test_function: The current value of the test function
             * \param &interpolation_function: The current value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the current value of the spatial gradient of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the current value of the spatial gradient of the interpolation function
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the material response dof vector
             * \param &dUDotdU: The derivative of the phase velocity w.r.t. the phase displacement dof
             * \param &dVolumeFractionDotdVolumeFraction: The derivative of the partial time derivative of the volume fraction phase w.r.t. the volume fraction
             * \param &result_begin: The starting iterator of the value of the balance of volume fraction
             * \param &result_end: The stopping iterator of the value of the balance of volume fraction
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdU_begin: The starting iterator of the derivative of the residual w.r.t. the phase spatial DOF associated with velocity
             * \param &dRdU_end: The stopping iterator of the derivative of the residual w.r.t. the phase spatial DOF associated with velocity
             * \param &dRdW_begin: The starting iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdW_end: The stopping iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdTheta_begin: The starting iterator of the derivative of the residual w.r.t. the phase temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the residual w.r.t. the phase temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the residual w.r.t. the phase internal energy
             * \param &dRdE_end: The stopping iterator of the derivative of the residual w.r.t. the phase internal energy
             * \param &dRdZ_begin: The starting iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdZ_end: The stopping iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdVolumeFraction_begin: The starting iterator of the derivative of the residual w.r.t. the volume fraction
             * \param &dRdVolumeFraction_end: The stopping iterator of the derivative of the residual w.r.t. the volume fraction
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param volume_fraction_tolerance: The tolerance of the volume fraction where if it is less than the tolerance, the
             *     true density is assumed to be the rest density.
             */

            const unsigned int nphases = ( unsigned int )( result_end - result_begin );
            const unsigned int material_response_size = ( unsigned int )( material_response_end - material_response_begin ) / nphases;
            constexpr unsigned int num_phase_dof = 3 + 2 * material_response_dim;
            const unsigned int num_additional_dof = material_response_num_dof - num_phase_dof;

            using density_type             = typename std::iterator_traits<density_iter>::value_type;
            using volume_fraction_type     = typename std::iterator_traits<volume_fraction_iter>::value_type;
            using volume_fraction_dot_type = typename std::iterator_traits<volume_fraction_dot_iter>::value_type;
            using rest_density_type        = typename std::iterator_traits<rest_density_iter>::value_type;
            using result_type              = typename std::iterator_traits<result_iter>::value_type;

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( density_end - density_begin ), "The density must have the same size as the result vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( velocity_end - velocity_begin ), "The velocity must be consistent with the size of the result vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( volume_fraction_end - volume_fraction_begin ), "The volume fraction must have the same size as the result vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( volume_fraction_dot_end - volume_fraction_dot_begin ), "The volume fraction dot must have the same size as the result vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim * nphases == ( unsigned int )( volume_fraction_gradient_end - volume_fraction_gradient_begin ), "The patial gradient of the volume fraction must be consistent with the size of the result vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size == ( unsigned int )( material_response_end - material_response_begin ), "The material response vector size must be an integer multiple of the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size * ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) == ( unsigned int )( material_response_jacobian_end - material_response_jacobian_begin ), "The material response Jacobian vector size must be consistent with the number of phases and the number of DOF in the material repsonse"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( rest_density_end - rest_density_begin ), "The rest density must have the same size as the result vector"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                ( nphases * num_phase_dof + num_additional_dof ) * material_response_dim == ( unsigned int )( full_material_response_dof_gradient_end - full_material_response_dof_gradient_begin ), "The full material response dof spatial gradient vector must be the material response dimension times the number of DOF in the material response in size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim * nphases == ( dRdU_end - dRdU_begin ), "The dRdU must be a consistent size with the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim * nphases == ( dRdW_end - dRdW_begin ), "The dRdW must be a consistent size with the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases == ( dRdTheta_end - dRdTheta_begin ), "The dRdTheta must be the same size as the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases == ( dRdE_end - dRdE_begin ), "The dRdE must be the same size as the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases == ( dRdVolumeFraction_end - dRdVolumeFraction_begin ), "The dRdVolumeFraction must be the same size as the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim == ( dRdUMesh_end - dRdUMesh_begin ), "The dRdUMesh must be the same size as the dimension"
            )

            for ( auto v = std::pair< unsigned int, result_iter >( 0, result_begin ); v.second != result_end; ++v.first, ++v.second ){

                computeBalanceOfVolumeFraction<
                    dim, material_response_dim,
                    mass_change_rate_index, trace_mass_change_velocity_gradient_index,
                    material_response_num_dof,
                    density_type, velocity_iter,
                    volume_fraction_type, volume_fraction_dot_type,
                    volume_fraction_gradient_iter,
                    material_response_iter,
                    material_response_jacobian_iter,
                    rest_density_type,
                    test_function_type,
                    interpolation_function_type, interpolation_function_gradient_iter,
                    full_material_response_dof_gradient_iter,
                    dUDotdU_type, dVolumeFractionDotdVolumeFraction_type,
                    result_type,
                    dRdRho_iter, dRdU_iter, dRdW_iter, dRdTheta_iter,
                    dRdE_iter, dRdZ_iter, dRdVolumeFraction_iter, dRdUMesh_iter,
                    density_index,
                    displacement_index,
                    velocity_index,
                    temperature_index,
                    internal_energy_index,
                    additional_dof_index
                >(
                    *( density_begin + v.first ),
                    velocity_begin + dim * v.first,                 velocity_begin + dim * ( v.first + 1 ),
                    *( volume_fraction_begin + v.first ),
                    *( volume_fraction_dot_begin + v.first ),
                    volume_fraction_gradient_begin + dim * v.first, volume_fraction_gradient_begin + dim * ( v.first + 1 ),
                    material_response_begin + material_response_size * v.first, material_response_begin + material_response_size * ( v.first + 1 ),
                    material_response_jacobian_begin + material_response_size * ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * v.first,
                    material_response_jacobian_begin + material_response_size * ( nphases * num_phase_dof + num_additional_dof ) * ( 1 + material_response_dim ) * ( v.first + 1 ),
                    *( rest_density_begin + v.first ),
                    test_function,
                    interpolation_function, interpolation_function_gradient_begin, interpolation_function_gradient_end,
                    full_material_response_dof_gradient_begin,
                    full_material_response_dof_gradient_end,
                    dUDotdU, dVolumeFractionDotdVolumeFraction,
                    v.first,
                    *v.second,
                    dRdRho_begin            +            nphases * v.first, dRdRho_begin            +            nphases * ( v.first + 1 ),
                    dRdU_begin              +      nphases * dim * v.first, dRdU_begin              +      nphases * dim * ( v.first + 1 ),
                    dRdW_begin              +      nphases * dim * v.first, dRdW_begin              +      nphases * dim * ( v.first + 1 ),
                    dRdTheta_begin          +            nphases * v.first, dRdTheta_begin          +            nphases * ( v.first + 1 ),
                    dRdE_begin              +            nphases * v.first, dRdE_begin              +            nphases * ( v.first + 1 ),
                    dRdZ_begin              + num_additional_dof * v.first, dRdZ_begin              + num_additional_dof * ( v.first + 1 ),
                    dRdVolumeFraction_begin +            nphases * v.first, dRdVolumeFraction_begin +            nphases * ( v.first + 1 ),
                    dRdUMesh_begin          +                dim * v.first, dRdUMesh_begin          +                dim * ( v.first + 1 ),
                    volume_fraction_tolerance
                );

            }

        }

    }

}
