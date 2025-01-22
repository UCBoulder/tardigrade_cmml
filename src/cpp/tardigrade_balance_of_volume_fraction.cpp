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

    }

}
