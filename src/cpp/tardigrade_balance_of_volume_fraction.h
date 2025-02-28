/**
  ******************************************************************************
  * \file tardigrade_balance_of_volume_fraction.h
  ******************************************************************************
  * The header file for the equations associated with the balance of the
  * volume fraction
  ******************************************************************************
  */

#ifndef TARDIGRADE_BALANCE_OF_VOLUME_FRACTION_H
#define TARDIGRADE_BALANCE_OF_VOLUME_FRACTION_H

#include<array>

#define USE_EIGEN
#include "tardigrade_error_tools.h"

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
        inline void computeBalanceOfVolumeFraction(
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
            const double volume_fraction_tolerance = 1e-8
        );

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
        inline void computeBalanceOfVolumeFraction(
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
            const double volume_fraction_tolerance = 1e-8
        );

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
        inline void computeBalanceOfVolumeFraction(
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
            const double volume_fraction_tolerance = 1e-8
        );

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
        inline void computeBalanceOfVolumeFraction(
            const density_iter                  &density_begin,                  const density_iter                  &density_end,
            const velocity_iter                 &velocity_begin,                 const velocity_iter                 &velocity_end,
            const volume_fraction_iter          &volume_fraction_begin,          const volume_fraction_iter          &volume_fraction_end,
            const volume_fraction_dot_iter      &volume_fraction_dot_begin,      const volume_fraction_dot_iter      &volume_fraction_dot_end,
            const volume_fraction_gradient_iter &volume_fraction_gradient_begin, const volume_fraction_gradient_iter &volume_fraction_gradient_end,
            const material_response_iter        &material_response_begin,        const material_response_iter        &material_response_end,
            const rest_density_iter             &rest_density_begin,             const rest_density_iter             &rest_density_end,
            const test_function_type &test_function,
            result_iter result_begin, result_iter result_end,
            const double volume_fraction_tolerance = 1e-8
        );

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
        inline void computeBalanceOfVolumeFraction(
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
            const double volume_fraction_tolerance = 1e-8
        );

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
            int density_index         = 0,
            int displacement_index    = 1,
            int velocity_index        = 4,
            int temperature_index     = 7,
            int internal_energy_index = 8,
            int additional_dof_index  = 9
        >
        inline void computeBalanceOfVolumeFraction(
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
            const double volume_fraction_tolerance = 1e-8
        );

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
        inline void computeBalanceOfVolumeFraction(
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
            const double volume_fraction_tolerance = 1e-8
        );

    }

}

#include "tardigrade_balance_of_volume_fraction.cpp"

#endif
