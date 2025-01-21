/**
  ******************************************************************************
  * \file tardigrade_balance_of_energy.h
  ******************************************************************************
  * The header file for the equations associated with the balance of energy
  ******************************************************************************
  */

#ifndef TARDIGRADE_BALANCE_OF_ENERGY_H
#define TARDIGRADE_BALANCE_OF_ENERGY_H

#include<array>

#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace balanceOfEnergy{

        template<
            int dim,
            typename density_type, typename density_dot_type,
            class density_gradient_iter,
            typename internal_energy_type, typename internal_energy_dot_type,
            class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class cauchy_stress_iter,
            typename volume_fraction_type,
            typename internal_heat_generation_type,
            class net_interphase_force_iter,
            class heat_flux_iter,
            typename test_function_type, class test_function_gradient_iter,
            typename result_type
        >
        void computeBalanceOfEnergy(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const internal_energy_type &internal_energy, const internal_energy_dot_type &internal_energy_dot,
            const internal_energy_gradient_iter &internal_energy_gradient_begin, const internal_energy_gradient_iter &internal_energy_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_type &volume_fraction,
            const internal_heat_generation_type &internal_heat_generation,
            const net_interphase_force_iter &net_interphase_force_begin, const net_interphase_force_iter &net_interphase_force_end,
            const heat_flux_iter &heat_flux_begin, const heat_flux_iter &heat_flux_end,
            const test_function_type &test_function,
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            result_type &result
        );

        template<
            int dim,
            typename density_type, typename density_dot_type,
            class density_gradient_iter,
            typename internal_energy_type, typename internal_energy_dot_type,
            class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class cauchy_stress_iter,
            typename volume_fraction_type,
            typename internal_heat_generation_type,
            class net_interphase_force_iter,
            class heat_flux_iter,
            typename test_function_type, class test_function_gradient_iter,
            typename interpolation_function_type, class interpolation_function_gradient_iter,
            typename dRhoDotdRho_type, typename dEDotdE_type, typename dUDotdU_type,
            typename result_type,
            typename dRdRho_type, typename dRdE_type,
            class dRdU_iter,
            class dRdCauchy_iter, typename dRdVolumeFraction_type,
            typename dRdr_type, class dRdpi_iter, class dRdq_iter,
            class dRdUMesh_iter
        >
        void computeBalanceOfEnergy(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const internal_energy_type &internal_energy, const internal_energy_dot_type &internal_energy_dot,
            const internal_energy_gradient_iter &internal_energy_gradient_begin, const internal_energy_gradient_iter &internal_energy_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_type &volume_fraction,
            const internal_heat_generation_type &internal_heat_generation,
            const net_interphase_force_iter &net_interphase_force_begin, const net_interphase_force_iter &net_interphase_force_end,
            const heat_flux_iter &heat_flux_begin, const heat_flux_iter &heat_flux_end,
            const test_function_type &test_function,
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin, const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const dRhoDotdRho_type dRhoDotdRho, const dEDotdE_type dEDotdE, const dUDotdU_type dUDotdU,
            result_type &result,
            dRdRho_type &dRdRho, dRdE_type &dRdE,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end,
            dRdCauchy_iter dRdCauchy_begin, dRdCauchy_iter dRdCauchy_end,
            dRdVolumeFraction_type &dRdVolumeFraction,
            dRdr_type &dRdr,
            dRdpi_iter dRdpi_begin, dRdpi_iter dRdpi_end,
            dRdq_iter dRdq_begin, dRdq_iter dRdq_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end
        );

        template<
            int dim,
            class density_iter, class density_dot_iter,
            class density_gradient_iter,
            class internal_energy_iter, class internal_energy_dot_iter,
            class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class cauchy_stress_iter,
            class volume_fraction_iter,
            class internal_heat_generation_iter,
            class net_interphase_force_iter,
            class heat_flux_iter,
            typename test_function_type, class test_function_gradient_iter,
            class result_iter
        >
        void computeBalanceOfEnergy(
            const density_iter &density_begin,                                   const density_iter &density_end,
            const density_dot_iter &density_dot_begin,                           const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin,                 const density_gradient_iter &density_gradient_end,
            const internal_energy_iter &internal_energy_begin,                   const internal_energy_iter &internal_energy_end,
            const internal_energy_dot_iter &internal_energy_dot_begin,           const internal_energy_dot_iter &internal_energy_dot_end,
            const internal_energy_gradient_iter &internal_energy_gradient_begin, const internal_energy_gradient_iter &internal_energy_gradient_end,
            const velocity_iter &velocity_begin,                                 const velocity_iter &velocity_end,
            const velocity_gradient_iter &velocity_gradient_begin,               const velocity_gradient_iter &velocity_gradient_end,
            const cauchy_stress_iter &cauchy_stress_begin,                       const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_iter &volume_fraction_begin,                   const volume_fraction_iter &volume_fraction_end,
            const internal_heat_generation_iter &internal_heat_generation_begin, const internal_heat_generation_iter &internal_heat_generation_end,
            const net_interphase_force_iter &net_interphase_force_begin,         const net_interphase_force_iter &net_interphase_force_end,
            const heat_flux_iter &heat_flux_begin,                               const heat_flux_iter &heat_flux_end,
            const test_function_type &test_function,
            const test_function_gradient_iter &test_function_gradient_begin,     const test_function_gradient_iter &test_function_gradient_end,
            result_iter result_begin,                                            result_iter result_end
        );

        template<
            int dim,
            class density_iter, class density_dot_iter,
            class density_gradient_iter,
            class internal_energy_iter, class internal_energy_dot_iter,
            class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class cauchy_stress_iter,
            class volume_fraction_iter,
            class internal_heat_generation_iter,
            class net_interphase_force_iter,
            class heat_flux_iter,
            typename test_function_type, class test_function_gradient_iter,
            typename interpolation_function_type, class interpolation_function_gradient_iter,
            typename dRhoDotdRho_type, typename dEDotdE_type, typename dUDotdU_type,
            class result_iter,
            class dRdRho_iter, class dRdE_iter,
            class dRdU_iter,
            class dRdCauchy_iter, class dRdVolumeFraction_iter,
            class dRdr_iter, class dRdpi_iter, class dRdq_iter,
            class dRdUMesh_iter
        >
        void computeBalanceOfEnergy(
            const density_iter &density_begin,                                   const density_iter &density_end,
            const density_dot_iter &density_dot_begin,                           const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin,                 const density_gradient_iter &density_gradient_end,
            const internal_energy_iter &internal_energy_begin,                   const internal_energy_iter &internal_energy_end,
            const internal_energy_dot_iter &internal_energy_dot_begin,           const internal_energy_dot_iter &internal_energy_dot_end,
            const internal_energy_gradient_iter &internal_energy_gradient_begin, const internal_energy_gradient_iter &internal_energy_gradient_end,
            const velocity_iter &velocity_begin,                                 const velocity_iter &velocity_end,
            const velocity_gradient_iter &velocity_gradient_begin,               const velocity_gradient_iter &velocity_gradient_end,
            const cauchy_stress_iter &cauchy_stress_begin,                       const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_iter &volume_fraction_begin,                   const volume_fraction_iter &volume_fraction_end,
            const internal_heat_generation_iter &internal_heat_generation_begin, const internal_heat_generation_iter &internal_heat_generation_end,
            const net_interphase_force_iter &net_interphase_force_begin,         const net_interphase_force_iter &net_interphase_force_end,
            const heat_flux_iter &heat_flux_begin,                               const heat_flux_iter &heat_flux_end,
            const test_function_type &test_function,
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin, const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const dRhoDotdRho_type dRhoDotdRho, const dEDotdE_type dEDotdE, const dUDotdU_type dUDotdU,
            result_iter result_begin,                       result_iter result_end,
            dRdRho_iter dRdRho_begin,                       dRdRho_iter dRdRho_end,
            dRdE_iter dRdE_begin,                           dRdE_iter dRdE_end,
            dRdU_iter dRdU_begin,                           dRdU_iter dRdU_end,
            dRdCauchy_iter dRdCauchy_begin,                 dRdCauchy_iter dRdCauchy_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin, dRdVolumeFraction_iter dRdVolumeFraction_end,
            dRdr_iter dRdr_begin,                           dRdr_iter dRdr_end,
            dRdpi_iter dRdpi_begin,                         dRdpi_iter dRdpi_end,
            dRdq_iter dRdq_begin,                           dRdq_iter dRdq_end,
            dRdUMesh_iter dRdUMesh_begin,                   dRdUMesh_iter dRdUMesh_end
        );

        template<
            int dim,
            typename density_type, typename density_dot_type,
            class density_gradient_iter,
            typename internal_energy_type, typename internal_energy_dot_type,
            class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class cauchy_stress_iter,
            typename volume_fraction_type,
            typename internal_heat_generation_type,
            class net_interphase_force_iter,
            typename result_type
        >
        void computeBalanceOfEnergyNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const internal_energy_type &internal_energy, const internal_energy_dot_type &internal_energy_dot,
            const internal_energy_gradient_iter &internal_energy_gradient_begin, const internal_energy_gradient_iter &internal_energy_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_type &volume_fraction,
            const internal_heat_generation_type &internal_heat_generation,
            const net_interphase_force_iter &net_interphase_force_begin, const net_interphase_force_iter &net_interphase_force_end,
            result_type &result
        );

        template<
            int dim,
            typename density_type, typename density_dot_type,
            class density_gradient_iter,
            typename internal_energy_type, typename internal_energy_dot_type,
            class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class cauchy_stress_iter,
            typename volume_fraction_type,
            typename internal_heat_generation_type,
            class net_interphase_force_iter,
            typename result_type,
            typename dRdRho_type, typename dRdRhoDot_type,
            class dRdGradRho_iter,
            typename dRdE_type, typename dRdEDot_type, class dRdGradE_iter,
            class dRdV_iter, class dRdGradV_iter,
            class dRdCauchy_iter,
            typename dRdVolumeFraction_type, typename dRdr_type,
            class dRdpi_iter
        >
        void computeBalanceOfEnergyNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const internal_energy_type &internal_energy, const internal_energy_dot_type &internal_energy_dot,
            const internal_energy_gradient_iter &internal_energy_gradient_begin, const internal_energy_gradient_iter &internal_energy_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_type &volume_fraction,
            const internal_heat_generation_type &internal_heat_generation,
            const net_interphase_force_iter &net_interphase_force_begin, const net_interphase_force_iter &net_interphase_force_end,
            result_type &result,
            dRdRho_type &dRdRho, dRdRhoDot_type &dRdRhoDot,
            dRdGradRho_iter dRdGradRho_begin, dRdGradRho_iter dRdGradRho_end,
            dRdE_type &dRdE, dRdEDot_type &dRdEDot, dRdGradE_iter dRdGradE_begin, dRdGradE_iter dRdGradE_end,
            dRdV_iter dRdV_begin, dRdV_iter dRdV_end,
            dRdGradV_iter dRdGradV_begin, dRdGradV_iter dRdGradV_end,
            dRdCauchy_iter dRdCauchy_begin, dRdCauchy_iter dRdCauchy_end,
            dRdVolumeFraction_type &dRdVolumeFraction, dRdr_type &dRdr,
            dRdpi_iter dRdpi_begin, dRdpi_iter dRdpi_end
        );

        template<
            int dim,
            class density_iter, class density_dot_iter, class density_gradient_iter,
            class internal_energy_iter, class internal_energy_dot_iter, class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class cauchy_stress_iter,
            class volume_fraction_iter,
            class internal_heat_generation_iter,
            class net_interphase_force_iter,
            class result_iter
        >
        void computeBalanceOfEnergyNonDivergence(
            const density_iter &density_begin, const density_iter &density_end,
            const density_dot_iter &density_dot_begin, const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const internal_energy_iter &internal_energy_begin, const internal_energy_iter &internal_energy_end,
            const internal_energy_dot_iter &internal_energy_dot_begin, const internal_energy_dot_iter &internal_energy_dot_end,
            const internal_energy_gradient_iter &internal_energy_gradient_begin, const internal_energy_gradient_iter &internal_energy_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const cauchy_stress_iter &cauchy_stress_begin, const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_iter &volume_fraction_begin, const volume_fraction_iter &volume_fraction_end,
            const internal_heat_generation_iter &internal_heat_generation_begin, const internal_heat_generation_iter &internal_heat_generation_end,
            const net_interphase_force_iter &net_interphase_force_begin, const net_interphase_force_iter &net_interphase_force_end,
            result_iter result_begin, result_iter result_end
        );

        template<
            int dim,
            class density_iter, class density_dot_iter, class density_gradient_iter,
            class internal_energy_iter, class internal_energy_dot_iter, class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class cauchy_stress_iter,
            class volume_fraction_iter,
            class internal_heat_generation_iter,
            class net_interphase_force_iter,
            class result_iter,
            class dRdRho_iter, class dRdRhoDot_iter, class dRdGradRho_iter,
            class dRdE_iter, class dRdEDot_iter, class dRdGradE_iter,
            class dRdV_iter, class dRdGradV_iter,
            class dRdCauchy_iter,
            class dRdVolumeFraction_iter, class dRdr_iter,
            class dRdpi_iter
        >
        void computeBalanceOfEnergyNonDivergence(
            const density_iter &density_begin,                      const density_iter &density_end,
            const density_dot_iter &density_dot_begin,              const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin,    const density_gradient_iter &density_gradient_end,
            const internal_energy_iter &internal_energy_begin,      const internal_energy_iter &internal_energy_end,
            const internal_energy_dot_iter &internal_energy_dot_begin, const internal_energy_dot_iter &internal_energy_dot_end,
            const internal_energy_gradient_iter &internal_energy_gradient_begin, const internal_energy_gradient_iter &internal_energy_gradient_end,
            const velocity_iter &velocity_begin,                    const velocity_iter &velocity_end,
            const velocity_gradient_iter &velocity_gradient_begin,  const velocity_gradient_iter &velocity_gradient_end,
            const cauchy_stress_iter &cauchy_stress_begin,          const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_iter &volume_fraction_begin,      const volume_fraction_iter &volume_fraction_end,
            const internal_heat_generation_iter &internal_heat_generation_begin, const internal_heat_generation_iter &internal_heat_generation_end,
            const net_interphase_force_iter &net_interphase_force_begin, const net_interphase_force_iter &net_interphase_force_end,
            result_iter result_begin,                               result_iter result_end,
            dRdRho_iter dRdRho_begin,                               dRdRho_iter dRdRho_end,
            dRdRhoDot_iter dRdRhoDot_begin,                         dRdRhoDot_iter dRdRhoDot_end,
            dRdGradRho_iter dRdGradRho_begin,                       dRdGradRho_iter dRdGradRho_end,
            dRdE_iter dRdE_begin,                                   dRdE_iter dRdE_end,
            dRdEDot_iter dRdEDot_begin,                             dRdEDot_iter dRdEDot_end,
            dRdGradE_iter dRdGradE_begin,                           dRdGradE_iter dRdGradE_end,
            dRdV_iter dRdV_begin,                                   dRdV_iter dRdV_end,
            dRdGradV_iter dRdGradV_begin,                           dRdGradV_iter dRdGradV_end,
            dRdCauchy_iter dRdCauchy_begin,                         dRdCauchy_iter dRdCauchy_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin,         dRdVolumeFraction_iter dRdVolumeFraction_end,
            dRdr_iter dRdr_begin,                                   dRdr_iter dRdr_end,
            dRdpi_iter dRdpi_begin,                                 dRdpi_iter dRdpi_end
        );

        template<
            int dim,
            class test_function_gradient_iter, class heat_flux_iter,
            typename result_type
        >
        void computeBalanceOfEnergyDivergence(
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            const heat_flux_iter &heat_flux_begin,                           const heat_flux_iter &heat_flux_end,
            result_type &result
        );

        template<
            int dim,
            class test_function_gradient_iter, class heat_flux_iter,
            typename result_type,
            class dRdGradTestFunction_iter, class dRdq_iter
        >
        void computeBalanceOfEnergyDivergence(
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            const heat_flux_iter &heat_flux_begin,                           const heat_flux_iter &heat_flux_end,
            result_type &result,
            dRdGradTestFunction_iter dRdGradTestFunction_begin,              dRdGradTestFunction_iter dRdGradTestFunction_end,
            dRdq_iter dRdq_begin,                                            dRdq_iter dRdq_end
        );

        template<
            int dim,
            class test_function_gradient_iter, class heat_flux_iter,
            class result_iter
        >
        void computeBalanceOfEnergyDivergence(
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            const heat_flux_iter &heat_flux_begin,                           const heat_flux_iter &heat_flux_end,
            result_iter result_begin,                                        result_iter result_end
        );

        template<
            int dim,
            class test_function_gradient_iter, class heat_flux_iter,
            class result_iter,
            class dRdGradTestFunction_iter, class dRdq_iter
        >
        void computeBalanceOfEnergyDivergence(
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            const heat_flux_iter &heat_flux_begin,                           const heat_flux_iter &heat_flux_end,
            result_iter result_begin,              result_iter result_end,
            dRdGradTestFunction_iter dRdGradTestFunction_begin,      dRdGradTestFunction_iter dRdGradTestFunction_end,
            dRdq_iter dRdq_begin,                  dRdq_iter dRdq_end
        );

    }

}

#include "tardigrade_balance_of_energy.cpp"

#endif
