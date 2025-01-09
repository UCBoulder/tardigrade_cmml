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

        constexpr unsigned int global_dim = 3; //!< Set the dimension as 3D by default

        constexpr unsigned int global_sot_dim = global_dim * global_dim; //!< Set the dimensions of a standard second order tensor

        typedef double floatType; //!< Define the float type as a double

        typedef std::array< floatType, global_dim > floatVector; //!< Define a standard vector

        typedef std::array< floatType, global_sot_dim > secondOrderTensor; //!< Define a standard second-order tensor

        template<
            int dim, typename density_type, typename density_dot_type,
            typename internal_energy_type, typename internal_energy_dot_type,
            typename volume_fraction_type,
            typename internal_heat_generation_type,
            typename result_type,
            class floatVector_iter, class secondOrderTensor_iter
        >
        void computeBalanceOfEnergyNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const floatVector_iter &density_gradient_begin, const floatVector_iter &density_gradient_end,
            const internal_energy_type &internal_energy, const internal_energy_dot_type &internal_energy_dot,
            const floatVector_iter &internal_energy_gradient_begin, const floatVector_iter &internal_energy_gradient_end,
            const floatVector_iter &velocity_begin, const floatVector_iter &velocity_end,
            const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
            const secondOrderTensor_iter &cauchy_stress_begin, const secondOrderTensor_iter &cauchy_stress_end,
            const volume_fraction_type &volume_fraction,
            const internal_heat_generation_type &internal_heat_generation,
            const floatVector_iter &net_interphase_force_begin, const floatVector_iter &net_interphase_force_end,
            result_type &result
        );

        template<
            int dim, typename density_type, typename density_dot_type,
            typename internal_energy_type, typename internal_energy_dot_type,
            typename volume_fraction_type,
            typename internal_heat_generation_type,
            typename result_type,
            class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out
        >
        void computeBalanceOfEnergyNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const floatVector_iter &density_gradient_begin, const floatVector_iter &density_gradient_end,
            const internal_energy_type &internal_energy, const internal_energy_type &internal_energy_dot,
            const floatVector_iter &internal_energy_gradient_begin, const floatVector_iter &internal_energy_gradient_end,
            const floatVector_iter &velocity_begin, const floatVector_iter &velocity_end,
            const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
            const secondOrderTensor_iter &cauchy_stress_begin, const secondOrderTensor_iter &cauchy_stress_end,
            const volume_fraction_type &volume_fraction,
            const internal_heat_generation_type &internal_heat_generation,
            const floatVector_iter &net_interphase_force_begin, const floatVector_iter &net_interphase_force_end,
            result_type &result,
            floatType &dRdRho, floatType &dRdRhoDot, floatVector_iter_out dRdGradRho_begin, floatVector_iter_out dRdGradRho_end,
            floatType &dRdE, floatType &dRdEDot, floatVector_iter_out dRdGradE_begin, floatVector_iter_out dRdGradE_end,
            floatVector_iter_out dRdV_begin, floatVector_iter_out dRdV_end,
            secondOrderTensor_iter_out dRdGradV_begin, secondOrderTensor_iter_out dRdGradV_end,
            secondOrderTensor_iter_out dRdCauchy_begin, secondOrderTensor_iter_out dRdCauchy_end,
            floatType &dRdPhi, floatType &dRdr,
            floatVector_iter_out dRdpi_begin, floatVector_iter_out dRdpi_end
        );

        template<
            int dim, class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class scalarArray_iter_out
        >
        void computeBalanceOfEnergyNonDivergence(
            const scalarArray_iter &density_begin, const scalarArray_iter &density_end,
            const scalarArray_iter &density_dot_begin, const scalarArray_iter &density_dot_end,
            const floatVector_iter &density_gradient_begin, const floatVector_iter &density_gradient_end,
            const scalarArray_iter &internal_energy_begin, const scalarArray_iter &internal_energy_end,
            const scalarArray_iter &internal_energy_dot_begin, const scalarArray_iter &internal_energy_dot_end,
            const floatVector_iter &internal_energy_gradient_begin, const floatVector_iter &internal_energy_gradient_end,
            const floatVector_iter &velocity_begin, const floatVector_iter &velocity_end,
            const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
            const secondOrderTensor_iter &cauchy_stress_begin, const secondOrderTensor_iter &cauchy_stress_end,
            const scalarArray_iter &volume_fraction_begin, const scalarArray_iter &volume_fraction_end,
            const scalarArray_iter &internal_heat_generation_begin, const scalarArray_iter &internal_heat_generation_end,
            const floatVector_iter &net_interphase_force_begin, const floatVector_iter &net_interphase_force_end,
            scalarArray_iter_out result_begin, scalarArray_iter_out result_end
        );

        template<
            int dim, class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class scalarArray_iter_out, class floatVector_iter_out, class secondOrderTensor_iter_out
        >
        void computeBalanceOfEnergyNonDivergence(
            const scalarArray_iter &density_begin,                  const scalarArray_iter &density_end,
            const scalarArray_iter &density_dot_begin,              const scalarArray_iter &density_dot_end,
            const floatVector_iter &density_gradient_begin,         const floatVector_iter &density_gradient_end,
            const scalarArray_iter &internal_energy_begin,          const scalarArray_iter &internal_energy_end,
            const scalarArray_iter &internal_energy_dot_begin,      const scalarArray_iter &internal_energy_dot_end,
            const floatVector_iter &internal_energy_gradient_begin, const floatVector_iter &internal_energy_gradient_end,
            const floatVector_iter &velocity_begin,                 const floatVector_iter &velocity_end,
            const secondOrderTensor_iter &velocity_gradient_begin,  const secondOrderTensor_iter &velocity_gradient_end,
            const secondOrderTensor_iter &cauchy_stress_begin,      const secondOrderTensor_iter &cauchy_stress_end,
            const scalarArray_iter &volume_fraction_begin,          const scalarArray_iter &volume_fraction_end,
            const scalarArray_iter &internal_heat_generation_begin, const scalarArray_iter &internal_heat_generation_end,
            const floatVector_iter &net_interphase_force_begin,     const floatVector_iter &net_interphase_force_end,
            scalarArray_iter_out result_begin,                      scalarArray_iter_out result_end,
            scalarArray_iter_out dRdRho_begin,                      scalarArray_iter_out dRdRho_end,
            scalarArray_iter_out dRdRhoDot_begin,                   scalarArray_iter_out dRdRhoDot_end,
            floatVector_iter_out dRdGradRho_begin,                  floatVector_iter_out dRdGradRho_end,
            scalarArray_iter_out dRdE_begin,                        scalarArray_iter_out dRdE_end,
            scalarArray_iter_out dRdEDot_begin,                     scalarArray_iter_out dRdEDot_end,
            floatVector_iter_out dRdGradE_begin,                    floatVector_iter_out dRdGradE_end,
            floatVector_iter_out dRdV_begin,                        floatVector_iter_out dRdV_end,
            secondOrderTensor_iter_out dRdGradV_begin,              secondOrderTensor_iter_out dRdGradV_end,
            secondOrderTensor_iter_out dRdCauchy_begin,             secondOrderTensor_iter_out dRdCauchy_end,
            scalarArray_iter_out dRdPhi_begin,                      scalarArray_iter_out dRdPhi_end,
            scalarArray_iter_out dRdr_begin,                        scalarArray_iter_out dRdr_end,
            floatVector_iter_out dRdpi_begin,                       floatVector_iter_out dRdpi_end
        );

        template<
            int dim,
            typename result_type, class floatVector_iter
        >
        void computeBalanceOfEnergyDivergence(
            const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
            const floatVector_iter &heat_flux_begin,              const floatVector_iter &heat_flux_end,
            result_type &result
        );

        template<
            int dim,
            typename result_type, class floatVector_iter, class floatVector_iter_out
        >
        void computeBalanceOfEnergyDivergence(
            const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
            const floatVector_iter &heat_flux_begin,              const floatVector_iter &heat_flux_end,
            result_type &result,
            floatVector_iter_out dRdGradPsi_begin,                floatVector_iter_out dRdGradPsi_end,
            floatVector_iter_out dRdq_begin,                      floatVector_iter_out dRdq_end
        );

        template<
            int dim, class floatVector_iter, class scalarArray_iter_out
        >
        void computeBalanceOfEnergyDivergence(
            const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
            const floatVector_iter &heat_flux_begin,              const floatVector_iter &heat_flux_end,
            scalarArray_iter_out result_begin,                    scalarArray_iter_out result_end
        );

        template<
            int dim, class floatVector_iter, class scalarArray_iter_out, class floatVector_iter_out
        >
        void computeBalanceOfEnergyDivergence(
            const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
            const floatVector_iter &heat_flux_begin,              const floatVector_iter &heat_flux_end,
            scalarArray_iter_out result_begin,                    scalarArray_iter_out result_end,
            floatVector_iter_out dRdGradPsi_begin,                floatVector_iter_out dRdGradPsi_end,
            floatVector_iter_out dRdq_begin,                      floatVector_iter_out dRdq_end
        );

    }

}

#include "tardigrade_balance_of_energy.cpp"

#endif
