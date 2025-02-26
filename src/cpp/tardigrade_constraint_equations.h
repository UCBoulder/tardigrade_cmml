/**
  ******************************************************************************
  * \file tardigrade_constraint_equations.h
  ******************************************************************************
  * The header file for the equations associated with the balance of mass
  ******************************************************************************
  */

#ifndef TARDIGRADE_CONSTRAINT_EQUATIONS_H
#define TARDIGRADE_CONSTRAINT_EQUATIONS_H

#define USE_EIGEN
#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace constraintEquations{

        template<
            int predicted_internal_energy_index,
            typename internal_energy_type,
            class material_response_iter,
            typename test_function_type,
            typename result_type
        >
        void computeInternalEnergyConstraint(
            const internal_energy_type &internal_energy,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const test_function_type &test_function,
            result_type &result
        );

        template<
            int predicted_internal_energy_index,
            class internal_energy_iter,
            class material_response_iter,
            typename test_function_type,
            class result_iter
        >
        void computeInternalEnergyConstraint(
            const internal_energy_iter &internal_energy_begin,
            const internal_energy_iter &internal_energy_end,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const test_function_type &test_function,
            result_iter result_begin,
            result_iter result_end
        );

    }

}

#include "tardigrade_constraint_equations.cpp"

#endif
