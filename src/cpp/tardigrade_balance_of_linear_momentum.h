/**
  ******************************************************************************
  * \file tardigrade_balance_of_linear_momentum.h
  ******************************************************************************
  * The header file for the equations associated with the balance of linear
  * momentum
  ******************************************************************************
  */

#ifndef TARDIGRADE_BALANCE_OF_LINEAR_MOMENTUM_H
#define TARDIGRADE_BALANCE_OF_LINEAR_MOMENTUM_H

#include<array>

#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace balanceOfLinearMomentum{

        constexpr unsigned int global_dim = 3; //!< Set the dimension as 3D by default

        constexpr unsigned int global_sot_dim = global_dim * global_dim; //!< Set the dimensions of a standard second order tensor

        typedef double floatType; //!< Define the float type as a double

        typedef std::array< floatType, global_dim > floatVector; //!< Define a standard vector

        typedef std::array< floatType, global_sot_dim > secondOrderTensor; //!< Define a standard second-order tensor

        template<
            int dim, typename density_type, typename density_dot_type,
            class density_gradient_iter, class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter, class result_iter
        >
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            result_iter result_begin,                              result_iter result_end
        );

        template<
            int dim, typename density_type, typename density_dot_type, typename testFunction_type,
            class density_gradient_iter, class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter, class result_iter, class testFunctionGradient_iter
        >
        void computeBalanceOfLinearMomentum(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            const testFunction_type &psi,
            const testFunctionGradient_iter &psi_gradient_begin,   const testFunctionGradient_iter &psi_gradient_end,
            result_iter result_begin,                              result_iter result_end
        );

        template<
            int dim, typename density_type, typename density_dot_type,
            class density_gradient_iter, class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter, class result_iter, class dRdRho_iter, class dRdRhoDot_iter, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out
        >
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            result_iter result_begin,                              result_iter result_end,
            dRdRho_iter dRdRho_begin,                              dRdRho_iter dRdRho_end,
            dRdRhoDot_iter dRdRhoDot_begin,                        dRdRhoDot_iter dRdRhoDot_end,
            secondOrderTensor_iter_out dRdGradRho_begin,           secondOrderTensor_iter_out dRdGradRho_end,
            secondOrderTensor_iter_out dRdV_begin,                 secondOrderTensor_iter_out dRdV_end,
            secondOrderTensor_iter_out dRdVDot_begin,              secondOrderTensor_iter_out dRdVDot_end,
            thirdOrderTensor_iter_out dRdGradV_begin,              thirdOrderTensor_iter_out dRdGradV_end,
            secondOrderTensor_iter_out dRdB_begin,                 secondOrderTensor_iter_out dRdB_end
        );

        template<
            int dim, class density_iter, class density_dot_iter,
            class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out
        >
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_iter &density_begin,                     const density_iter &density_end,
            const density_dot_iter &density_dot_begin,             const density_dot_iter &density_dot_end,
            const floatVector_iter &density_gradient_begin,        const floatVector_iter &density_gradient_end,
            const floatVector_iter &velocity_begin,                const floatVector_iter &velocity_end,
            const floatVector_iter &velocity_dot_begin,            const floatVector_iter &velocity_dot_end,
            const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
            const floatVector_iter &body_force_begin,              const floatVector_iter &body_force_end,
            floatVector_iter_out result_begin,                     floatVector_iter_out result_end
        );

        template<
            int dim, class density_iter, class density_dot_iter,
            class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out
        >
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_iter &density_begin,                     const density_iter &density_end,
            const density_dot_iter &density_dot_begin,             const density_dot_iter &density_dot_end,
            const floatVector_iter &density_gradient_begin,        const floatVector_iter &density_gradient_end,
            const floatVector_iter &velocity_begin,                const floatVector_iter &velocity_end,
            const floatVector_iter &velocity_dot_begin,            const floatVector_iter &velocity_dot_end,
            const secondOrderTensor_iter &velocity_gradient_begin, const secondOrderTensor_iter &velocity_gradient_end,
            const floatVector_iter &body_force_begin,              const floatVector_iter &body_force_end,
            floatVector_iter_out result_begin,                     floatVector_iter_out result_end,
            floatVector_iter_out dRdRho_begin,                     floatVector_iter_out dRdRho_end,
            floatVector_iter_out dRdRhoDot_begin,                  floatVector_iter_out dRdRhoDot_end,
            secondOrderTensor_iter_out dRdGradRho_begin,           secondOrderTensor_iter_out dRdGradRho_end,
            secondOrderTensor_iter_out dRdV_begin,                 secondOrderTensor_iter_out dRdV_end,
            secondOrderTensor_iter_out dRdVDot_begin,              secondOrderTensor_iter_out dRdVDot_end,
            thirdOrderTensor_iter_out dRdGradV_begin,              thirdOrderTensor_iter_out dRdGradV_end,
            secondOrderTensor_iter_out dRdB_begin,                 secondOrderTensor_iter_out dRdB_end
        );

        template<int dim, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out>
        void computeBalanceOfLinearMomentumDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                                       const secondOrderTensor_iter &cauchy_stress_begin,    const secondOrderTensor_iter &cauchy_stress_end,
                                                       const floatType &volume_fraction,
                                                       floatVector_iter_out result_begin,                   floatVector_iter_out result_end );

        template<int dim, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out>
        void computeBalanceOfLinearMomentumDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                                       const secondOrderTensor_iter &cauchy_stress_begin,    const secondOrderTensor_iter &cauchy_stress_end,
                                                       const floatType &volume_fraction,
                                                       floatVector_iter_out result_begin,                    floatVector_iter_out result_end,
                                                       secondOrderTensor_iter_out dRdGradPsi_begin,          secondOrderTensor_iter_out dRdGradPsi_end,
                                                       thirdOrderTensor_iter_out dRdCauchy_begin,            thirdOrderTensor_iter_out dRdCauchy_end,
                                                       floatVector_iter_out dRdPhi_begin,                    floatVector_iter_out dRdPhi_end );

        template<int dim, class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out>
        void computeBalanceOfLinearMomentumDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                                       const secondOrderTensor_iter &cauchy_stress_begin,    const secondOrderTensor_iter &cauchy_stress_end,
                                                       const scalarArray_iter &volume_fraction_begin,        const scalarArray_iter &volume_fraction_end,
                                                       floatVector_iter_out result_begin,                    floatVector_iter_out result_end );

        template<int dim, class scalarArray_iter, class floatVector_iter, class secondOrderTensor_iter, class floatVector_iter_out, class secondOrderTensor_iter_out, class thirdOrderTensor_iter_out>
        void computeBalanceOfLinearMomentumDivergence( const floatVector_iter &test_function_gradient_begin, const floatVector_iter &test_function_gradient_end,
                                                       const secondOrderTensor_iter &cauchy_stress_begin,    const secondOrderTensor_iter &cauchy_stress_end,
                                                       const scalarArray_iter &volume_fraction_begin,        const scalarArray_iter &volume_fraction_end,
                                                       floatVector_iter_out result_begin,                    floatVector_iter_out result_end,
                                                       secondOrderTensor_iter_out dRdGradPsi_begin,          secondOrderTensor_iter_out dRdGradPsi_end,
                                                       thirdOrderTensor_iter_out dRdCauchy_begin,            thirdOrderTensor_iter_out dRdCauchy_end,
                                                       floatVector_iter_out dRdPhi_begin,                    floatVector_iter_out dRdPhi_end );

    }

}

#include "tardigrade_balance_of_linear_momentum.cpp"

#endif
