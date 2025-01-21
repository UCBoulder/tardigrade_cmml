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

        template<
            int dim,
            typename density_type, typename density_dot_type, class density_gradient_iter,
            class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter,
            class result_iter
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
            int dim,
            typename density_type, typename density_dot_type, class density_gradient_iter,
            class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter,
            class cauchy_stress_iter, typename volume_fraction_type,
            typename testFunction_type, class testFunctionGradient_iter,
            class result_iter
        >
        void computeBalanceOfLinearMomentum(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            const cauchy_stress_iter &cauchy_stress_begin,         const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_type &volume_fraction,
            const testFunction_type &test_function,
            const testFunctionGradient_iter &test_function_gradient_begin,
            const testFunctionGradient_iter &test_function_gradient_end,
            result_iter result_begin,                              result_iter result_end
        );

        template<
            int dim,
            typename density_type, typename density_dot_type, class density_gradient_iter,
            class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter,
            class cauchy_stress_iter, typename volume_fraction_type,
            typename testFunction_type, class testFunctionGradient_iter,
            typename interpolationFunction_type, class interpolationFunctionGradient_iter,
            class result_iter,
            class dRdRho_iter, class dRdU_iter, class dRdB_iter, class dRdCauchy_iter,
            class dRdVolumeFraction_iter, class dRdUMesh_iter,
            typename dRhoDotdRho_type, typename dUDotdU_type, typename dUDDotdU_type
        >
        void computeBalanceOfLinearMomentum(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            const cauchy_stress_iter &cauchy_stress_begin,         const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_type &volume_fraction,
            const testFunction_type &test_function,
            const testFunctionGradient_iter &test_function_gradient_begin,
            const testFunctionGradient_iter &test_function_gradient_end,
            const interpolationFunction_type &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const dRhoDotdRho_type dRhoDotdRho,
            const dUDotdU_type dUDotdU, const dUDDotdU_type dUDDotdU,
            result_iter result_begin,                              result_iter result_end,
            dRdRho_iter dRdRho_begin,                              dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin,                                  dRdU_iter dRdU_end,
            dRdB_iter dRdB_begin,                                  dRdB_iter dRdB_end,
            dRdCauchy_iter dRdCauchy_begin,                        dRdCauchy_iter dRdCauchy_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin,        dRdVolumeFraction_iter dRdVolumeFraction_end,
            dRdUMesh_iter dRdUMesh_begin,                          dRdUMesh_iter dRdUMesh_end
        );

        template<
            int dim,
            class density_iter, class density_dot_iter, class density_gradient_iter,
            class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter,
            class cauchy_stress_iter, class volume_fraction_iter,
            typename testFunction_type, class testFunctionGradient_iter,
            class result_iter
        >
        void computeBalanceOfLinearMomentum(
            const density_iter &density_begin,                     const density_iter &density_end,
            const density_dot_iter &density_dot_begin,             const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            const cauchy_stress_iter &cauchy_stress_begin,         const cauchy_stress_iter &cauchy_stress_end,
            const volume_fraction_iter &volume_fraction_begin,     const volume_fraction_iter &volume_fraction_end,
            const testFunction_type &test_function,
            const testFunctionGradient_iter &test_function_gradient_begin,
            const testFunctionGradient_iter &test_function_gradient_end,
            result_iter result_begin,                              result_iter result_end
        );

        template<
            int dim,
            typename density_type, typename density_dot_type, class density_gradient_iter,
            class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter,
            class result_iter,
            class dRdRho_iter, class dRdRhoDot_iter, class dRdGradRho_iter,
            class dRdV_iter, class dRdVDot_iter, class dRdGradV_iter,
            class dRdB_iter
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
            dRdGradRho_iter dRdGradRho_begin,                      dRdGradRho_iter dRdGradRho_end,
            dRdV_iter dRdV_begin,                                  dRdV_iter dRdV_end,
            dRdVDot_iter dRdVDot_begin,                            dRdVDot_iter dRdVDot_end,
            dRdGradV_iter dRdGradV_begin,                          dRdGradV_iter dRdGradV_end,
            dRdB_iter dRdB_begin,                                  dRdB_iter dRdB_end
        );

        template<
            int dim,
            class density_iter, class density_dot_iter, class density_gradient_iter,
            class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter,
            class result_iter
        >
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_iter &density_begin,                     const density_iter &density_end,
            const density_dot_iter &density_dot_begin,             const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            result_iter result_begin,                              result_iter result_end
        );

        template<
            int dim,
            class density_iter, class density_dot_iter, class density_gradient_iter,
            class velocity_iter, class velocity_dot_iter, class velocity_gradient_iter,
            class body_force_iter,
            class result_iter,
            class dRdRho_iter, class dRdRhoDot_iter, class dRdGradRho_iter,
            class dRdV_iter, class dRdVDot_iter, class dRdGradV_iter,
            class dRdB_iter
        >
        void computeBalanceOfLinearMomentumNonDivergence(
            const density_iter &density_begin,                     const density_iter &density_end,
            const density_dot_iter &density_dot_begin,             const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin,   const density_gradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocity_dot_iter &velocity_dot_begin,           const velocity_dot_iter &velocity_dot_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const body_force_iter &body_force_begin,               const body_force_iter &body_force_end,
            result_iter result_begin,                              result_iter result_end,
            dRdRho_iter dRdRho_begin,                              dRdRho_iter dRdRho_end,
            dRdRhoDot_iter dRdRhoDot_begin,                        dRdRhoDot_iter dRdRhoDot_end,
            dRdGradRho_iter dRdGradRho_begin,                      dRdGradRho_iter dRdGradRho_end,
            dRdV_iter dRdV_begin,                                  dRdV_iter dRdV_end,
            dRdVDot_iter dRdVDot_begin,                            dRdVDot_iter dRdVDot_end,
            dRdGradV_iter dRdGradV_begin,                          dRdGradV_iter dRdGradV_end,
            dRdB_iter dRdB_begin,                                  dRdB_iter dRdB_end
        );

        template<
            int dim,
            class testFunctionGradient_iter, class cauchyStress_iter, typename volumeFraction_type,
            class result_iter
        >
        void computeBalanceOfLinearMomentumDivergence(
            const testFunctionGradient_iter &grad_test_function_begin,
            const testFunctionGradient_iter &grad_test_function_end,
            const cauchyStress_iter &cauchy_stress_begin,    const cauchyStress_iter &cauchy_stress_end,
            const volumeFraction_type &volume_fraction,
            result_iter result_begin,                        result_iter result_end
        );

        template<
            int dim,
            class testFunctionGradient_iter, class cauchyStress_iter, typename volumeFraction_type,
            class result_iter,
            class dRdGradPsi_iter, class dRdCauchy_iter, class dRdVolumeFraction_iter
        >
        void computeBalanceOfLinearMomentumDivergence(
            const testFunctionGradient_iter &grad_test_function_begin,
            const testFunctionGradient_iter &grad_test_function_end,
            const cauchyStress_iter &cauchy_stress_begin,    const cauchyStress_iter &cauchy_stress_end,
            const volumeFraction_type &volume_fraction,
            result_iter result_begin,                        result_iter result_end,
            dRdGradPsi_iter dRdGradPsi_begin,                dRdGradPsi_iter dRdGradPsi_end,
            dRdCauchy_iter dRdCauchy_begin,                  dRdCauchy_iter dRdCauchy_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin,  dRdVolumeFraction_iter dRdVolumeFraction_end
        );

        template<
            int dim,
            class testFunctionGradient_iter, class cauchyStress_iter, class volumeFraction_iter,
            class result_iter
        >
        void computeBalanceOfLinearMomentumDivergence(
            const testFunctionGradient_iter &grad_test_function_begin,
            const testFunctionGradient_iter &grad_test_function_end,
            const cauchyStress_iter &cauchy_stress_begin,     const cauchyStress_iter &cauchy_stress_end,
            const volumeFraction_iter &volume_fraction_begin, const volumeFraction_iter &volume_fraction_end,
            result_iter result_begin,                          result_iter result_end
        );

        template<
            int dim,
            class testFunctionGradient_iter, class cauchyStress_iter, class volumeFraction_iter,
            class result_iter,
            class dRdGradPsi_iter, class dRdCauchy_iter, class dRdVolumeFraction_iter
        >
        void computeBalanceOfLinearMomentumDivergence(
            const testFunctionGradient_iter &grad_test_function_begin,
            const testFunctionGradient_iter &grad_test_function_end,
            const cauchyStress_iter &cauchy_stress_begin,     const cauchyStress_iter &cauchy_stress_end,
            const volumeFraction_iter &volume_fraction_begin, const volumeFraction_iter &volume_fraction_end,
            result_iter result_begin,                         result_iter result_end,
            dRdGradPsi_iter dRdGradPsi_begin,                 dRdGradPsi_iter dRdGradPsi_end,
            dRdCauchy_iter dRdCauchy_begin,                   dRdCauchy_iter dRdCauchy_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin,   dRdVolumeFraction_iter dRdVolumeFraction_end
        );

    }

}

#include "tardigrade_balance_of_linear_momentum.cpp"

#endif
