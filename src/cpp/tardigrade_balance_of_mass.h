/**
  ******************************************************************************
  * \file tardigrade_balance_of_mass.h
  ******************************************************************************
  * The header file for the equations associated with the balance of mass
  ******************************************************************************
  */

#ifndef TARDIGRADE_BALANCE_OF_MASS_H
#define TARDIGRADE_BALANCE_OF_MASS_H

#include<array>

#define USE_EIGEN
#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace balanceOfMass{

        constexpr unsigned int global_dim = 3; //!< Set the dimension as 3D by default

        constexpr unsigned int global_sot_dim = global_dim * global_dim; //!< Set the dimensions of a standard second order tensor

        typedef double floatType; //!< Define the float type as a double

        typedef std::array< floatType, global_dim > floatVector; //!< Define a standard vector

        typedef std::array< floatType, global_sot_dim > secondOrderTensor; //!< Define a standard second-order tensor

        void computeBalanceOfMass(
            const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
            const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &result
        );

        void computeBalanceOfMass(
            const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
            const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &result,
            floatType   &dRdRho, floatType         &dRdRhoDot, floatVector &dRdGradRho,
            floatVector &dRdV,   secondOrderTensor &dRdGradV
        );

        template<
            int dim, typename density_type, typename densityDot_type, typename result_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin, const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            result_type &result
        );

        template<
            int dim, typename density_type, typename densityDot_type, typename result_type,
            typename testFunction_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin, const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &psi,
            result_type &result
        );

        template<
            int dim, int matresponse_dim, typename density_type, typename densityDot_type, typename result_type,
            typename testFunction_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter,
            class material_response_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin, const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin,      const material_response_iter &material_response_end,
            const testFunction_type &psi,
            result_type &result
        );

        template<
            int dim, int matresponse_dim, int matresponse_values, class density_iter, class densityDot_iter, class result_iter,
            typename testFunction_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter,
            class material_response_iter
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                     const density_iter &density_end,
            const densityDot_iter &density_dot_begin,              const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,    const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                   const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,  const velocityGradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const testFunction_type &test_function,
            result_iter result_begin, result_iter result_end
        );

        template<
            int dim, typename density_type, typename densityDot_type, typename result_type,
            typename dRdRho_type, typename dRdRhoDot_type, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter,
            class dRdGradRho_iter_out, class dRdV_iter_out, class dRdGradV_iter_out
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            result_type &result,
            dRdRho_type &dRdRho, dRdRhoDot_type &dRdRhoDot,
            dRdGradRho_iter_out dRdGradRho_begin, dRdGradRho_iter_out dRdGradRho_end,
            dRdV_iter_out dRdV_begin,             dRdV_iter_out dRdV_end,
            dRdGradV_iter_out dRdGradV_begin,     dRdGradV_iter_out dRdGradV_end
        );

        template<
            int dim, typename density_type, typename densityDot_type, typename result_type,
            typename testFunction_type, typename interpolationFunction_type,
            typename dRdRho_type, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class interpolationFunctionGradient_iter,
            class dRdU_iter_out, class dRdUMesh_iter_out,
            typename dDensityDotdDensity_type, typename dUDotdU_type
        >
        void computeBalanceOfMass(
            const density_type &density,                                  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,           const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                          const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,         const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &psi,                                 const interpolationFunction_type &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const dDensityDotdDensity_type &dDensityDotdDensity,          const dUDotdU_type &dUDotdU,
            result_type &result,
            dRdRho_type &dRdRho,
            dRdU_iter_out dRdU_begin,             dRdU_iter_out dRdU_end,
            dRdUMesh_iter_out dRdUMesh_begin,     dRdUMesh_iter_out dRdUMesh_end
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class result_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            result_iter_out result_start,                               result_iter_out result_stop
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            typename testFunction_type,
            class velocity_iter, class velocityGradient_iter, class result_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &psi,
            result_iter_out result_start,                               result_iter_out result_stop
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter,
            class result_iter_out,
            class dRdRho_iter_out, class dRdRhoDot_iter_out, class dRdGradRho_iter_out,
            class dRdV_iter_out, class dRdGradV_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                    const density_iter &density_end,
            const densityDot_iter &density_dot_begin,             const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            result_iter_out result_begin,         result_iter_out result_end,
            dRdRho_iter_out dRdRho_begin,         dRdRho_iter_out dRdRho_end,
            dRdRhoDot_iter_out dRdRhoDot_begin,   dRdRhoDot_iter_out dRdRhoDot_end,
            dRdGradRho_iter_out dRdGradRho_begin, dRdGradRho_iter_out dRdGradRho_end,
            dRdV_iter_out dRdV_begin,             dRdV_iter_out dRdV_end,
            dRdGradV_iter_out dRdGradV_begin,     dRdGradV_iter_out dRdGradV_end
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            typename testFunction_type, typename interpolationFunction_type,
            class velocity_iter, class velocityGradient_iter, class interpolationFunctionGradient_iter,
            class result_iter_out,
            class dRdRho_iter_out, class dRdU_iter_out, class dRdUMesh_iter_out,
            class dDensityDotdDensity_iter, class dUDotdU_iter
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                            const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                     const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,           const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                          const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,         const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &psi,                                 const interpolationFunction_type &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const dDensityDotdDensity_iter &dDensityDotdDensity_begin,    const dDensityDotdDensity_iter &dDensityDotdDensity_end,
            const dUDotdU_iter &dUDotdU_begin,                            const dUDotdU_iter &dUDotdU_end,
            result_iter_out result_begin,                                 result_iter_out result_end,
            dRdRho_iter_out dRdRho_begin,         dRdRho_iter_out dRdRho_end,
            dRdU_iter_out dRdU_begin,             dRdU_iter_out dRdU_end,
            dRdUMesh_iter_out dRdUMesh_begin,     dRdUMesh_iter_out dRdUMesh_end
        );

    }

}

#include "tardigrade_balance_of_mass.cpp"

#endif
