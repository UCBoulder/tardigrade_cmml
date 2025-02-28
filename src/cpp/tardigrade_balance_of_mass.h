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

        inline void computeBalanceOfMass(
            const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
            const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &result
        );

        inline void computeBalanceOfMass(
            const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
            const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &result,
            floatType   &dRdRho, floatType         &dRdRhoDot, floatVector &dRdGradRho,
            floatVector &dRdV,   secondOrderTensor &dRdGradV
        );

        template<
            int dim, typename density_type, typename densityDot_type, typename result_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter
        >
        inline void computeBalanceOfMass(
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
        inline void computeBalanceOfMass(
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
        inline void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin, const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin,      const material_response_iter &material_response_end,
            const testFunction_type &psi,
            result_type &result
        );

        template<
            int dim, int matertial_response_dim, class density_iter, class densityDot_iter, class result_iter,
            typename testFunction_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter,
            class material_response_iter
        >
        inline void computeBalanceOfMass(
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
            class dRdGradRho_iter, class dRdV_iter, class dRdGradV_iter
        >
        inline void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            result_type &result,
            dRdRho_type &dRdRho, dRdRhoDot_type &dRdRhoDot,
            dRdGradRho_iter dRdGradRho_begin, dRdGradRho_iter dRdGradRho_end,
            dRdV_iter dRdV_begin,             dRdV_iter dRdV_end,
            dRdGradV_iter dRdGradV_begin,     dRdGradV_iter dRdGradV_end
        );

        template<
            int dim, typename density_type, typename densityDot_type, typename result_type,
            typename testFunction_type, typename interpolationFunction_type,
            typename dRdRho_type, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class interpolationFunctionGradient_iter,
            class dRdU_iter, class dRdUMesh_iter,
            typename dDensityDotdDensity_type, typename dUDotdU_type
        >
        inline void computeBalanceOfMass(
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
            dRdU_iter dRdU_begin,             dRdU_iter dRdU_end,
            dRdUMesh_iter dRdUMesh_begin,     dRdUMesh_iter dRdUMesh_end
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class result_iter
        >
        inline void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            result_iter result_start,                               result_iter result_stop
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            typename testFunction_type,
            class velocity_iter, class velocityGradient_iter, class result_iter
        >
        inline void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &psi,
            result_iter result_start,                               result_iter result_stop
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter,
            class result_iter,
            class dRdRho_iter, class dRdRhoDot_iter, class dRdGradRho_iter,
            class dRdV_iter, class dRdGradV_iter
        >
        inline void computeBalanceOfMass(
            const density_iter &density_begin,                    const density_iter &density_end,
            const densityDot_iter &density_dot_begin,             const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            result_iter result_begin,         result_iter result_end,
            dRdRho_iter dRdRho_begin,         dRdRho_iter dRdRho_end,
            dRdRhoDot_iter dRdRhoDot_begin,   dRdRhoDot_iter dRdRhoDot_end,
            dRdGradRho_iter dRdGradRho_begin, dRdGradRho_iter dRdGradRho_end,
            dRdV_iter dRdV_begin,             dRdV_iter dRdV_end,
            dRdGradV_iter dRdGradV_begin,     dRdGradV_iter dRdGradV_end
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            typename testFunction_type, typename interpolationFunction_type,
            class velocity_iter, class velocityGradient_iter, class interpolationFunctionGradient_iter,
            class result_iter,
            class dRdRho_iter, class dRdU_iter, class dRdUMesh_iter,
            class dDensityDotdDensity_iter, class dUDotdU_iter
        >
        inline void computeBalanceOfMass(
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
            result_iter result_begin,                                 result_iter result_end,
            dRdRho_iter dRdRho_begin,         dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin,             dRdU_iter dRdU_end,
            dRdUMesh_iter dRdUMesh_begin,     dRdUMesh_iter dRdUMesh_end
        );

        template<
            int dim, int mass_change_index, int material_response_dim, int material_response_num_dof,
            typename density_type, typename densityDot_type, typename result_type,
            typename testFunction_type, typename interpolationFunction_type,
            class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter,
            class material_response_iter, class material_response_jacobian_iter,
            class interpolationFunctionGradient_iter,
            class full_material_response_dof_gradient_iter,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter, class dRdTheta_iter, class dRdE_iter, class dRdZ_iter,
            class dRdUMesh_iter,
            typename dDensityDotdDensity_type, typename dUDotdU_type,
            int density_index         = 0,
            int displacement_index    = 1,
            int velocity_index        = 4,
            int temperature_index     = 7,
            int internal_energy_index = 8,
            int additional_dof_index  = 9
        >
        inline void computeBalanceOfMass(
            const density_type &density,                                  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,           const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                          const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,         const velocityGradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin,        const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const testFunction_type &test_function,                       const interpolationFunction_type &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dDensityDotdDensity_type &dDensityDotdDensity,          const dUDotdU_type &dUDotdU,
            const unsigned int phase,
            result_type &result,
            dRdRho_iter dRdRho_begin,         dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin,             dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin,             dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin,     dRdTheta_iter dRdTheta_end,
            dRdE_iter dRdE_begin,             dRdE_iter dRdE_end,
            dRdZ_iter dRdZ_begin,             dRdZ_iter dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin,     dRdUMesh_iter dRdUMesh_end
        );

        template<
            int dim, int mass_change_index, int material_response_dim, int material_response_num_dof,
            class density_iter, class densityDot_iter, class result_iter,
            typename testFunction_type, typename interpolationFunction_type,
            class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter,
            class material_response_iter, class material_response_jacobian_iter,
            class interpolationFunctionGradient_iter,
            class full_material_response_dof_gradient_iter,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter, class dRdTheta_iter, class dRdE_iter, class dRdZ_iter,
            class dRdUMesh_iter,
            typename dDensityDotdDensity_type, typename dUDotdU_type,
            int density_index         = 0,
            int displacement_index    = 1,
            int velocity_index        = 4,
            int temperature_index     = 7,
            int internal_energy_index = 8,
            int additional_dof_index  = 9
        >
        inline void computeBalanceOfMass(
            const density_iter &density_begin,                            const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                     const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,           const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                          const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,         const velocityGradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin,        const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const testFunction_type &test_function,                       const interpolationFunction_type &interpolation_function,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_begin,
            const interpolationFunctionGradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dDensityDotdDensity_type &dDensityDotdDensity,          const dUDotdU_type &dUDotdU,
            result_iter result_begin,         result_iter result_end,
            dRdRho_iter dRdRho_begin,         dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin,             dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin,             dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin,     dRdTheta_iter dRdTheta_end,
            dRdE_iter dRdE_begin,             dRdE_iter dRdE_end,
            dRdZ_iter dRdZ_begin,             dRdZ_iter dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin,     dRdUMesh_iter dRdUMesh_end
        );

    }

}

#include "tardigrade_balance_of_mass.cpp"

#endif
