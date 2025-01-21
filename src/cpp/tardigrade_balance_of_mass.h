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
            const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &mass_change_rate
        );

        void computeBalanceOfMass(
            const floatType   &density,  const floatType         &density_dot,       const floatVector &density_gradient,
            const floatVector &velocity, const secondOrderTensor &velocity_gradient, floatType &mass_change_rate,
            floatType   &dCdRho, floatType         &dCdRhoDot, floatVector &dCdGradRho,
            floatVector &dCdV,   secondOrderTensor &dCdGradV
        );

        template<
            int dim, typename density_type, typename densityDot_type, typename c_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin, const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            c_type &mass_change_rate
        );

        template<
            int dim, typename density_type, typename densityDot_type, typename c_type,
            typename testFunction_type,
            class densityGradient_iter, class velocity_iter, class velocityGradient_iter
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin, const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &psi,
            c_type &mass_change_rate
        );

        template<
            int dim, typename density_type, typename densityDot_type, typename c_type,
            typename dCdRho_type, typename dCdRhoDot_type, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter,
            class dCdGradRho_iter_out, class dCdV_iter_out, class dCdGradV_iter_out
        >
        void computeBalanceOfMass(
            const density_type &density,  const densityDot_type &density_dot,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            c_type &mass_change_rate,
            dCdRho_type &dCdRho, dCdRhoDot_type &dCdRhoDot,
            dCdGradRho_iter_out dCdGradRho_begin, dCdGradRho_iter_out dCdGradRho_end,
            dCdV_iter_out dCdV_begin,             dCdV_iter_out dCdV_end,
            dCdGradV_iter_out dCdGradV_begin,     dCdGradV_iter_out dCdGradV_end
        );

        template<
            int dim, typename density_type, typename densityDot_type, typename c_type,
            typename testFunction_type, typename interpolationFunction_type,
            typename dCdRho_type, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class interpolationFunctionGradient_iter,
            class dCdU_iter_out, class dCdUMesh_iter_out,
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
            c_type &mass_change_rate,
            dCdRho_type &dCdRho,
            dCdU_iter_out dCdU_begin,             dCdU_iter_out dCdU_end,
            dCdUMesh_iter_out dCdUMesh_begin,     dCdUMesh_iter_out dCdUMesh_end
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter, class mass_change_rate_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            mass_change_rate_iter_out mass_change_rate_start,           mass_change_rate_iter_out mass_change_rate_stop
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            typename testFunction_type,
            class velocity_iter, class velocityGradient_iter, class mass_change_rate_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                          const density_iter &density_end,
            const densityDot_iter &density_dot_begin,                   const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,         const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                        const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin,       const velocityGradient_iter &velocity_gradient_end,
            const testFunction_type &psi,
            mass_change_rate_iter_out mass_change_rate_start,           mass_change_rate_iter_out mass_change_rate_stop
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            class velocity_iter, class velocityGradient_iter,
            class mass_change_rate_iter_out,
            class dCdRho_iter_out, class dCdRhoDot_iter_out, class dCdGradRho_iter_out,
            class dCdV_iter_out, class dCdGradV_iter_out
        >
        void computeBalanceOfMass(
            const density_iter &density_begin,                    const density_iter &density_end,
            const densityDot_iter &density_dot_begin,             const densityDot_iter &density_dot_end,
            const densityGradient_iter &density_gradient_begin,   const densityGradient_iter &density_gradient_end,
            const velocity_iter &velocity_begin,                  const velocity_iter &velocity_end,
            const velocityGradient_iter &velocity_gradient_begin, const velocityGradient_iter &velocity_gradient_end,
            mass_change_rate_iter_out mass_change_rate_begin,     mass_change_rate_iter_out mass_change_rate_end,
            dCdRho_iter_out dCdRho_begin,         dCdRho_iter_out dCdRho_end,
            dCdRhoDot_iter_out dCdRhoDot_begin,   dCdRhoDot_iter_out dCdRhoDot_end,
            dCdGradRho_iter_out dCdGradRho_begin, dCdGradRho_iter_out dCdGradRho_end,
            dCdV_iter_out dCdV_begin,             dCdV_iter_out dCdV_end,
            dCdGradV_iter_out dCdGradV_begin,     dCdGradV_iter_out dCdGradV_end
        );

        template<
            int dim, class density_iter, class densityDot_iter, class densityGradient_iter,
            typename testFunction_type, typename interpolationFunction_type,
            class velocity_iter, class velocityGradient_iter, class interpolationFunctionGradient_iter,
            class mass_change_rate_iter_out,
            class dCdRho_iter_out, class dCdU_iter_out, class dCdUMesh_iter_out,
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
            mass_change_rate_iter_out mass_change_rate_begin,     mass_change_rate_iter_out mass_change_rate_end,
            dCdRho_iter_out dCdRho_begin,         dCdRho_iter_out dCdRho_end,
            dCdU_iter_out dCdU_begin,             dCdU_iter_out dCdU_end,
            dCdUMesh_iter_out dCdUMesh_begin,     dCdUMesh_iter_out dCdUMesh_end
        );

    }

}

#include "tardigrade_balance_of_mass.cpp"

#endif
