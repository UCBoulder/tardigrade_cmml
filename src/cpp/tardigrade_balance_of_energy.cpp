/**
  ******************************************************************************
  * \file tardigrade_balance_of_energy.cpp
  ******************************************************************************
  * The source file for the equations associated with the balance of energy
  ******************************************************************************
  */

#include "tardigrade_balance_of_energy.h"
#include "tardigrade_balance_of_mass.h"
#include<numeric>
#include<algorithm>
#include<functional>

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
        ){
            /*!
             * Compute the full balance of energy in a variational context
             * 
             * \param &density: The apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot: The partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy: The internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot: The partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction: The starting iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation: The internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &result: The result of the non-divergence part of the balance of energy
             */

            computeBalanceOfEnergyNonDivergence<dim>(
                density, density_dot, density_gradient_begin, density_gradient_end,
                internal_energy, internal_energy_dot, internal_energy_gradient_begin, internal_energy_gradient_end,
                velocity_begin, velocity_end, velocity_gradient_begin, velocity_gradient_end,
                cauchy_stress_begin, cauchy_stress_end,
                volume_fraction,
                internal_heat_generation,
                net_interphase_force_begin, net_interphase_force_end,
                result
            );

            result *= test_function;

            result_type result_divergence;

            computeBalanceOfEnergyDivergence<dim>(
                test_function_gradient_begin, test_function_gradient_end,
                heat_flux_begin, heat_flux_end,
                result_divergence
            );

            result += result_divergence;

        }

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
        ){
            /*!
             * Compute the balance of energy.
             * 
             * \param &density: The apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot: The partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy: The internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot: The partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction: The starting iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation: The internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &interpolation_function: The value of the interpolation function \f$ \left( \phi \right) \f$
             * \param &interpolation_function_gradient_begin: The starting iterator of the spatial gradient of the interpolation function \f$ \left( \phi_{,i} \right) \f$
             * \param &interpolation_function_gradient_end: The stopping iterator of the spatial gradient of the interpolation function \f$ \left( \phi_{,i} \right) \f$
             * \param &dRhoDotdRho: The derivative of the time rate of change of the density w.r.t. the density
             * \param &dEDotdE: The derivative of the time rate of change of the internal energy w.r.t. the internal energy
             * \param &dUDotdU: The derivative of the time rate of change of the displacement degree of freedom w.r.t. the displacement degree of freedom
             * \param &result: The result of the non-divergence part of the balance of energy
             * \param &dRdRho: The derivative of the residual w.r.t. the density
             * \param &dRdE: The derivative of the residual w.r.t. the internal energy
             * \param &dRdU_begin: The starting iterator derivative of the residual w.r.t. the displacement degree of freedom
             * \param &dRdU_end: The stopping iterator derivative of the residual w.r.t. the displacement degree of freedom
             * \param &dRdCauchy_begin: The starting iterator of the derivative of the residual w.r.t. the Cauchy stress
             * \param &dRdCauchy_end: The stopping iterator of the derivative of the residual w.r.t. the Cauchy stress
             * \param &dRdVolumeFraction: The derivative of the residual w.r.t. the volume fraction
             * \param &dRdr: The derivative of the residual w.r.t. the internal heat generation
             * \param &dRdpi_begin: The starting iterator of the derivative of the residual w.r.t. the net interphase force
             * \param &dRdpi_end: The stopping iterator of the derivative of the residual w.r.t. the net interphase force
             * \param &dRdq_begin: The starting iterator of the derivative of the residual w.r.t. the heat flux
             * \param &dRdq_end: The stopping iterator of the derivative of the residual w.r.t. the heat flux
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh displacement
             */

            using dRdRhoDot_type  = result_type;
            using dRdGradRho_type = result_type;

            using dRdEDot_type    = result_type;
            using dRdGradE_type   = result_type;

            using dRdU_type       = typename std::iterator_traits<dRdU_iter>::value_type;
            using dRdGradV_type   = result_type;

            using dRdpi_type       = typename std::iterator_traits<dRdpi_iter>::value_type;

            using dRdGradTestFunction_type = result_type;

            dRdRhoDot_type dRdRhoDot;
            std::array< dRdGradRho_type, dim > dRdGradRho;

            dRdEDot_type dRdEDot;
            std::array< dRdGradE_type, dim > dRdGradE;

            std::array< dRdGradV_type, dim * dim > dRdGradV;

            std::array< dRdGradTestFunction_type, dim > dRdGradTestFunction;

            computeBalanceOfEnergyNonDivergence<dim>(
                density, density_dot, density_gradient_begin, density_gradient_end,
                internal_energy, internal_energy_dot, internal_energy_gradient_begin, internal_energy_gradient_end,
                velocity_begin, velocity_end, velocity_gradient_begin, velocity_gradient_end,
                cauchy_stress_begin, cauchy_stress_end,
                volume_fraction,
                internal_heat_generation,
                net_interphase_force_begin, net_interphase_force_end,
                result,
                dRdRho, dRdRhoDot, std::begin( dRdGradRho ), std::end( dRdGradRho ),
                dRdE, dRdEDot, std::begin( dRdGradE ), std::end( dRdGradE ),
                dRdU_begin, dRdU_end, std::begin( dRdGradV ), std::end( dRdGradV ),
                dRdCauchy_begin, dRdCauchy_end, dRdVolumeFraction, dRdr,
                dRdpi_begin, dRdpi_end
            );

            result *= test_function;

            result_type result_divergence;

            computeBalanceOfEnergyDivergence<dim>(
                test_function_gradient_begin, test_function_gradient_end,
                heat_flux_begin, heat_flux_end,
                result_divergence,
                std::begin( dRdGradTestFunction ), std::end( dRdGradTestFunction ),
                dRdq_begin,                        dRdq_end
            );

            result += result_divergence;

            // Assemble the Jacobians

            // Jacobian w.r.t. density
            dRdRho *= test_function * interpolation_function;
            dRdRho += test_function * dRdRhoDot * dRhoDotdRho * interpolation_function;
            dRdRho += test_function * std::inner_product( std::begin( dRdGradRho ), std::end( dRdGradRho ), interpolation_function_gradient_begin, result_type( ) );

            // Jacobian w.r.t. internal energy
            dRdE *= test_function * interpolation_function;
            dRdE += test_function * dRdEDot * dEDotdE * interpolation_function;
            dRdE += test_function * std::inner_product( std::begin( dRdGradE ), std::end( dRdGradE ), interpolation_function_gradient_begin, result_type( ) );

            // Jacobian w.r.t. deformation DOF
            std::transform(
                dRdU_begin, dRdU_end, dRdU_begin,
                std::bind(
                    std::multiplies< dRdU_type >( ),
                    std::placeholders::_1,
                    test_function * dUDotdU * interpolation_function
                )
            );

            for ( unsigned int i = 0; i < dim; ++i ){
                for ( unsigned int j = 0; j < dim; ++j ){
                    *( dRdU_begin + i ) += test_function * dRdGradV[ dim * i + j ] * dUDotdU * ( *( interpolation_function_gradient_begin + j ) );
                }
            }

            // Jacobian w.r.t. Cauchy stress
            std::transform(
                dRdCauchy_begin, dRdCauchy_end, dRdCauchy_begin,
                std::bind(
                    std::multiplies< dRdU_type >( ),
                    std::placeholders::_1,
                    test_function
                )
            );

            // Jacobian w.r.t. the volume fraction
            dRdVolumeFraction *= test_function;

            // Jacobian w.r.t. the internal heat generation
            dRdr *= test_function;

            // Jacobian w.r.t. the net interphase force
            std::transform(
                dRdpi_begin, dRdpi_end, dRdpi_begin,
                std::bind(
                    std::multiplies< dRdpi_type >( ),
                    std::placeholders::_1,
                    test_function
                )
            );

            // Jacobian w.r.t. the mesh deformation
            std::fill( dRdUMesh_begin, dRdUMesh_end, result_type( ) );
            for ( unsigned int a = 0; a < dim; ++a ){

                *( dRdUMesh_begin + a ) += result * ( *( interpolation_function_gradient_begin + a ) );

                for ( unsigned int i = 0; i < dim; ++i ){

                    *( dRdUMesh_begin + a ) -= test_function *
                    (
                        dRdGradRho[ i ] * ( *( density_gradient_begin         + a ) ) +
                        dRdGradE[ i ]   * ( *( internal_energy_gradient_begin + a ) )
                    ) * ( *( interpolation_function_gradient_begin + i ) );

                    *( dRdUMesh_begin + a ) -= dRdGradTestFunction[ i ] * ( *( test_function_gradient_begin + a ) ) * ( *( interpolation_function_gradient_begin + i ) );

                    for ( unsigned int j = 0; j < dim; ++j ){

                        *( dRdUMesh_begin + a ) -= test_function *
                        (
                            dRdGradV[ dim * i + j ] * ( *( velocity_gradient_begin + dim * i + a ) ) * ( *( interpolation_function_gradient_begin + j ) )
                        );

                    }

                }

            }

        }

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
        ){
            /*!
             * Compute the multiphase balance of energy
             * 
             * \param &density_begin: The starting iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_end: The stopping iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot_begin: The starting iterator of the partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_dot_end: The stopping iterator of the partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_begin: The starting iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_end: The stopping iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot_begin: The starting iterator of the partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_dot_end: The stopping iterator of the partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation_begin: The starting iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &internal_heat_generation_end: The stopping iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &result_begin: The starting iterator of the result of the non-divergence part of the balance of energy
             * \param &result_end: The stopping iterator of the result of the non-divergence part of the balance of energy
             */

            TARDIGRADE_ERROR_TOOLS_EVAL( unsigned int nphases = ( unsigned int )( density_end - density_begin ); )

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density_dot must be the same size" )

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases * dim == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient terms are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( internal_energy_end - internal_energy_begin ), "The internal energy and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( internal_energy_dot_end - internal_energy_dot_begin ), "The internal energy dot and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases * dim == ( unsigned int )( internal_energy_gradient_end - internal_energy_gradient_begin ), "The density and internal energy gradient terms are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases * dim == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity terms are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases * dim * dim == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient terms are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases * dim * dim == ( unsigned int )( cauchy_stress_end - cauchy_stress_begin ), "The density and Cauchy stress terms are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( volume_fraction_end - volume_fraction_begin ), "The volume fraction and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( internal_heat_generation_end - internal_heat_generation_begin ), "The internal heat generation and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * nphases == ( unsigned int )( net_interphase_force_end - net_interphase_force_begin ), "The net interphase force and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * nphases == ( unsigned int )( heat_flux_end - heat_flux_begin ), "The heat flux and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( result_end - result_begin ), "The result and density must be the same size" );

            for ( auto rho = std::pair< unsigned int, density_iter >( 0, density_begin ); rho.second != density_end; ++rho.first, ++rho.second ){

                computeBalanceOfEnergy<dim>(
                    *rho.second,                                      *( density_dot_begin + rho.first ),
                    density_gradient_begin + dim * rho.first,         density_gradient_begin + dim * ( rho.first + 1 ),
                    *( internal_energy_begin + rho.first ),           *( internal_energy_dot_begin + rho.first ),
                    internal_energy_gradient_begin + dim * rho.first, internal_energy_gradient_begin + dim * ( rho.first + 1 ),
                    velocity_begin + dim * rho.first,                 velocity_begin + dim * ( rho.first + 1 ),
                    velocity_gradient_begin + dim * dim * rho.first,  velocity_gradient_begin + dim * dim * ( rho.first + 1 ),
                    cauchy_stress_begin + dim * dim * rho.first,      cauchy_stress_begin + dim * dim * ( rho.first + 1 ),
                    *( volume_fraction_begin + rho.first ),
                    *( internal_heat_generation_begin + rho.first ),
                    net_interphase_force_begin + dim * rho.first,     net_interphase_force_begin + dim * ( rho.first + 1 ),
                    heat_flux_begin + dim * rho.first,                heat_flux_begin + dim * ( rho.first + 1 ),
                    test_function,
                    test_function_gradient_begin,                     test_function_gradient_end,
                    *( result_begin + rho.first )
                );

            }

        }

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
        ){
            /*!
             * Compute the multiphase balance of energy.
             * 
             * \param &density_begin: The starting iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_end: The stopping iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot_begin: The starting iterator of the partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_dot_end: The stopping iterator of the partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_begin: The starting iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_end: The stopping iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot_begin: The starting iterator of the partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_dot_end: The stopping iterator of the partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation_begin: The starting iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &internal_heat_generation_end: The stopping iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &interpolation_function: The value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the spatial gradient of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the spatial gradient of the interpolation function
             * \param &dRhoDotdRho: The derivative of the time rate of change of the density w.r.t. the density
             * \param &dEDotdE: The derivative of the time rate of change of the internal energy w.r.t. the internal energy
             * \param &dUDotdU: The derivative of the time rate of change of the displacement degree of freedom w.r.t. the displacement degree of freedom
             * \param &result_begin: The starting iterator of the result of the non-divergence part of the balance of energy
             * \param &result_end: The stopping iterator of the result of the non-divergence part of the balance of energy
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the density
             * \param &dRdE_begin: The starting iterator of the derivative of the residual w.r.t. the internal energy
             * \param &dRdE_end: The stopping iterator of the derivative of the residual w.r.t. the internal energy
             * \param &dRdU_begin: The starting iterator derivative of the residual w.r.t. the displacement degree of freedom
             * \param &dRdU_end: The stopping iterator derivative of the residual w.r.t. the displacement degree of freedom
             * \param &dRdCauchy_begin: The starting iterator of the derivative of the residual w.r.t. the Cauchy stress
             * \param &dRdCauchy_end: The stopping iterator of the derivative of the residual w.r.t. the Cauchy stress
             * \param &dRdVolumeFraction_begin: The starting iterator of the derivative of the residual w.r.t. the volume fraction
             * \param &dRdVolumeFraction_end: The stopping iterator of the derivative of the residual w.r.t. the volume fraction
             * \param &dRdr_begin: The starting iterator of the derivative of the residual w.r.t. the internal heat generation
             * \param &dRdr_end: The stopping iterator of the derivative of the residual w.r.t. the internal heat generation
             * \param &dRdpi_begin: The starting iterator of the derivative of the residual w.r.t. the net interphase force
             * \param &dRdpi_end: The stopping iterator of the derivative of the residual w.r.t. the net interphase force
             * \param &dRdq_begin: The starting iterator of the derivative of the residual w.r.t. the heat flux
             * \param &dRdq_end: The stopping iterator of the derivative of the residual w.r.t. the heat flux
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh displacement
             */

            TARDIGRADE_ERROR_TOOLS_EVAL( unsigned int nphases = ( unsigned int )( density_end - density_begin ); )

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density_dot must be the same size" )

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases * dim == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient terms are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( internal_energy_end - internal_energy_begin ), "The internal energy and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( internal_energy_dot_end - internal_energy_dot_begin ), "The internal energy dot and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases * dim == ( unsigned int )( internal_energy_gradient_end - internal_energy_gradient_begin ), "The density and internal energy gradient terms are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases * dim == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity terms are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases * dim * dim == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient terms are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases * dim * dim == ( unsigned int )( cauchy_stress_end - cauchy_stress_begin ), "The density and Cauchy stress terms are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( volume_fraction_end - volume_fraction_begin ), "The volume fraction and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( internal_heat_generation_end - internal_heat_generation_begin ), "The internal heat generation and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * nphases == ( unsigned int )( net_interphase_force_end - net_interphase_force_begin ), "The net interphase force and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * nphases == ( unsigned int )( heat_flux_end - heat_flux_begin ), "The heat flux and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( result_end - result_begin ), "The result and density must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( dRdRho_end - dRdRho_begin ), "The result and dRdRho must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( dRdRho_end - dRdRho_begin ), "The result and dRdE must be the same size" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * nphases == ( unsigned int )( dRdU_end - dRdU_begin ), "The result and dRdU are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( dim * dim * nphases == ( unsigned int )( dRdCauchy_end - dRdCauchy_begin ), "The result and dRdCauchy are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( dRdVolumeFraction_end - dRdVolumeFraction_begin ), "The result and dRdVolumeFraction are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( nphases == ( unsigned int )( dRdr_end - dRdr_begin ), "The result and dRdr are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( 3 * nphases == ( unsigned int )( dRdpi_end - dRdpi_begin ), "The result and dRdpi are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( 3 * nphases == ( unsigned int )( dRdq_end - dRdq_begin ), "The result and dRdq are of inconsistent sizes" );

            TARDIGRADE_ERROR_TOOLS_CHECK( 3 * nphases == ( unsigned int )( dRdUMesh_end - dRdUMesh_begin ), "The result and dRdUMesh are of inconsistent sizes" );

            for ( auto rho = std::pair< unsigned int, density_iter >( 0, density_begin ); rho.second != density_end; ++rho.first, ++rho.second ){

                computeBalanceOfEnergy<dim>(
                    *rho.second,                                      *( density_dot_begin + rho.first ),
                    density_gradient_begin + dim * rho.first,         density_gradient_begin + dim * ( rho.first + 1 ),
                    *( internal_energy_begin + rho.first ),           *( internal_energy_dot_begin + rho.first ),
                    internal_energy_gradient_begin + dim * rho.first, internal_energy_gradient_begin + dim * ( rho.first + 1 ),
                    velocity_begin + dim * rho.first,                 velocity_begin + dim * ( rho.first + 1 ),
                    velocity_gradient_begin + dim * dim * rho.first,  velocity_gradient_begin + dim * dim * ( rho.first + 1 ),
                    cauchy_stress_begin + dim * dim * rho.first,      cauchy_stress_begin + dim * dim * ( rho.first + 1 ),
                    *( volume_fraction_begin + rho.first ),
                    *( internal_heat_generation_begin + rho.first ),
                    net_interphase_force_begin + dim * rho.first,     net_interphase_force_begin + dim * ( rho.first + 1 ),
                    heat_flux_begin + dim * rho.first,                heat_flux_begin + dim * ( rho.first + 1 ),
                    test_function,
                    test_function_gradient_begin,                     test_function_gradient_end,
                    interpolation_function,
                    interpolation_function_gradient_begin,            interpolation_function_gradient_end,
                    dRhoDotdRho, dEDotdE, dUDotdU,
                    *( result_begin + rho.first ),
                    *( dRdRho_begin + rho.first ),
                    *( dRdE_begin + rho.first ),
                    dRdU_begin + dim * rho.first,                     dRdU_begin + dim * ( rho.first + 1 ),
                    dRdCauchy_begin + dim * dim * rho.first,          dRdCauchy_begin + dim * dim * ( rho.first + 1 ),
                    *( dRdVolumeFraction_begin + rho.first ),
                    *( dRdr_begin + rho.first ),
                    dRdpi_begin + dim * rho.first,                    dRdpi_begin + dim * ( rho.first + 1 ),
                    dRdq_begin + dim * rho.first,                     dRdq_begin + dim * ( rho.first + 1 ),
                    dRdUMesh_begin + dim * rho.first,                 dRdUMesh_begin + dim * ( rho.first + 1 )
                );

            }

        }

        template<
            int dim, int material_response_dim, int cauchy_stress_index, int internal_heat_generation_index, int heat_flux_index,
            int interphasic_force_index, int interphasic_heat_transfer_index,
            typename density_type, typename density_dot_type,
            class density_gradient_iter,
            typename internal_energy_type, typename internal_energy_dot_type,
            class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class material_response_iter,
            typename volume_fraction_type,
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
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const volume_fraction_type &volume_fraction,
            const test_function_type &test_function,
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            result_type &result
        ){
            /*!
             * Compute the full balance of energy in a variational context using a generalized material response vector
             * 
             * \param &density: The apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot: The partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy: The internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot: The partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &volume_fraction: The the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &result: The result of the non-divergence part of the balance of energy
             */

            computeBalanceOfEnergy<dim>(
                density,         density_dot,         density_gradient_begin,         density_gradient_end,
                internal_energy, internal_energy_dot, internal_energy_gradient_begin, internal_energy_gradient_end,
                velocity_begin,  velocity_end,        velocity_gradient_begin,        velocity_gradient_end,
                material_response_begin + cauchy_stress_index,
                material_response_begin + cauchy_stress_index + material_response_dim * material_response_dim,
                volume_fraction,
                *( material_response_begin + internal_heat_generation_index ),
                material_response_begin + interphasic_force_index, material_response_begin + interphasic_force_index + material_response_dim,
                material_response_begin + heat_flux_index,         material_response_begin + heat_flux_index + material_response_dim,
                test_function,
                test_function_gradient_begin, test_function_gradient_end,
                result
            );

//            result -= test_function * ( *( material_response_begin + interphasic_heat_transfer_index ) );

        }

        template<
            int dim, int material_response_dim, int cauchy_stress_index, int internal_heat_generation_index, int heat_flux_index,
            int interphasic_force_index, int interphasic_heat_transfer_index, int material_response_num_dof,
            typename density_type, typename density_dot_type,
            class density_gradient_iter,
            typename internal_energy_type, typename internal_energy_dot_type,
            class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class material_response_iter,
            class material_response_jacobian_iter,
            typename volume_fraction_type,
            typename test_function_type, class test_function_gradient_iter,
            typename interpolation_function_type, class interpolation_function_gradient_iter,
            class full_material_response_dof_gradient_iter,
            typename dRhoDotdRho_type, typename dEDotdE_type, typename dUDotdU_type,
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
        void computeBalanceOfEnergy(
            const density_type &density, const density_dot_type &density_dot,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const internal_energy_type &internal_energy, const internal_energy_dot_type &internal_energy_dot,
            const internal_energy_gradient_iter &internal_energy_gradient_begin, const internal_energy_gradient_iter &internal_energy_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const volume_fraction_type &volume_fraction,
            const test_function_type &test_function,
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin, const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dRhoDotdRho_type &dRhoDotdRho, const dEDotdE_type dEDotdE, const dUDotdU_type &dUDotdU,
            const unsigned int &phase,
            result_type &result,
            dRdRho_iter dRdRho_begin,                       dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin,                           dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin,                           dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin,                   dRdTheta_iter dRdTheta_end,
            dRdE_iter dRdE_begin,                           dRdE_iter dRdE_end,
            dRdZ_iter dRdZ_begin,                           dRdZ_iter dRdZ_end,
            dRdVolumeFraction_iter dRdVolumeFraction_begin, dRdVolumeFraction_iter dRdVolumeFraction_end,
            dRdUMesh_iter dRdUMesh_begin,                   dRdUMesh_iter dRdUMesh_end
        ){
            /*!
             * Compute the full balance of energy in a variational context using a generalized material response vector
             * 
             * \param &density: The apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot: The partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy: The internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot: The partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &volume_fraction: The the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &interpolation_function: The value of the interpolation function \f$ \left( \phi \right) \f$
             * \param &interpolation_function_gradient_begin: The starting iterator of the spatial gradient of the interpolation function \f$ \left( \phi_{,i} \right) \f$
             * \param &interpolation_function_gradient_end: The stopping iterator of the spatial gradient of the interpolation function \f$ \left( \phi_{,i} \right) \f$
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the material response dof vector
             * \param &dRhoDotdRho: The derivative of the time rate of change of the density w.r.t. the density
             * \param &dEDotdE: The derivative of the time rate of change of the internal energy w.r.t. the internal energy
             * \param &dUDotdU: The derivative of the time rate of change of the displacement degree of freedom w.r.t. the displacement degree of freedom
             * \param &phase: The phase that is having the balance of energy computed on
             * \param &result: The result of the non-divergence part of the balance of energy
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the density
             * \param &dRdU_begin: The starting iterator derivative of the residual w.r.t. the spatial degree of freedom
             * \param &dRdU_end: The stopping iterator derivative of the residual w.r.t. the spatial degree of freedom
             * \param &dRdW_begin: The starting iterator derivative of the residual w.r.t. the displacement degree of freedom
             * \param &dRdW_end: The stopping iterator derivative of the residual w.r.t. the displacement degree of freedom
             * \param &dRdTheta_begin: The starting iterator of the derivative of the residual w.r.t. the temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the residual w.r.t. the temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the residual w.r.t. the internal energy
             * \param &dRdE_end: The stopping iterator of the derivative of the residual w.r.t. the internal energy
             * \param &dRdZ_begin: The starting iterator of the derivative of the residual w.r.t. the additional degrees of freedom
             * \param &dRdZ_end: The stopping iterator of the derivative of the residual w.r.t. the additional degrees of freedom
             * \param &dRdVolumeFraction_begin: The starting iterator of the derivative of the residual w.r.t. the volume fraction
             * \param &dRdVolumeFraction_end: The stopping iterator of the derivative of the residual w.r.t. the volume fraction
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh displacement
             */

            std::array< result_type, material_response_dim * material_response_dim > dRdCauchy_phase;

            result_type dRdr_phase;

            std::array< result_type, material_response_dim > dRdpi_phase, dRdq_phase;

            std::fill( dRdRho_begin,            dRdRho_end,            0 );
            std::fill( dRdU_begin,              dRdU_end,              0 );
            std::fill( dRdW_begin,              dRdW_end,              0 );
            std::fill( dRdTheta_begin,          dRdTheta_end,          0 );
            std::fill( dRdE_begin,              dRdE_end,              0 );
            std::fill( dRdZ_begin,              dRdZ_end,              0 );
            std::fill( dRdVolumeFraction_begin, dRdVolumeFraction_end, 0 );
            std::fill( dRdUMesh_begin,          dRdUMesh_end,          0 );

            computeBalanceOfEnergy<dim>(
                density,         density_dot,         density_gradient_begin,         density_gradient_end,
                internal_energy, internal_energy_dot, internal_energy_gradient_begin, internal_energy_gradient_end,
                velocity_begin,  velocity_end,        velocity_gradient_begin,        velocity_gradient_end,
                material_response_begin + cauchy_stress_index,
                material_response_begin + cauchy_stress_index + material_response_dim * material_response_dim,
                volume_fraction,
                *( material_response_begin + internal_heat_generation_index ),
                material_response_begin + interphasic_force_index,
                material_response_begin + interphasic_force_index + material_response_dim,
                material_response_begin + heat_flux_index,
                material_response_begin + heat_flux_index + material_response_dim,
                test_function,          test_function_gradient_begin,          test_function_gradient_end,
                interpolation_function, interpolation_function_gradient_begin, interpolation_function_gradient_end,
                dRhoDotdRho, dEDotdE, dUDotdU,
                result,
                *( dRdRho_begin + phase ),            *( dRdE_begin + phase ),
                dRdU_begin   + dim * phase,           dRdU_begin + dim * ( phase + 1 ),
                std::begin( dRdCauchy_phase ),        std::end( dRdCauchy_phase ),
                *( dRdVolumeFraction_begin + phase ), dRdr_phase,
                std::begin( dRdpi_phase ),            std::end( dRdpi_phase ),
                std::begin( dRdq_phase ),             std::end( dRdq_phase ),
                dRdUMesh_begin,                       dRdUMesh_end
            );

//            result -= test_function * ( *( material_response_begin + interphasic_heat_transfer_index ) );

            const unsigned int nphases = ( unsigned int )( dRdRho_end - dRdRho_begin );

            // Scale the volume fraction by the interpolation function
            *( dRdVolumeFraction_begin + phase ) *= interpolation_function;

            // CAUCHY STRESS CONTRIBUTIONS
            for ( unsigned int j = 0; j < dim * dim; ++j ){

                // density
                for ( auto p = std::pair< unsigned int, dRdRho_iter >( 0, dRdRho_begin ); p.second != dRdRho_end; ++p.first, ++p.second ){

                    *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + nphases * density_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + material_response_num_dof + material_response_dim * ( nphases * density_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // velocity
                for ( auto p = std::pair< unsigned int, dRdU_iter >( 0, dRdU_begin ); p.second != dRdU_end; ++p.first, ++p.second ){

                    *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + nphases * velocity_index + p.first ) ) * interpolation_function * dUDotdU;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + material_response_num_dof + material_response_dim * ( nphases * velocity_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) ) * dUDotdU;

                    }

                }

                // displacement
                for ( auto p = std::pair< unsigned int, dRdW_iter >( 0, dRdW_begin ); p.second != dRdW_end; ++p.first, ++p.second ){

                    *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + nphases * displacement_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + material_response_num_dof + material_response_dim * ( nphases * displacement_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // temperature
                for ( auto p = std::pair< unsigned int, dRdTheta_iter >( 0, dRdTheta_begin ); p.second != dRdTheta_end; ++p.first, ++p.second ){

                    *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + nphases * temperature_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + material_response_num_dof + material_response_dim * ( nphases * temperature_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // internal energy
                for ( auto p = std::pair< unsigned int, dRdE_iter >( 0, dRdE_begin ); p.second != dRdE_end; ++p.first, ++p.second ){

                    *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + nphases * internal_energy_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + material_response_num_dof + material_response_dim * ( nphases * internal_energy_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // additional dof
                for ( auto p = std::pair< unsigned int, dRdZ_iter >( 0, dRdZ_begin ); p.second != dRdZ_end; ++p.first, ++p.second ){

                    *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + nphases * additional_dof_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdCauchy_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + material_response_num_dof + material_response_dim * ( nphases * additional_dof_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // mesh displacement
                for ( unsigned int I = 0; I < material_response_num_dof; ++I ){

                    for ( unsigned int k = 0; k < material_response_dim; ++k ){

                        for ( unsigned int a = 0; a < material_response_dim; ++a ){

                            *( dRdUMesh_begin + a ) -=
                                dRdCauchy_phase[ j ]
                                * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( cauchy_stress_index + j ) + material_response_num_dof + material_response_dim * I + k ) )
                                * ( *( full_material_response_dof_gradient_begin + material_response_dim * I + a ) )
                                * ( *( interpolation_function_gradient_begin + k ) );

                        }

                    }

                }

            }

            // INTERNAL HEAT GENERATION
            // density
            for ( auto p = std::pair< unsigned int, dRdRho_iter >( 0, dRdRho_begin ); p.second != dRdRho_end; ++p.first, ++p.second ){

                *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + nphases * density_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + material_response_num_dof + material_response_dim * ( nphases * density_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // velocity
            for ( auto p = std::pair< unsigned int, dRdU_iter >( 0, dRdU_begin ); p.second != dRdU_end; ++p.first, ++p.second ){

                *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + nphases * velocity_index + p.first ) ) * interpolation_function * dUDotdU;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + material_response_num_dof + material_response_dim * ( nphases * velocity_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) ) * dUDotdU;

                }

            }

            // displacement
            for ( auto p = std::pair< unsigned int, dRdW_iter >( 0, dRdW_begin ); p.second != dRdW_end; ++p.first, ++p.second ){

                *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + nphases * displacement_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + material_response_num_dof + material_response_dim * ( nphases * displacement_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // temperature
            for ( auto p = std::pair< unsigned int, dRdTheta_iter >( 0, dRdTheta_begin ); p.second != dRdTheta_end; ++p.first, ++p.second ){

                *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + nphases * temperature_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + material_response_num_dof + material_response_dim * ( nphases * temperature_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // internal energy
            for ( auto p = std::pair< unsigned int, dRdE_iter >( 0, dRdE_begin ); p.second != dRdE_end; ++p.first, ++p.second ){

                *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + nphases * internal_energy_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + material_response_num_dof + material_response_dim * ( nphases * internal_energy_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // additional dof
            for ( auto p = std::pair< unsigned int, dRdZ_iter >( 0, dRdZ_begin ); p.second != dRdZ_end; ++p.first, ++p.second ){

                *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + nphases * additional_dof_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += dRdr_phase * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + material_response_num_dof + material_response_dim * ( nphases * additional_dof_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // mesh displacement
            for ( unsigned int I = 0; I < material_response_num_dof; ++I ){

                for ( unsigned int k = 0; k < material_response_dim; ++k ){

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *( dRdUMesh_begin + a ) -=
                            dRdr_phase
                            * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( internal_heat_generation_index ) + material_response_num_dof + material_response_dim * I + k ) )
                            * ( *( full_material_response_dof_gradient_begin + material_response_dim * I + a ) )
                            * ( *( interpolation_function_gradient_begin + k ) );

                    }

                }

            }

            // INTERPHASIC FORCE CONTRIBUTIONS
            for ( unsigned int j = 0; j < dim; ++j ){

                // density
                for ( auto p = std::pair< unsigned int, dRdRho_iter >( 0, dRdRho_begin ); p.second != dRdRho_end; ++p.first, ++p.second ){

                    *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + nphases * density_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + material_response_num_dof + material_response_dim * ( nphases * density_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // velocity
                for ( auto p = std::pair< unsigned int, dRdU_iter >( 0, dRdU_begin ); p.second != dRdU_end; ++p.first, ++p.second ){

                    *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + nphases * velocity_index + p.first ) ) * interpolation_function * dUDotdU;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + material_response_num_dof + material_response_dim * ( nphases * velocity_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) ) * dUDotdU;

                    }

                }

                // displacement
                for ( auto p = std::pair< unsigned int, dRdW_iter >( 0, dRdW_begin ); p.second != dRdW_end; ++p.first, ++p.second ){

                    *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + nphases * displacement_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + material_response_num_dof + material_response_dim * ( nphases * displacement_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // temperature
                for ( auto p = std::pair< unsigned int, dRdTheta_iter >( 0, dRdTheta_begin ); p.second != dRdTheta_end; ++p.first, ++p.second ){

                    *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + nphases * temperature_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + material_response_num_dof + material_response_dim * ( nphases * temperature_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // internal energy
                for ( auto p = std::pair< unsigned int, dRdE_iter >( 0, dRdE_begin ); p.second != dRdE_end; ++p.first, ++p.second ){

                    *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + nphases * internal_energy_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + material_response_num_dof + material_response_dim * ( nphases * internal_energy_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // additional dof
                for ( auto p = std::pair< unsigned int, dRdZ_iter >( 0, dRdZ_begin ); p.second != dRdZ_end; ++p.first, ++p.second ){

                    *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + nphases * additional_dof_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdpi_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + material_response_num_dof + material_response_dim * ( nphases * additional_dof_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // mesh displacement
                for ( unsigned int I = 0; I < material_response_num_dof; ++I ){

                    for ( unsigned int k = 0; k < material_response_dim; ++k ){

                        for ( unsigned int a = 0; a < material_response_dim; ++a ){

                            *( dRdUMesh_begin + a ) -=
                                dRdpi_phase[ j ]
                                * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( interphasic_force_index + j ) + material_response_num_dof + material_response_dim * I + k ) )
                                * ( *( full_material_response_dof_gradient_begin + material_response_dim * I + a ) )
                                * ( *( interpolation_function_gradient_begin + k ) );

                        }

                    }

                }

            }

            // HEAT FLUX CONTRIBUTIONS
            for ( unsigned int j = 0; j < dim; ++j ){

                // density
                for ( auto p = std::pair< unsigned int, dRdRho_iter >( 0, dRdRho_begin ); p.second != dRdRho_end; ++p.first, ++p.second ){

                    *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + nphases * density_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + material_response_num_dof + material_response_dim * ( nphases * density_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // velocity
                for ( auto p = std::pair< unsigned int, dRdU_iter >( 0, dRdU_begin ); p.second != dRdU_end; ++p.first, ++p.second ){

                    *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + nphases * velocity_index + p.first ) ) * interpolation_function * dUDotdU;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + material_response_num_dof + material_response_dim * ( nphases * velocity_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) ) * dUDotdU;

                    }

                }

                // displacement
                for ( auto p = std::pair< unsigned int, dRdW_iter >( 0, dRdW_begin ); p.second != dRdW_end; ++p.first, ++p.second ){

                    *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + nphases * displacement_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + material_response_num_dof + material_response_dim * ( nphases * displacement_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // temperature
                for ( auto p = std::pair< unsigned int, dRdTheta_iter >( 0, dRdTheta_begin ); p.second != dRdTheta_end; ++p.first, ++p.second ){

                    *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + nphases * temperature_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + material_response_num_dof + material_response_dim * ( nphases * temperature_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // internal energy
                for ( auto p = std::pair< unsigned int, dRdE_iter >( 0, dRdE_begin ); p.second != dRdE_end; ++p.first, ++p.second ){

                    *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + nphases * internal_energy_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + material_response_num_dof + material_response_dim * ( nphases * internal_energy_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // additional dof
                for ( auto p = std::pair< unsigned int, dRdZ_iter >( 0, dRdZ_begin ); p.second != dRdZ_end; ++p.first, ++p.second ){

                    *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + nphases * additional_dof_index + p.first ) ) * interpolation_function;

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *p.second += dRdq_phase[ j ] * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + material_response_num_dof + material_response_dim * ( nphases * additional_dof_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                    }

                }

                // mesh displacement
                for ( unsigned int I = 0; I < material_response_num_dof; ++I ){

                    for ( unsigned int k = 0; k < material_response_dim; ++k ){

                        for ( unsigned int a = 0; a < material_response_dim; ++a ){

                            *( dRdUMesh_begin + a ) -=
                                dRdq_phase[ j ]
                                * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( heat_flux_index + j ) + material_response_num_dof + material_response_dim * I + k ) )
                                * ( *( full_material_response_dof_gradient_begin + material_response_dim * I + a ) )
                                * ( *( interpolation_function_gradient_begin + k ) );

                        }

                    }

                }

            }

            // INTERPHASTIC HEAT TRANSFER

        }

        template<
            int dim, int material_response_dim, int cauchy_stress_index, int internal_heat_generation_index, int heat_flux_index,
            int interphasic_force_index, int interphasic_heat_transfer_index,
            class density_iter, class density_dot_iter,
            class density_gradient_iter,
            class internal_energy_iter, class internal_energy_dot_iter,
            class internal_energy_gradient_iter,
            class velocity_iter, class velocity_gradient_iter,
            class material_response_iter,
            class volume_fraction_iter,
            typename test_function_type, class test_function_gradient_iter,
            class result_iter
        >
        void computeBalanceOfEnergy(
            const density_iter &density_begin, const density_iter &density_end,
            const density_dot_iter &density_dot_begin, const density_dot_iter &density_dot_end,
            const density_gradient_iter &density_gradient_begin, const density_gradient_iter &density_gradient_end,
            const internal_energy_iter &internal_energy_begin, const internal_energy_iter &internal_energy_end,
            const internal_energy_dot_iter &internal_energy_dot_begin, const internal_energy_dot_iter &internal_energy_dot_end,
            const internal_energy_gradient_iter &internal_energy_gradient_begin, const internal_energy_gradient_iter &internal_energy_gradient_end,
            const velocity_iter &velocity_begin, const velocity_iter &velocity_end,
            const velocity_gradient_iter &velocity_gradient_begin, const velocity_gradient_iter &velocity_gradient_end,
            const material_response_iter &material_response_begin, const material_response_iter &material_response_end,
            const volume_fraction_iter &volume_fraction_begin, const volume_fraction_iter &volume_fraction_end,
            const test_function_type &test_function,
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            result_iter result_begin, result_iter result_end
        ){
            /*!
             * Compute the full balance of energy in a variational context using a generalized material response vector for a multiphasic problem
             * 
             * \param &density_begin: The starting iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_end: The stopping iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot_begin: The starting iterator of the partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_dot_end: The stopping iterator of the partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_begin: The starting iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_end: The stopping iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot_begin: The starting iterator of the partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_dot_end: The stopping iterator of the partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &test_function: The value of the test function \f$ \left( \psi \right) \f$
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param result_begin: The starting iterator of the result of the non-divergence part of the balance of energy
             * \param result_end: The stopping iterator of the result of the non-divergence part of the balance of energy
             */

            const unsigned int nphases = ( unsigned int )( density_end - density_begin );
            const unsigned int material_response_size = ( unsigned int )( material_response_end - material_response_begin ) / nphases;

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( density_dot_end - density_dot_begin ), "The density and density dot vectors must be the same size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim == ( unsigned int )( density_gradient_end - density_gradient_begin ), "The density and density gradient vectors must be of consistent sizes"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( internal_energy_end - internal_energy_begin ), "The density and internal energy vectors must be the same size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( internal_energy_dot_end - internal_energy_dot_begin ), "The density and internal energy dot vectors must be the same size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim == ( unsigned int )( internal_energy_gradient_end - internal_energy_gradient_begin ), "The density and internal energy gradient vectors must be a consistent size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim == ( unsigned int )( velocity_end - velocity_begin ), "The density and velocity vectors must be a consistent size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * dim * dim == ( unsigned int )( velocity_gradient_end - velocity_gradient_begin ), "The density and velocity gradient vectors must be a consistent size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size == ( unsigned int )( material_response_end - material_response_begin ), "The density and material response vectors must be a consistent size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( volume_fraction_end - volume_fraction_begin ), "The density and volume fraction vectors must be the same size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                dim == ( unsigned int )( test_function_gradient_end - test_function_gradient_begin ), "The test function gradient must be of size dim"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( result_end - result_begin ), "The density and result vectors must be the same size"
            )

            for ( auto v = std::pair< unsigned int, density_iter >( 0, density_begin ); v.second != density_end; ++v.first, ++v.second ){

                computeBalanceOfEnergy<
                    dim, material_response_dim,
                    cauchy_stress_index, internal_heat_generation_index,
                    heat_flux_index, interphasic_force_index, interphasic_heat_transfer_index
                >
                (
                    *( density_begin + v.first ), *( density_dot_begin + v.first ),
                    density_gradient_begin + dim * v.first, density_gradient_begin + dim * ( v.first + 1 ),
                    *( internal_energy_begin + v.first ), *( internal_energy_dot_begin + v.first ),
                    internal_energy_gradient_begin + dim * v.first, internal_energy_gradient_begin + dim * ( v.first + 1 ),
                    velocity_begin + dim * v.first,                velocity_begin + dim * ( v.first + 1 ),
                    velocity_gradient_begin + dim * dim * v.first, velocity_gradient_begin + dim * dim * ( v.first + 1 ),
                    material_response_begin + material_response_size * v.first,
                    material_response_begin + material_response_size * ( v.first + 1 ),
                    *( volume_fraction_begin + v.first ),
                    test_function,
                    test_function_gradient_begin, test_function_gradient_end,
                    *( result_begin + v.first )
                );

            }

        }

        template<
            int dim, typename density_type, typename density_dot_type,
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
        ){
            /*!
             * Compute the non-divergence parts of the balance of energy i.e.
             * 
             * \f$ \frac{\partial}{\partial t}\left( \rho^{\alpha} e^{\alpha} \right) + \left( \rho^{\alpha} v_i^{\alpha} e^{\alpha} \right)_{,i} - \frac{1}{2} c^{\alpha} v_i^{\alpha} v_i^{\alpha} + \sum_{\beta} \pi_i^{\alpha \beta} v_i^{\alpha} - \phi^{\alpha}\sigma_{ji}^{\alpha}v_{i,j}^{\alpha} - \rho^{\alpha} r^{\alpha} \f$
             * 
             * \param &density: The apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot: The partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy: The internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot: The partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction: The volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation: The internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &result: The result of the non-divergence part of the balance of energy
             */
            
            using velocity_gradient_type = typename std::iterator_traits<velocity_gradient_iter>::value_type;

            // Compute the mass change rate 
            density_dot_type mass_change_rate;

            tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( density, density_dot, density_gradient_begin, density_gradient_end, velocity_begin, velocity_end,
                                                                                  velocity_gradient_begin, velocity_gradient_end, mass_change_rate );

            // Compute the trace of the velocity gradient
            velocity_gradient_type trace_velocity_gradient = 0.;
            for ( unsigned int i = 0; i < dim; i++ ){
                trace_velocity_gradient += *( velocity_gradient_begin + dim * i + i );
            }

            result = density_dot * internal_energy + density * internal_energy_dot
                   + internal_energy * std::inner_product( density_gradient_begin, density_gradient_end, velocity_begin, 0. )
                   + density * internal_energy * trace_velocity_gradient
                   + density * std::inner_product( velocity_begin, velocity_end, internal_energy_gradient_begin, 0. )
                   - 0.5 * mass_change_rate * std::inner_product( velocity_begin, velocity_end, velocity_begin, 0. )
                   + std::inner_product( net_interphase_force_begin, net_interphase_force_end, velocity_begin, 0. )
                   - density * internal_heat_generation;

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    result -= volume_fraction * ( *( cauchy_stress_begin + dim * j + i ) ) * ( *( velocity_gradient_begin + dim * i + j ) );

                }

            }

        }

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
        ){
            /*!
             * Compute the non-divergence parts of the balance of energy including the Jacobians i.e.
             * 
             * \f$ \frac{\partial}{\partial t}\left( \rho^{\alpha} e^{\alpha} \right) + \left( \rho^{\alpha} v_i^{\alpha} e^{\alpha} \right)_{,i} - \frac{1}{2} c^{\alpha} v_i^{\alpha} v_i^{\alpha} + \sum_{\beta} \pi_i^{\alpha \beta} v_i^{\alpha} - \phi^{\alpha}\sigma_{ji}^{\alpha}v_{i,j}^{\alpha} - \rho^{\alpha} r^{\alpha} \f$
             * 
             * \param &density: The apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot: The partial temporal derivative of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy: The internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot: The partial temporal derivative of the internal energy of phase \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction: The volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation: The internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &result: The result of the non-divergence part of the balance of energy
             * \param &dRdRho: The Jacobian w.r.t. the apparent density
             * \param &dRdRhoDot: The Jacobian w.r.t. the partial temporal derivative of the apparent density
             * \param &dRdGradRho_begin: The starting iterator of the Jacobian w.r.t. the spatial gradient of the apparent density
             * \param &dRdGradRho_end: The stopping iterator of the Jacobian w.r.t. the spatial gradient of the apparent density
             * \param &dRdE: The Jacobian w.r.t. the internal energy
             * \param &dRdEDot: The Jacobian w.r.t. the partial temporal derivative of the internal energy
             * \param &dRdGradE_begin: The starting iterator of the Jacobian w.r.t. the spatial gradient of the internal energy
             * \param &dRdGradE_end: The stopping iterator of the Jacobian w.r.t. the spatial gradient of the internal energy
             * \param &dRdV_begin: The starting iterator of the Jacobian w.r.t. the velocity
             * \param &dRdV_end: The stopping iterator of the Jacobian w.r.t. the velocity
             * \param &dRdGradV_begin: The starting iterator of the Jacobian w.r.t. the velocity gradient
             * \param &dRdGradV_end: The stopping iterator of the Jacobian w.r.t. the velocity gradient
             * \param &dRdCauchy_begin: The starting iterator of the Jacobian w.r.t. the Cauchy stress
             * \param &dRdCauchy_end: The stopping iterator of the Jacobian w.r.t. the Cauchy stress
             * \param &dRdVolumeFraction: The Jacobian w.r.t. the volume fraction
             * \param &dRdr: The Jacobian w.r.t. the internal heat generation
             * \param &dRdpi_begin: The starting iterator of the Jacobian w.r.t. the net interphase force
             * \param &dRdpi_end: The stopping iterator of the Jacobian w.r.t. the net interphase force
             */
            
            using velocity_gradient_type = typename std::iterator_traits<velocity_gradient_iter>::value_type;

            // Compute the mass change rate 
            density_dot_type mass_change_rate;

            density_dot_type dCdRho, dCdRhoDot;

            std::array< density_dot_type, dim > dCdGradRho, dCdV;

            std::array< density_dot_type, dim * dim > dCdGradV;

            tardigradeBalanceEquations::balanceOfMass::computeBalanceOfMass<dim>( density, density_dot, density_gradient_begin, density_gradient_end, velocity_begin, velocity_end,
                                                                                  velocity_gradient_begin, velocity_gradient_end, mass_change_rate,
                                                                                  dCdRho, dCdRhoDot, std::begin( dCdGradRho ), std::end( dCdGradRho ),
                                                                                  std::begin( dCdV ), std::end( dCdV ), std::begin( dCdGradV ), std::end( dCdGradV ) );

            // Compute the trace of the velocity gradient
            velocity_gradient_type trace_velocity_gradient = 0.;
            for ( unsigned int i = 0; i < dim; i++ ){
                trace_velocity_gradient += *( velocity_gradient_begin + dim * i + i );
            }

            velocity_gradient_type v_dot_v = std::inner_product( velocity_begin, velocity_end, velocity_begin, 0. );
            velocity_gradient_type grad_rho_dot_v = std::inner_product( density_gradient_begin, density_gradient_end, velocity_begin, 0. );

            result = density_dot * internal_energy + density * internal_energy_dot
                   + internal_energy * grad_rho_dot_v
                   + density * internal_energy * trace_velocity_gradient
                   + density * std::inner_product( velocity_begin, velocity_end, internal_energy_gradient_begin, 0. )
                   - 0.5 * mass_change_rate * v_dot_v
                   + std::inner_product( net_interphase_force_begin, net_interphase_force_end, velocity_begin, 0. )
                   - density * internal_heat_generation;

            dRdRho = internal_energy_dot + internal_energy * trace_velocity_gradient + std::inner_product( velocity_begin, velocity_end, internal_energy_gradient_begin, 0. )
                   - internal_heat_generation - 0.5 * dCdRho * v_dot_v;

            dRdRhoDot = internal_energy - 0.5 * dCdRhoDot * v_dot_v;

            std::transform( velocity_begin, velocity_end, dRdGradRho_begin, std::bind( std::multiplies<result_type>(), std::placeholders::_1, internal_energy ) );

            dRdE = density_dot + std::inner_product( density_gradient_begin, density_gradient_end, velocity_begin, 0. )
                 + density * trace_velocity_gradient;

            dRdEDot = density;

            std::transform( velocity_begin, velocity_end, dRdGradE_begin, std::bind( std::multiplies<result_type>(), std::placeholders::_1, density ) );

            std::fill( dRdV_begin, dRdV_end, 0. );

            std::fill( dRdGradV_begin, dRdGradV_end, 0. );

            std::fill( dRdCauchy_begin, dRdCauchy_end, 0. );

            dRdVolumeFraction = 0;

            dRdr = -density;

            std::copy( velocity_begin, velocity_end, dRdpi_begin );

            for ( unsigned int i = 0; i < dim; i++ ){

                *( dRdGradRho_begin + i ) -= 0.5 * dCdGradRho[ i ] * v_dot_v;

                *( dRdV_begin + i ) += internal_energy * ( *( density_gradient_begin + i ) )
                                     + density * ( *( internal_energy_gradient_begin + i ) )
                                     - 0.5 * dCdV[ i ] * v_dot_v
                                     - mass_change_rate * ( *( velocity_begin + i ) )
                                     + *( net_interphase_force_begin + i );

                *( dRdGradV_begin + dim * i + i ) += density * internal_energy;

                for ( unsigned int j = 0; j < dim; j++ ){

                    dRdVolumeFraction -= ( *( cauchy_stress_begin + dim * j + i ) ) * ( *( velocity_gradient_begin + dim * i + j ) );

                    *( dRdGradV_begin + dim * i + j ) -= 0.5 * dCdGradV[ dim * i + j ] * v_dot_v + volume_fraction * ( *( cauchy_stress_begin + dim * j + i ) );

                    *( dRdCauchy_begin + dim * i + j ) -= volume_fraction * ( *( velocity_gradient_begin + dim * j + i ) );

                }

            }

            result += volume_fraction * dRdVolumeFraction;

        }

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
        ){
            /*!
             * Compute the non-divergence parts of the balance of energy i.e.
             * 
             * \f$ \frac{\partial}{\partial t}\left( \rho^{\alpha} e^{\alpha} \right) + \left( \rho^{\alpha} v_i^{\alpha} e^{\alpha} \right)_{,i} - \frac{1}{2} c^{\alpha} v_i^{\alpha} v_i^{\alpha} + \sum_{\beta} \pi_i^{\alpha \beta} v_i^{\alpha} - \phi^{\alpha}\sigma_{ji}^{\alpha}v_{i,j}^{\alpha} - \rho^{\alpha} r^{\alpha} \f$
             * 
             * \param &density_begin: The starting iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_end: The stopping iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot_begin: The starting iterator of the partial temporal derivative of the apparent density (dm / dv) of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_dot_end: The partial temporal derivative of the apparent density (dm / dv) of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_begin: The starting iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_end: The stopping iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot_begin: The starting iterator of the partial temporal derivative of the internal energy of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_dot_end: The stopping iterator of the partial temporal derivative of the internal energy of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation_begin: The starting iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &internal_heat_generation_end: The stopping iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &result_begin: The starting iterator of the result of the non-divergence part of the balance of energy
             * \param &result_end: The stopping iterator of the result of the non-divergence part of the balance of energy
             */

            for ( auto rho = density_begin; rho != density_end; rho++ ){

                unsigned int phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfEnergyNonDivergence<dim>(
                    *( density_begin + phase ), *( density_dot_begin + phase ), density_gradient_begin + dim * phase, density_gradient_begin + dim * ( phase + 1 ),
                    *( internal_energy_begin + phase ), *( internal_energy_dot_begin + phase ),
                    internal_energy_gradient_begin + dim * phase, internal_energy_gradient_begin + dim * ( phase + 1 ),
                    velocity_begin + dim * phase, velocity_begin + dim * ( phase + 1 ),
                    velocity_gradient_begin + dim * dim * phase, velocity_gradient_begin + dim * dim * ( phase + 1 ),
                    cauchy_stress_begin + dim * dim * phase, cauchy_stress_begin + dim * dim * ( phase + 1 ),
                    *( volume_fraction_begin + phase ), *( internal_heat_generation_begin + phase ),
                    net_interphase_force_begin + dim * phase, net_interphase_force_begin + dim * ( phase + 1 ),
                    *( result_begin + phase )
                );

            }

        }

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
        ){
            /*!
             * Compute the non-divergence parts of the balance of energy i.e.
             * 
             * \f$ \frac{\partial}{\partial t}\left( \rho^{\alpha} e^{\alpha} \right) + \left( \rho^{\alpha} v_i^{\alpha} e^{\alpha} \right)_{,i} - \frac{1}{2} c^{\alpha} v_i^{\alpha} v_i^{\alpha} + \sum_{\beta} \pi_i^{\alpha \beta} v_i^{\alpha} - \phi^{\alpha}\sigma_{ji}^{\alpha}v_{i,j}^{\alpha} - \rho^{\alpha} r^{\alpha} \f$
             * 
             * \param &density_begin: The starting iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_end: The stopping iterator of the apparent density (dm / dv) of phase \f$ \alpha \f$ \f$\left(\rho^{\alpha}\right)\f$
             * \param &density_dot_begin: The starting iterator of the partial temporal derivative of the apparent density (dm / dv) of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_dot_end: The partial temporal derivative of the apparent density (dm / dv) of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} \rho^{\alpha}\right)\f$
             * \param &density_gradient_begin: The starting iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &density_gradient_end: The stopping iterator of the spatial gradient of the apparent density (dm/dv) of phase \f$ \alpha \f$ \f$\left( \rho^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_begin: The starting iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_end: The stopping iterator of the internal energy of phase \f$ \alpha \f$ \f$\left(e^{\alpha}\right)\f$
             * \param &internal_energy_dot_begin: The starting iterator of the partial temporal derivative of the internal energy of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_dot_end: The stopping iterator of the partial temporal derivative of the internal energy of phase
             *     \f$ \alpha \f$ \f$\left(\frac{\partial}{\partial t} e^{\alpha}\right)\f$
             * \param &internal_energy_gradient_begin: The starting iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &internal_energy_gradient_end: The stopping iterator of the spatial gradient of the internal energy of phase \f$ \alpha \f$ \f$\left( e^{\alpha}_{,i} \right) \f$
             * \param &velocity_begin: The starting iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_end: The stopping iterator of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_i \right) \f$
             * \param &velocity_gradient_begin: The starting iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &velocity_gradient_end: The stopping iterator of the spatial gradient of the velocity of phase \f$ \alpha \f$ \f$\left( v^{\alpha}_{i,j} \right) \f$
             * \param &cauchy_stress_begin: The starting iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &cauchy_stress_end: The stopping iterator of the true Cauchy stress \f$ \bar{\bf{\sigma}} \f$ where \f$ \bf{\sigma} = \phi \bar{\bf{\sigma}} \f$
             * \param &volume_fraction_begin: The starting iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &volume_fraction_end: The stopping iterator of the volume fraction of phase \f$ \alpha \f$ \f$ \left(\phi^{\alpha}\right) \f$
             * \param &internal_heat_generation_begin: The starting iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &internal_heat_generation_end: The stopping iterator of the internal heat generation per unit mass of phase \f$ \alpha \f$ \f$\left( r^{\alpha} \right)\f$
             * \param &net_interphase_force_begin: The starting iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &net_interphase_force_end: The stopping iterator of the net interphase force acting on phase \f$ \alpha \f$ \f$\left( \sum_{\beta} \phi^{\alpha \beta}_i \right) \f$
             * \param &result_begin: The starting iterator of the result of the non-divergence part of the balance of energy
             * \param &result_end: The stopping iterator of the result of the non-divergence part of the balance of energy
             * \param &dRdRho_begin: The starting iterator of the Jacobian w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the Jacobian w.r.t. the apparent density
             * \param &dRdRhoDot_begin: The starting iterator of the Jacobian w.r.t. the partial temporal derivative of the apparent density
             * \param &dRdRhoDot_end: The stopping iterator of the Jacobian w.r.t. the partial temporal derivative of the apparent density
             * \param &dRdGradRho_begin: The starting iterator of the Jacobian w.r.t. the spatial gradient of the apparent density
             * \param &dRdGradRho_end: The stopping iterator of the Jacobian w.r.t. the spatial gradient of the apparent density
             * \param &dRdE_begin: The starting iterator of the Jacobian w.r.t. the internal energy
             * \param &dRdE_end: The stopping iterator of the Jacobian w.r.t. the internal energy
             * \param &dRdEDot_begin: The starting iterator of the Jacobian w.r.t. the partial temporal derivative of the internal energy
             * \param &dRdEDot_end: The stopping iterator of the Jacobian w.r.t. the partial temporal derivative of the internal energy
             * \param &dRdGradE_begin: The starting iterator of the Jacobian w.r.t. the spatial gradient of the internal energy
             * \param &dRdGradE_end: The stopping iterator of the Jacobian w.r.t. the spatial gradient of the internal energy
             * \param &dRdV_begin: The starting iterator of the Jacobian w.r.t. the velocity
             * \param &dRdV_end: The stopping iterator of the Jacobian w.r.t. the velocity
             * \param &dRdGradV_begin: The starting iterator of the Jacobian w.r.t. the velocity gradient
             * \param &dRdGradV_end: The stopping iterator of the Jacobian w.r.t. the velocity gradient
             * \param &dRdCauchy_begin: The starting iterator of the Jacobian w.r.t. the Cauchy stress
             * \param &dRdCauchy_end: The stopping iterator of the Jacobian w.r.t. the Cauchy stress
             * \param &dRdVolumeFraction_begin: The starting iterator of the Jacobian w.r.t. the volume fraction
             * \param &dRdVolumeFraction_end: The stopping iterator of the Jacobian w.r.t. the volume fraction
             * \param &dRdr_begin: The starting iterator of the Jacobian w.r.t. the internal heat generation
             * \param &dRdr_end: The stopping iterator of the Jacobian w.r.t. the internal heat generation
             * \param &dRdpi_begin: The starting iterator of the Jacobian w.r.t. the net interphase force
             * \param &dRdpi_end: The stopping iterator of the Jacobian w.r.t. the net interphase force
             */

            for ( auto rho = density_begin; rho != density_end; rho++ ){

                unsigned int phase = ( unsigned int )( rho - density_begin );

                computeBalanceOfEnergyNonDivergence<dim>(
                    *( density_begin + phase ), *( density_dot_begin + phase ), density_gradient_begin + dim * phase, density_gradient_begin + dim * ( phase + 1 ),
                    *( internal_energy_begin + phase ), *( internal_energy_dot_begin + phase ),
                    internal_energy_gradient_begin + dim * phase, internal_energy_gradient_begin + dim * ( phase + 1 ),
                    velocity_begin + dim * phase, velocity_begin + dim * ( phase + 1 ),
                    velocity_gradient_begin + dim * dim * phase, velocity_gradient_begin + dim * dim * ( phase + 1 ),
                    cauchy_stress_begin + dim * dim * phase, cauchy_stress_begin + dim * dim * ( phase + 1 ),
                    *( volume_fraction_begin + phase ), *( internal_heat_generation_begin + phase ),
                    net_interphase_force_begin + dim * phase, net_interphase_force_begin + dim * ( phase + 1 ),
                    *( result_begin + phase ),
                    *( dRdRho_begin + phase ), *( dRdRhoDot_begin + phase ), dRdGradRho_begin + dim * phase, dRdGradRho_begin + dim * ( phase + 1 ),
                    *( dRdE_begin + phase ), *( dRdEDot_begin + phase ), dRdGradE_begin + dim * phase, dRdGradE_begin + dim * ( phase + 1 ),
                    dRdV_begin + dim * phase, dRdV_begin + dim * ( phase + 1 ),
                    dRdGradV_begin + dim * dim * phase, dRdGradV_begin + dim * dim * ( phase + 1 ),
                    dRdCauchy_begin + dim * dim * phase, dRdCauchy_begin + dim * dim * ( phase + 1 ),
                    *( dRdVolumeFraction_begin + phase ), *( dRdr_begin + phase ),
                    dRdpi_begin + dim * phase, dRdpi_begin + dim * ( phase + 1 )
                );

            }

        }

        template<
            int dim,
            class test_function_gradient_iter, class heat_flux_iter,
            typename result_type
        >
        void computeBalanceOfEnergyDivergence(
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            const heat_flux_iter &heat_flux_begin,                           const heat_flux_iter &heat_flux_end,
            result_type &result
        ){
            /*!
             * Compute the divergence part of the balance of energy for a variationally based method i.e.:
             * 
             * \f$ \psi q_{i,i} = \left( \psi q_i \right)_{,i} - \psi_{,i} q_i \f$
             *
             * where we will compute the second term as the first term is the heat flux boundary condition.
             * 
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &result: The resulting divergence term
             */

            result = -std::inner_product( test_function_gradient_begin, test_function_gradient_end, heat_flux_begin, 0. );

        }

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
            dRdGradTestFunction_iter dRdGradTestFunction_begin, dRdGradTestFunction_iter dRdGradTestFunction_end,
            dRdq_iter dRdq_begin,                               dRdq_iter dRdq_end
        ){
            /*!
             * Compute the divergence part of the balance of energy for a variationally based method i.e.:
             * 
             * \f$ \psi q_{i,i} = \left( \psi q_i \right)_{,i} - \psi_{,i} q_i \f$
             *
             * where we will compute the second term as the first term is the heat flux boundary condition.
             * 
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &result: The resulting divergence term
             * \param &dRdGradTestFunction_begin: The starting iterator of the Jacobian of the result w.r.t. the test function gradient
             * \param &dRdGradTestFunction_end: The stopping iterator of the Jacobian of the result w.r.t. the test function gradient
             * \param &dRdq_begin: The starting iterator of the Jacobian of the result w.r.t. the heat flux
             * \param &dRdq_end: The stopping iterator of the Jacobian of the result w.r.t. the heat flux
             */

            using test_function_gradient_type = typename std::iterator_traits<test_function_gradient_iter>::value_type;

            result = -std::inner_product( test_function_gradient_begin, test_function_gradient_end, heat_flux_begin, 0. );

            std::transform( heat_flux_begin, heat_flux_end, dRdGradTestFunction_begin, std::bind( std::multiplies<result_type>(), std::placeholders::_1, -1. ) );

            std::transform( test_function_gradient_begin, test_function_gradient_end, dRdq_begin, std::bind( std::multiplies<test_function_gradient_type>(), std::placeholders::_1, -1. ) );

        }

        template<
            int dim,
            class test_function_gradient_iter, class heat_flux_iter,
            class result_iter
        >
        void computeBalanceOfEnergyDivergence(
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            const heat_flux_iter &heat_flux_begin,                           const heat_flux_iter &heat_flux_end,
            result_iter result_begin,                                        result_iter result_end
        ){
            /*!
             * Compute the multiphase divergence part of the balance of energy for a variationally based method i.e.:
             * 
             * \f$ \psi q_{i,i} = \left( \psi q_i \right)_{,i} - \psi_{,i} q_i \f$
             *
             * where we will compute the second term as the first term is the heat flux boundary condition.
             * 
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &result_begin: The starting iterator of the resulting divergence term
             * \param &result_end: The stopping iterator of the resulting divergence term
             */

            for ( auto r = result_begin; r != result_end; r++ ){

                const unsigned int phase = ( unsigned int )( r - result_begin );

                computeBalanceOfEnergyDivergence<dim>(
                    test_function_gradient_begin,  test_function_gradient_end,
                    heat_flux_begin + dim * phase, heat_flux_begin + dim * ( phase + 1 ),
                    *r
                );

            }

        }

        template<
            int dim,
            class test_function_gradient_iter, class heat_flux_iter,
            class result_iter,
            class dRdGradTestFunction_iter, class dRdq_iter
        >
        void computeBalanceOfEnergyDivergence(
            const test_function_gradient_iter &test_function_gradient_begin, const test_function_gradient_iter &test_function_gradient_end,
            const heat_flux_iter &heat_flux_begin,                           const heat_flux_iter &heat_flux_end,
            result_iter result_begin, result_iter result_end,
            dRdGradTestFunction_iter dRdGradTestFunction_begin, dRdGradTestFunction_iter dRdGradTestFunction_end,
            dRdq_iter dRdq_begin,             dRdq_iter dRdq_end
        ){
            /*!
             * Compute the multiphase divergence part of the balance of energy for a variationally based method i.e.:
             * 
             * \f$ \psi q_{i,i} = \left( \psi q_i \right)_{,i} - \psi_{,i} q_i \f$
             *
             * where we will compute the second term as the first term is the heat flux boundary condition.
             * 
             * \param &test_function_gradient_begin: The starting iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &test_function_gradient_end: The stopping iterator of the spatial gradient of the test function \f$ \left( \psi_{,i} \right) \f$
             * \param &heat_flux_begin: The starting iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &heat_flux_end: The stopping iterator of the heat flux vector \f$ \left( q_i \right) \f$
             * \param &result_begin: The starting iterator of the resulting divergence term
             * \param &result_end: The stopping iterator of the resulting divergence term
             * \param &dRdGradTestFunction_begin: The starting iterator of the Jacobian of the result w.r.t. the test function gradient
             * \param &dRdGradTestFunction_end: The stopping iterator of the Jacobian of the result w.r.t. the test function gradient
             * \param &dRdq_begin: The starting iterator of the Jacobian of the result w.r.t. the heat flux
             * \param &dRdq_end: The stopping iterator of the Jacobian of the result w.r.t. the heat flux
             */

            for ( auto r = result_begin; r != result_end; r++ ){

                const unsigned int phase = ( unsigned int )( r - result_begin );

                computeBalanceOfEnergyDivergence<dim>(
                    test_function_gradient_begin,   test_function_gradient_end,
                    heat_flux_begin + dim * phase,  heat_flux_begin + dim * ( phase + 1 ),
                    *r,
                    dRdGradTestFunction_begin + dim * phase, dRdGradTestFunction_begin + dim * ( phase + 1 ),
                    dRdq_begin + dim * phase,       dRdq_begin + dim * ( phase + 1 )
                );

            }

        }

    }

}
