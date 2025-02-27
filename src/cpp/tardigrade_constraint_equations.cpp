/**
  ******************************************************************************
  * \file tardigrade_constraint_equations.cpp
  ******************************************************************************
  * The source file for the constraint equations associated with the balance
  * equations
  ******************************************************************************
  */

#include "tardigrade_constraint_equations.h"
#include<numeric>
#include<algorithm>
#include<functional>
#include<array>

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
        ){
            /*!
             * Compute the value of the constraint on the internal energy.
             * 
             * Traditionally, the internal energy is related to the temperature via
             * 
             * \f$ e = c \theta \f$
             * 
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             * 
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             * 
             * \param &internal_energy: The internal energy degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &result: The resulting error between the internal energy DOF and the material response
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(
                predicted_internal_energy_index < ( unsigned int )( material_response_end - material_response_begin ),
                "The index for the predicted internal energy is out of range for the material response"
            )

            result = ( *( material_response_begin + predicted_internal_energy_index ) - internal_energy ) * test_function;

        }

        template<
            int material_response_dim, int predicted_internal_energy_index,
            int material_response_num_dof,
            typename internal_energy_type,
            class material_response_iter,
            class material_response_jacobian_iter,
            typename test_function_type,
            typename interpolation_function_type,
            class interpolation_function_gradient_iter,
            class full_material_response_dof_gradient_iter,
            typename dUDotdU_type,
            typename result_type,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter,
            class dRdTheta_iter, class dRdE_iter, class dRdZ_iter,
            class dRdUMesh_iter,
            int density_index        ,
            int displacement_index   ,
            int velocity_index       ,
            int temperature_index    ,
            int internal_energy_index,
            int additional_dof_index 
        >
        void computeInternalEnergyConstraint(
            const internal_energy_type &internal_energy,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU,
            const unsigned int &phase,
            result_type &result,
            dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end,
            dRdE_iter     dRdE_begin,     dRdE_iter     dRdE_end,
            dRdZ_iter     dRdZ_begin,     dRdZ_iter     dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end
        ){
            /*!
             * Compute the value of the constraint on the internal energy.
             * 
             * Traditionally, the internal energy is related to the temperature via
             * 
             * \f$ e = c \theta \f$
             * 
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             * 
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             * 
             * \param &internal_energy: The internal energy degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &interpolation_function: The current value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the current value of the spatial gradient of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the current value of the spatial gradient of the interpolation function
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the material response dof vector
             * \param &dUDotdU: The derivative of the phase velocity w.r.t. the phase displacement dof
             * \param &phase: The current active phase
             * \param &result: The resulting error between the internal energy DOF and the material response
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdU_begin: The starting iterator of the derivative of the residual w.r.t. the phase spatial DOF associated with velocity
             * \param &dRdU_end: The stopping iterator of the derivative of the residual w.r.t. the phase spatial DOF associated with velocity
             * \param &dRdW_begin: The starting iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdW_end: The stopping iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdTheta_begin: The starting iterator of the derivative of the residual w.r.t. the phase temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the residual w.r.t. the phase temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the residual w.r.t. the phase internal energy
             * \param &dRdE_end: The stopping iterator of the derivative of the residual w.r.t. the phase internal energy
             * \param &dRdZ_begin: The starting iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdZ_end: The stopping iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh displacement
             */

            const unsigned int nphases = ( unsigned int )( dRdRho_end - dRdRho_begin );

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_dim == ( unsigned int )( dRdU_end - dRdU_begin ),
                "dRdU must be a consistent size with the material response dimension and the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_dim == ( unsigned int )( dRdW_end - dRdW_begin ),
                "dRdW must be a consistent size with the material response dimension and the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( dRdTheta_end - dRdTheta_begin ),
                "dRdTheta must be the same size as the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( dRdE_end - dRdE_begin ),
                "dRdE must be the same size as the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                material_response_dim == ( unsigned int )( dRdUMesh_end - dRdUMesh_begin ),
                "dRdUMesh must be the same size as the spatial dimension of the material response"
            )

            std::fill( dRdRho_begin,   dRdRho_end,   0 );
            std::fill( dRdU_begin,     dRdU_end,     0 );
            std::fill( dRdW_begin,     dRdW_end,     0 );
            std::fill( dRdTheta_begin, dRdTheta_end, 0 );
            std::fill( dRdE_begin,     dRdE_end,     0 );
            std::fill( dRdZ_begin,     dRdZ_end,     0 );
            std::fill( dRdUMesh_begin, dRdUMesh_end, 0 );

            computeInternalEnergyConstraint<
                predicted_internal_energy_index
            >(
                internal_energy,
                material_response_begin, material_response_end,
                test_function, result
            );

            *( dRdE_begin + phase ) -= test_function * interpolation_function;

            // PREDICTED INTERNAL ENERGY
            // density
            for ( auto p = std::pair< unsigned int, dRdRho_iter >( 0, dRdRho_begin ); p.second != dRdRho_end; ++p.first, ++p.second ){

                *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + nphases * density_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + material_response_num_dof + material_response_dim * ( nphases * density_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // velocity
            for ( auto p = std::pair< unsigned int, dRdU_iter >( 0, dRdU_begin ); p.second != dRdU_end; ++p.first, ++p.second ){

                *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + nphases * velocity_index + p.first ) ) * interpolation_function * dUDotdU;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + material_response_num_dof + material_response_dim * ( nphases * velocity_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) ) * dUDotdU;

                }

            }

            // displacement
            for ( auto p = std::pair< unsigned int, dRdW_iter >( 0, dRdW_begin ); p.second != dRdW_end; ++p.first, ++p.second ){

                *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + nphases * displacement_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + material_response_num_dof + material_response_dim * ( nphases * displacement_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // temperature
            for ( auto p = std::pair< unsigned int, dRdTheta_iter >( 0, dRdTheta_begin ); p.second != dRdTheta_end; ++p.first, ++p.second ){

                *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + nphases * temperature_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + material_response_num_dof + material_response_dim * ( nphases * temperature_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // internal energy
            for ( auto p = std::pair< unsigned int, dRdE_iter >( 0, dRdE_begin ); p.second != dRdE_end; ++p.first, ++p.second ){

                *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + nphases * internal_energy_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + material_response_num_dof + material_response_dim * ( nphases * internal_energy_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // additional dof
            for ( auto p = std::pair< unsigned int, dRdZ_iter >( 0, dRdZ_begin ); p.second != dRdZ_end; ++p.first, ++p.second ){

                *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + nphases * additional_dof_index + p.first ) ) * interpolation_function;

                for ( unsigned int a = 0; a < material_response_dim; ++a ){

                    *p.second += test_function * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + material_response_num_dof + material_response_dim * ( nphases * additional_dof_index + p.first ) + a ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

            // mesh displacement
            for ( unsigned int I = 0; I < material_response_num_dof; ++I ){

                for ( unsigned int k = 0; k < material_response_dim; ++k ){

                    for ( unsigned int a = 0; a < material_response_dim; ++a ){

                        *( dRdUMesh_begin + a ) -=
                            test_function
                            * ( *( material_response_jacobian_begin + material_response_num_dof * ( 1 + material_response_dim ) * ( predicted_internal_energy_index ) + material_response_num_dof + material_response_dim * I + k ) )
                            * ( *( full_material_response_dof_gradient_begin + material_response_dim * I + a ) )
                            * ( *( interpolation_function_gradient_begin + k ) );

                    }

                }

            }

            for ( unsigned int a = 0; a < material_response_dim; ++a ){

                *( dRdUMesh_begin + a ) += result * ( *( interpolation_function_gradient_begin + a ) );

            }

        }

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
        ){
            /*!
             * Compute the value of the constraint on the internal energy for a multiphase problem.
             * 
             * Traditionally, the internal energy is related to the temperature via
             * 
             * \f$ e = c \theta \f$
             * 
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             * 
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             * 
             * \param &internal_energy_begin: The starting iterator of the internal energy degree of freedom
             * \param &internal_energy_end: The stopping iterator of the internal energy degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &result_begin: The starting iterator of the resulting error between the internal energy DOF and the material response
             * \param &result_end: The starting iterator of the resulting error between the internal energy DOF and the material response
             */

            const unsigned int nphases = ( unsigned int )( result_end - result_begin );
            const unsigned int material_response_size = ( unsigned int )( material_response_end - material_response_begin ) / nphases;

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( internal_energy_end - internal_energy_begin ),
                "The number of internal energy values must be equal to the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size == ( unsigned int )( material_response_end - material_response_begin ),
                "The material response vector must be a scalar multiple of the number of phases"
            )

            for ( auto v = std::pair< unsigned int, result_iter >( 0, result_begin ); v.second != result_end; ++v.first, ++v.second ){

                computeInternalEnergyConstraint<
                    predicted_internal_energy_index
                >
                (
                    *( internal_energy_begin + v.first ),
                    material_response_begin + material_response_size * v.first,
                    material_response_begin + material_response_size * ( v.first + 1 ),
                    test_function,
                    *v.second
                );

            }

        }

        template<
            int material_response_dim, int predicted_internal_energy_index,
            int material_response_num_dof,
            class internal_energy_iter,
            class material_response_iter,
            class material_response_jacobian_iter,
            typename test_function_type,
            typename interpolation_function_type,
            class interpolation_function_gradient_iter,
            class full_material_response_dof_gradient_iter,
            typename dUDotdU_type,
            class result_iter,
            class dRdRho_iter, class dRdU_iter, class dRdW_iter,
            class dRdTheta_iter, class dRdE_iter, class dRdZ_iter,
            class dRdUMesh_iter,
            int density_index        ,
            int displacement_index   ,
            int velocity_index       ,
            int temperature_index    ,
            int internal_energy_index,
            int additional_dof_index 
        >
        void computeInternalEnergyConstraint(
            const internal_energy_iter &internal_energy_begin,
            const internal_energy_iter &internal_energy_end,
            const material_response_iter &material_response_begin,
            const material_response_iter &material_response_end,
            const material_response_jacobian_iter &material_response_jacobian_begin,
            const material_response_jacobian_iter &material_response_jacobian_end,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_begin,
            const full_material_response_dof_gradient_iter &full_material_response_dof_gradient_end,
            const dUDotdU_type dUDotdU,
            result_iter result_begin,
            result_iter result_end,
            dRdRho_iter dRdRho_begin, dRdRho_iter dRdRho_end,
            dRdU_iter dRdU_begin, dRdU_iter dRdU_end,
            dRdW_iter dRdW_begin, dRdW_iter dRdW_end,
            dRdTheta_iter dRdTheta_begin, dRdTheta_iter dRdTheta_end,
            dRdE_iter     dRdE_begin,     dRdE_iter     dRdE_end,
            dRdZ_iter     dRdZ_begin,     dRdZ_iter     dRdZ_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end
        ){
            /*!
             * Compute the value of the constraint on the internal energy.
             * 
             * Traditionally, the internal energy is related to the temperature via
             * 
             * \f$ e = c \theta \f$
             * 
             * where \f$ e \f$ is the internal energy, \f$ c \f$ is the specific heat,
             * and \f$ \theta \f$ is the temperature. The specific heat may not be a
             * constant value however and may depend on the temperature, pressure, or
             * other measures of the material state. By defining a constraint as
             * 
             * \f$ R = e^{p} - e \f$ where \f$ e^p \f$ is the predicted internal energy
             * from the material response calculation.
             * 
             * \param &internal_energy_begin: The starting iterator of the internal energy degree of freedom
             * \param &internal_energy_end: The stopping iterator internal energy degree of freedom
             * \param &material_response_begin: The starting iterator of the material response vector
             * \param &material_response_end: The stopping iterator of the material response vector
             * \param &material_response_jacobian_begin: The starting iterator of the material response Jacobian vector
             * \param &material_response_jacobian_end: The stopping iterator of the material response Jacobian vector
             * \param &test_function: The test function from the variational solution strategy
             * \param &interpolation_function: The current value of the interpolation function
             * \param &interpolation_function_gradient_begin: The starting iterator of the current value of the spatial gradient of the interpolation function
             * \param &interpolation_function_gradient_end: The stopping iterator of the current value of the spatial gradient of the interpolation function
             * \param &full_material_response_dof_gradient_begin: The starting iterator of the spatial gradient of the material response dof vector
             * \param &full_material_response_dof_gradient_end: The stopping iterator of the spatial gradient of the material response dof vector
             * \param &dUDotdU: The derivative of the phase velocity w.r.t. the phase displacement dof
             * \param &result_begin: The starting iterator of the resulting error between the internal energy DOF and the material response
             * \param &result_end: The stopping starting iterator of the resulting error between the internal energy DOF and the material response
             * \param &dRdRho_begin: The starting iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdRho_end: The stopping iterator of the derivative of the residual w.r.t. the apparent density
             * \param &dRdU_begin: The starting iterator of the derivative of the residual w.r.t. the phase spatial DOF associated with velocity
             * \param &dRdU_end: The stopping iterator of the derivative of the residual w.r.t. the phase spatial DOF associated with velocity
             * \param &dRdW_begin: The starting iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdW_end: The stopping iterator of the derivative of the residual w.r.t. the phase displacement dof
             * \param &dRdTheta_begin: The starting iterator of the derivative of the residual w.r.t. the phase temperature
             * \param &dRdTheta_end: The stopping iterator of the derivative of the residual w.r.t. the phase temperature
             * \param &dRdE_begin: The starting iterator of the derivative of the residual w.r.t. the phase internal energy
             * \param &dRdE_end: The stopping iterator of the derivative of the residual w.r.t. the phase internal energy
             * \param &dRdZ_begin: The starting iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdZ_end: The stopping iterator of the derivative of the residual w.r.t. the phase additional dof
             * \param &dRdUMesh_begin: The starting iterator of the derivative of the residual w.r.t. the mesh displacement
             * \param &dRdUMesh_end: The stopping iterator of the derivative of the residual w.r.t. the mesh displacement
             */

            using internal_energy_type     = typename std::iterator_traits<internal_energy_iter>::value_type;
            using result_type              = typename std::iterator_traits<result_iter>::value_type;

            const unsigned int nphases = ( unsigned int )( result_end - result_begin );
            const unsigned int material_response_size = ( unsigned int )( material_response_end - material_response_begin ) / nphases;
            const unsigned int num_additional_dof = ( unsigned int )( dRdZ_end - dRdZ_begin ) / nphases;

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases == ( unsigned int )( internal_energy_end - internal_energy_begin ),
                "The number of internal energy values must be equal to the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size == ( unsigned int )( material_response_end - material_response_begin ),
                "The material response vector must be a scalar multiple of the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_size * ( 1 + material_response_dim ) * material_response_num_dof == ( unsigned int )( material_response_jacobian_end - material_response_jacobian_begin ),
                "The material response jacobian vector must be a scalar multiple of the number of phases, the material response size and 1 + the material response dimension times the number of dof in the material response"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                material_response_dim * material_response_num_dof == ( unsigned int )( full_material_response_dof_gradient_end - full_material_response_dof_gradient_begin ),
                "The full material response dof gradient have a size of the material response dimension times the number of dof in the material resionse"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * material_response_dim == ( unsigned int )( dRdU_end - dRdU_begin ),
                "dRdU must be a consistent size with the material response dimension and the number of phases"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases * material_response_dim == ( unsigned int )( dRdW_end - dRdW_begin ),
                "dRdW must be a consistent size with the material response dimension and the number of phases squared"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases == ( unsigned int )( dRdTheta_end - dRdTheta_begin ),
                "dRdTheta must be the same size as the number of phases squared"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * nphases == ( unsigned int )( dRdE_end - dRdE_begin ),
                "dRdE must be the same size as the number of phases squared"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                nphases * material_response_dim == ( unsigned int )( dRdUMesh_end - dRdUMesh_begin ),
                "dRdUMesh must be the same size as the spatial dimension of the material response times the number of phases"
            )

            for ( auto v = std::pair< unsigned int, result_iter >( 0, result_begin ); v.second != result_end; ++v.first, ++v.second ){

                computeInternalEnergyConstraint<
                    material_response_dim, predicted_internal_energy_index,
                    material_response_num_dof,
                    internal_energy_type,
                    material_response_iter,
                    material_response_jacobian_iter,
                    test_function_type,
                    interpolation_function_type,
                    interpolation_function_gradient_iter,
                    full_material_response_dof_gradient_iter,
                    dUDotdU_type,
                    result_type,
                    dRdRho_iter, dRdU_iter, dRdW_iter,
                    dRdTheta_iter, dRdE_iter, dRdZ_iter,
                    dRdUMesh_iter,
                    density_index        ,
                    displacement_index   ,
                    velocity_index       ,
                    temperature_index    ,
                    internal_energy_index,
                    additional_dof_index 
                >
                (
                    *( internal_energy_begin + v.first ),
                    material_response_begin + material_response_size * v.first,
                    material_response_begin + material_response_size * ( v.first + 1 ),
                    material_response_jacobian_begin + material_response_size * material_response_num_dof * ( 1 + material_response_dim ) * v.first,
                    material_response_jacobian_begin + material_response_size * material_response_num_dof * ( 1 + material_response_dim ) * ( v.first + 1 ),
                    test_function,
                    interpolation_function,
                    interpolation_function_gradient_begin,
                    interpolation_function_gradient_end,
                    full_material_response_dof_gradient_begin,
                    full_material_response_dof_gradient_end,
                    dUDotdU,
                    v.first,
                    *v.second,
                    dRdRho_begin   +                                 nphases * v.first, dRdRho_begin   +                                 nphases * ( v.first + 1 ),
                    dRdU_begin     +         nphases * material_response_dim * v.first, dRdU_begin     +         nphases * material_response_dim * ( v.first + 1 ),
                    dRdW_begin     +         nphases * material_response_dim * v.first, dRdW_begin     +         nphases * material_response_dim * ( v.first + 1 ),
                    dRdTheta_begin +                                 nphases * v.first, dRdTheta_begin +                                 nphases * ( v.first + 1 ),
                    dRdE_begin     +                                 nphases * v.first, dRdE_begin     +                                 nphases * ( v.first + 1 ),
                    dRdZ_begin     +                      num_additional_dof * v.first, dRdZ_begin     +                      num_additional_dof * ( v.first + 1 ),
                    dRdUMesh_begin +                   material_response_dim * v.first, dRdUMesh_begin +                   material_response_dim * ( v.first + 1 )
                );

            }

        }

        template<
            class displacement_dot_iter,
            class velocity_iter,
            typename test_function_type,
            class result_iter
        >
        void computeDisplacementConstraint(
            const displacement_dot_iter &displacement_dot_begin, const displacement_dot_iter &displacement_dot_end,
            const velocity_iter         &velocity_begin,         const velocity_iter         &velocity_end,
            const test_function_type &test_function,
            result_iter result_begin, result_iter result_end
        ){
            /*!
             * Compute the constraint for the displacement to the velocity DOF. If this is integrated over the volume
             * then the test function should be provided as normal. If this is directly applied at the nodes, then
             * setting the test function to 1 will achieve the desired result.
             * 
             * \param displacement_dot_begin: The starting iterator for the time derivative of the displacement dof
             * \param displacement_dot_end: The stopping iterator for the time derivative of the displacement dof
             * \param velocity_begin: The starting iterator for the velocity dof
             * \param velocity_end: The stopping iterator for the velocity dof
             * \param test_function: The test function. If the constraint is being applied in a volume-integrated
             *     way then this should be set to the test function value. If it is being applied at a mesh's nodes
             *     then it should be set to 1.
             * \param result_begin: The starting iterator of the constraint violation result
             * \param result_end: The stopping iterator of the constraint violation result
             */

            using result_type = typename std::iterator_traits<result_iter>::value_type;

            TARDIGRADE_ERROR_TOOLS_EVAL(
                const unsigned int length = ( unsigned int )( displacement_dot_end - displacement_dot_begin );
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                length == ( unsigned int )( velocity_end - velocity_begin ),
                "The velocity and density dot vectors must be the same size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                length == ( unsigned int )( result_end - result_begin ),
                "The result and density dot vectors must be the same size"
            )

            std::transform(
                displacement_dot_begin,
                displacement_dot_end,
                velocity_begin,
                result_begin,
                std::minus< result_type >( )
            );

            std::transform(
                result_begin, result_end, result_begin,
                std::bind(
                    std::multiplies< result_type >( ),
                    std::placeholders::_1,
                    test_function
                )
            );

        }

        template<
            int dim,
            class displacement_dot_iter,
            class velocity_iter,
            typename test_function_type,
            typename interpolation_function_type,
            class interpolation_function_gradient_iter,
            typename dDDotdD_type,
            class result_iter,
            class dRdD_iter, class dRdV_iter, class dRdUMesh_iter
        >
        void computeDisplacementConstraint(
            const displacement_dot_iter &displacement_dot_begin, const displacement_dot_iter &displacement_dot_end,
            const velocity_iter         &velocity_begin,         const velocity_iter         &velocity_end,
            const test_function_type &test_function,
            const interpolation_function_type &interpolation_function,
            const interpolation_function_gradient_iter &interpolation_function_gradient_begin,
            const interpolation_function_gradient_iter &interpolation_function_gradient_end,
            const dDDotdD_type &dDDotdD,
            result_iter result_begin, result_iter result_end,
            dRdD_iter dRdD_begin,         dRdD_iter dRdD_end,
            dRdV_iter dRdV_begin,         dRdV_iter dRdV_end,
            dRdUMesh_iter dRdUMesh_begin, dRdUMesh_iter dRdUMesh_end
        ){
            /*!
             * Compute the constraint for the displacement to the velocity DOF. If this is integrated over the volume
             * then the test function, interpolation function and interpolation function gradient should be provided
             * as normal. If this is directly applied at the nodes, then setting the test function and interpolation
             * function to 1 and the interpolation function gradient to zero will achieve the desired result.
             * 
             * \param displacement_dot_begin: The starting iterator for the time derivative of the displacement dof
             * \param displacement_dot_end: The stopping iterator for the time derivative of the displacement dof
             * \param velocity_begin: The starting iterator for the velocity dof
             * \param velocity_end: The stopping iterator for the velocity dof
             * \param test_function: The test function. If the constraint is being applied in a volume-integrated
             *     way then this should be set to the test function value. If it is being applied at a mesh's nodes
             *     then it should be set to 1.
             * \param interpolation_function: The interpolation function. If the constraint is being applied in a volume-integrated
             *     way then this should be set to the interpolation function value. If it is being applied at a mesh's nodes
             *     then it should be set to 1.
             * \param interpolation_function_gradient_begin: The starting iterator of the interpolation function. If the constraint is being applied in a volume-integrated
             *     way then this should be set to the interpolation function gradient value. If it is being applied at a mesh's nodes
             *     then it should be set to 0.
             * \param interpolation_function_gradient_end: The starting iterator of the interpolation function. If the constraint is being applied in a volume-integrated
             *     way then this should be set to the interpolation function gradient value. If it is being applied at a mesh's nodes
             *     then it should be set to 0.
             * \param dDDotdD: The derivative of the displacement dot w.r.t. the displacement
             * \param result_begin: The starting iterator of the constraint violation result
             * \param result_end: The stopping iterator of the constraint violation result
             * \param dRdD_begin: The starting iterator for the derivative of the residual w.r.t. the displacement (only the diagonal is returned as everything else is zero)
             * \param dRdD_end:   The stopping iterator for the derivative of the residual w.r.t. the displacement (only the diagonal is returned as everything else is zero)
             * \param dRdV_begin: The starting iterator for the derivative of the residual w.r.t. the velocity (only the diagonal is returned as everything else is zero)
             * \param dRdV_end:   The stopping iterator for the derivative of the residual w.r.t. the velocity (only the diagonal is returned as everything else is zero)
             * \param dRdUMesh_begin: The starting iterator for the derivative of the residual w.r.t. the mesh displacement
             * \param dRdUMesh_end: The stopping iterator for the derivative of the residual w.r.t. the mesh displacement
             */

            TARDIGRADE_ERROR_TOOLS_EVAL(
                const unsigned int length = ( unsigned int )( displacement_dot_end - displacement_dot_begin );
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                length == ( unsigned int )( dRdD_end - dRdD_begin ),
                "The dRdD and density dot vectors must be the same size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                length == ( unsigned int )( dRdV_end - dRdV_begin ),
                "The dRdV and density dot vectors must be the same size"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                length * dim == ( unsigned int )( dRdUMesh_end - dRdUMesh_begin ),
                "The dRdUMesh vector must be the length of the density dot vector time the length of the interpolation function gradient vector"
            )

            computeDisplacementConstraint(
                displacement_dot_begin, displacement_dot_end,
                velocity_begin,         velocity_end,
                test_function,
                result_begin, result_end
            );

            std::fill( dRdD_begin, dRdD_end, 0 );
            std::fill( dRdV_begin, dRdV_end, 0 );
            std::fill( dRdUMesh_begin, dRdUMesh_end, 0 );

            for ( auto v = std::pair< unsigned int, dRdD_iter >( 0, dRdD_begin ); v.second != dRdD_end; ++v.first, ++v.second ){
                *v.second                 =  test_function * dDDotdD * interpolation_function;
                *( dRdV_begin + v.first ) = -test_function * interpolation_function;

                for ( unsigned int a = 0; a < dim; ++a ){

                    *( dRdUMesh_begin + dim * v.first + a ) = ( *( result_begin + v.first ) ) * ( *( interpolation_function_gradient_begin + a ) );

                }

            }

        }

    }

}
