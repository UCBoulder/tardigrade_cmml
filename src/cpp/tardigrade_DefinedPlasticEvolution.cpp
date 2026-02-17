/**
 *****************************************************************************
 * \file tardigrade_DefinedPlasticEvolution.cpp
 *****************************************************************************
 * A definition of a material where a permanent deformation configuration
 * evolves due to an externally defined velocity gradient
 *****************************************************************************
 */

#include "tardigrade_DefinedPlasticEvolution.h"

namespace tardigradeCMML {

    namespace DefinedPlasticEvolution {

        int DefinedPlasticEvolution::evaluate_model(const time_type &current_time, const time_type &dt,
                                                    const current_dof_type  *current_dof_begin,
                                                    const previous_dof_type *previous_dof_begin,
                                                    const unsigned int dof_size, const parameter_type *parameters_begin,
                                                    const unsigned int parameters_size, sdvs_type *sdvs_begin,
                                                    const unsigned int sdvs_size, result_type *result_begin,
                                                    const unsigned int result_size, std::string &output_message) {
            /*!
             * Evaluate the material model given the incoming data
             *
             * Returns:
             *     0: No errors. Solution converged.
             *     1: Convergence Error. Request timestep cutback.
             *     2: Fatal errors encountered. Terminate evaluation.
             *
             * \param &current_time: The current time of the simulation
             * \param &dt:           The change in time of the simulation
             * \param *current_dof_begin: A pointer to the starting element of the current values of the degrees of
             * freedom \param *previous_dof_begin: A pointer to the starting element of the previous values of the
             * degrees of freedom \param dof_size: The size of the dof arrays \param *parameters_begin: A pointer to the
             * starting element of the model parameters \param parameters_size: The size of the parameter vector \param
             * *sdvs_begin: A pointer to the starting element of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param sdvs_size: The size of the SDVS vector
             * \param *result_begin: A pointer to the starting element of the result vector.
             * \param result_size: The size of the result
             * \param &output_message: An output string containing messages from the code.
             */

            try {
                TARDIGRADE_ERROR_TOOLS_CHECK(getEvaluateModelResultSize() == result_size,
                                             "The output material size must be equal to " +
                                                 std::to_string(getEvaluateModelResultSize()))

                // Extract the parameters
                extract_parameters(parameters_begin, parameters_size);

                // Form the deformation gradients
                std::vector<current_dof_type>  deformationGradient;
                std::vector<previous_dof_type> previousDeformationGradient;
                formDeformationGradients(current_dof_begin, previous_dof_begin, deformationGradient,
                                         previousDeformationGradient);

                // Assemble the degrees of freedom
                tardigradeHydra::DOFStorageBase DOFStorage(
                    current_time, dt, *(current_dof_begin + getTemperatureIndex()),
                    *(previous_dof_begin + getTemperatureIndex()), deformationGradient, previousDeformationGradient,
                    std::vector<current_dof_type>(current_dof_begin, current_dof_begin + dof_size),
                    std::vector<previous_dof_type>(previous_dof_begin, previous_dof_begin + dof_size));

                // Assemble the model configuration
                tardigradeHydra::ModelConfigurationBase model_configuration(
                    std::vector<sdvs_type>(sdvs_begin, sdvs_begin + sdvs_size),
                    std::vector<parameter_type>(parameters_begin + 8, parameters_begin + parameters_size), 2, 9);

                // Form the hydra model
                definedPlasticEvolutionHydra hydra(DOFStorage, model_configuration);

                hydra.setDensityIndex(getDensityIndex());
                hydra.setDOFInternalEnergyIndex(getDOFInternalEnergyIndex());
                hydra.setDefinedVelocityGradientIndex(getDefinedVelocityGradientIndex());
                hydra.setInternalEnergyScaledByDensity(getInternalEnergyScaledByDensity());
                hydra.setDensityGradientIndex(getDensityGradientIndex());
                hydra.setTemperatureGradientIndex((unsigned int)(*(parameters_begin + 7) + 0.5));
                hydra.setStressParameters(parameters_begin + 8, parameters_begin + 10);
                hydra.setInternalEnergyParameters(parameters_begin + 10, parameters_begin + 11);
                hydra.setHeatConductionParameters(parameters_begin + 11, parameters_begin + 12);
                hydra.setDefinedDeformationParameters(parameters_begin + 12, parameters_begin + 14);
                hydra.setMassDiffusionParameters(parameters_begin + 14, parameters_begin + 15);

                hydra.evaluate();

                std::fill(result_begin, result_begin + result_size, 0);

                // Copy the stress response
                std::copy(std::begin(*hydra.getUnknownVector()), std::begin(*hydra.getUnknownVector()) + dim * dim,
                          result_begin + getStressIndex());

                // Copy the internal energy response
                *(result_begin + getInternalEnergyIndex()) = (*hydra.getUnknownVector())[dim * dim + dim * dim];
                *(sdvs_begin + dim * dim)                  = (*hydra.getUnknownVector())[dim * dim + dim * dim];

                // Copy the heat flux response
                std::copy(std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1,
                          std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim,
                          result_begin + getHeatFluxIndex());
                std::copy(std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1,
                          std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim,
                          sdvs_begin + dim * dim + 1);

                // Copy the mass diffusion response
                std::copy(std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim,
                          std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim + dim,
                          result_begin + getMassDiffusionIndex());
                std::copy(std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim,
                          std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim + dim,
                          sdvs_begin + dim * dim + 1 + dim);

                // Copy the mass change rate response
                *(result_begin + getMassChangeRateIndex()) =
                    (*hydra.getUnknownVector())[dim * dim + dim * dim + 1 + 2 * dim];
                *(sdvs_begin + dim * dim + 1 + 2 * dim) =
                    (*hydra.getUnknownVector())[dim * dim + dim * dim + 1 + 2 * dim];

                // Copy the internal heat generation rate response
                *(result_begin + getInternalHeatGenerationRateIndex()) =
                    (*hydra.getUnknownVector())[dim * dim + dim * dim + 1 + 2 * dim + 1];
                *(sdvs_begin + dim * dim + 1 + 2 * dim + 1) =
                    (*hydra.getUnknownVector())[dim * dim + dim * dim + 1 + 2 * dim + 1];

                // Copy the defined deformation response
                std::copy(std::begin(*hydra.getUnknownVector()) + dim * dim,
                          std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim, sdvs_begin);
                for (unsigned int i = 0; i < dim; ++i) {
                    *(sdvs_begin + 3 * i + i) -= 1;
                }

            } catch (tardigradeHydra::convergence_error &e) {
                return 1;

            } catch (std::exception &e) {
                tardigradeErrorTools::captureNestedExceptions(e, output_message);

                return 2;
            }

            return 0;
        }

        int DefinedPlasticEvolution::evaluate_model(const time_type &current_time, const time_type &dt,
                                                    const current_dof_type  *current_dof_begin,
                                                    const previous_dof_type *previous_dof_begin,
                                                    const unsigned int dof_size, const parameter_type *parameters_begin,
                                                    const unsigned int parameters_size, sdvs_type *sdvs_begin,
                                                    const unsigned int sdvs_size, result_type *result_begin,
                                                    const unsigned int result_size, jacobian_type *jacobian_begin,
                                                    additional_type   *additional_begin,
                                                    const unsigned int additional_size, std::string &output_message) {
            /*!
             * Evaluate the material model given the incoming data
             *
             * Returns:
             *     0: No errors. Solution converged.
             *     1: Convergence Error. Request timestep cutback.
             *     2: Fatal errors encountered. Terminate evaluation.
             *
             * \param &current_time: The current time of the simulation
             * \param &dt:           The change in time of the simulation
             * \param *current_dof_begin: A pointer to the starting element of the current values of the degrees of
             * freedom \param *previous_dof_begin: A pointer to the starting element of the previous values of the
             * degrees of freedom \param dof_size: The size of the dof arrays \param *parameters_begin: A pointer to the
             * starting element of the model parameters \param parameters_size: The size of the parameter vector \param
             * *sdvs_begin: A pointer to the starting element of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param sdvs_size: The size of the SDVS vector
             * \param *result_begin: A pointer to the starting element of the result vector.
             * \param result_size: The size of the result
             * \param *jacobian_begin: A pointer to the starting element of the jacobian of the result w.r.t. the dof.
             * (must have a size of result_size * dof_size) \param *additional_begin: A pointer to the starting element
             * of the additional information returned by the model \param additional_size: The size of the additional
             * information array \param &output_message: An output string containing messages from the code.
             */

            try {
                TARDIGRADE_ERROR_TOOLS_CHECK(getEvaluateModelResultSize() == result_size,
                                             "The output material size must be equal to " +
                                                 std::to_string(getEvaluateModelResultSize()))

                // Extract the parameters
                extract_parameters(parameters_begin, parameters_size);

                // Form the deformation gradients
                std::vector<current_dof_type>  deformationGradient;
                std::vector<previous_dof_type> previousDeformationGradient;
                std::vector<current_dof_type>  dFdGradU;
                formDeformationGradients(current_dof_begin, previous_dof_begin, deformationGradient,
                                         previousDeformationGradient, dFdGradU);

                // Assemble the degrees of freedom
                tardigradeHydra::DOFStorageBase DOFStorage(
                    current_time, dt, *(current_dof_begin + getTemperatureIndex()),
                    *(previous_dof_begin + getTemperatureIndex()), deformationGradient, previousDeformationGradient,
                    std::vector<current_dof_type>(current_dof_begin, current_dof_begin + dof_size),
                    std::vector<previous_dof_type>(previous_dof_begin, previous_dof_begin + dof_size));

                // Assemble the model configuration
                tardigradeHydra::ModelConfigurationBase model_configuration(
                    std::vector<sdvs_type>(sdvs_begin, sdvs_begin + sdvs_size),
                    std::vector<parameter_type>(parameters_begin + 8, parameters_begin + parameters_size), 2, 9);

                // Form the hydra model
                definedPlasticEvolutionHydra hydra(DOFStorage, model_configuration);

                hydra.setDensityIndex(getDensityIndex());
                hydra.setDOFInternalEnergyIndex(getDOFInternalEnergyIndex());
                hydra.setDefinedVelocityGradientIndex(getDefinedVelocityGradientIndex());
                hydra.setInternalEnergyScaledByDensity(getInternalEnergyScaledByDensity());
                hydra.setDensityGradientIndex(getDensityGradientIndex());
                hydra.setTemperatureGradientIndex((unsigned int)(*(parameters_begin + 7) + 0.5));
                hydra.setStressParameters(parameters_begin + 8, parameters_begin + 10);
                hydra.setInternalEnergyParameters(parameters_begin + 10, parameters_begin + 11);
                hydra.setHeatConductionParameters(parameters_begin + 11, parameters_begin + 12);
                hydra.setDefinedDeformationParameters(parameters_begin + 12, parameters_begin + 14);
                hydra.setMassDiffusionParameters(parameters_begin + 14, parameters_begin + 15);

                hydra.evaluate();

                std::fill(result_begin, result_begin + result_size, 0);

                // Copy the stress response
                std::copy(std::begin(*hydra.getUnknownVector()), std::begin(*hydra.getUnknownVector()) + dim * dim,
                          result_begin + getStressIndex());

                // Copy the internal energy response
                *(result_begin + getInternalEnergyIndex()) = (*hydra.getUnknownVector())[dim * dim + dim * dim];
                *(sdvs_begin + dim * dim)                  = (*hydra.getUnknownVector())[dim * dim + dim * dim];

                // Copy the heat flux response
                std::copy(std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1,
                          std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim,
                          result_begin + getHeatFluxIndex());
                std::copy(std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1,
                          std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim,
                          sdvs_begin + dim * dim + 1);

                // Copy the mass diffusion response
                std::copy(std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim,
                          std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim + dim,
                          result_begin + getMassDiffusionIndex());
                std::copy(std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim,
                          std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim + 1 + dim + dim,
                          sdvs_begin + dim * dim + 1 + dim);

                // Copy the mass change rate response
                *(result_begin + getMassChangeRateIndex()) =
                    (*hydra.getUnknownVector())[dim * dim + dim * dim + 1 + 2 * dim];
                *(sdvs_begin + dim * dim + 1 + 2 * dim) =
                    (*hydra.getUnknownVector())[dim * dim + dim * dim + 1 + 2 * dim];

                // Copy the internal heat generation rate response
                *(result_begin + getInternalHeatGenerationRateIndex()) =
                    (*hydra.getUnknownVector())[dim * dim + dim * dim + 1 + 2 * dim + 1];
                *(sdvs_begin + dim * dim + 1 + 2 * dim + 1) =
                    (*hydra.getUnknownVector())[dim * dim + dim * dim + 1 + 2 * dim + 1];

                // Copy the defined deformation response
                std::copy(std::begin(*hydra.getUnknownVector()) + dim * dim,
                          std::begin(*hydra.getUnknownVector()) + dim * dim + dim * dim, sdvs_begin);
                for (unsigned int i = 0; i < dim; ++i) {
                    *(sdvs_begin + 3 * i + i) -= 1;
                }

                // Construct the Jacobian
                hydra.computeTangents();
                hydra.computedXdAdditionalDOF();

                std::fill(jacobian_begin, jacobian_begin + result_size * dof_size, 0);

                // Copy over the Jacobians of the stresses w.r.t. the additional DOF
                std::copy(std::begin(*hydra.getFlatdXdAdditionalDOF()),
                          std::begin(*hydra.getFlatdXdAdditionalDOF()) + dim * dim * dof_size, jacobian_begin);

                // Add the contributions from the deformation gradient
                for (unsigned int i = 0; i < dim * dim; ++i) {
                    for (unsigned int j = 0; j < dim * dim; ++j) {
                        for (unsigned int k = 0; k < dim * dim; ++k) {
                            *(jacobian_begin + dof_size * (i + getStressIndex()) + k +
                              getDisplacementGradientIndex()) +=
                                (*hydra.getFlatdXdF())[dim * dim * i + j] * dFdGradU[dim * dim * j + k];
                        }
                    }
                }

                // Add the contributions from the temperature
                for (unsigned int i = 0; i < dim * dim; ++i) {
                    *(jacobian_begin + dof_size * (i + getStressIndex()) + getTemperatureIndex()) +=
                        (*hydra.getFlatdXdT())[i];
                }

                // Copy over the Jacobians of the internal energy w.r.t. the additional DOF
                std::copy(std::begin(*hydra.getFlatdXdAdditionalDOF()) + (dim * dim + dim * dim) * dof_size,
                          std::begin(*hydra.getFlatdXdAdditionalDOF()) + (dim * dim + dim * dim + 1) * dof_size,
                          jacobian_begin + getInternalEnergyIndex() * dof_size);

                // Add the contributions from the deformation gradient
                for (unsigned int i = 0; i < 1; ++i) {
                    for (unsigned int j = 0; j < dim * dim; ++j) {
                        for (unsigned int k = 0; k < dim * dim; ++k) {
                            *(jacobian_begin + dof_size * (i + getInternalEnergyIndex()) + k +
                              getDisplacementGradientIndex()) +=
                                (*hydra.getFlatdXdF())[dim * dim * (i + dim * dim + dim * dim) + j] *
                                dFdGradU[dim * j + k];
                        }
                    }
                }

                // Add the contributions from the temperature
                for (unsigned int i = 0; i < 1; ++i) {
                    *(jacobian_begin + dof_size * (i + getInternalEnergyIndex()) + getTemperatureIndex()) +=
                        (*hydra.getFlatdXdT())[i + dim * dim + dim * dim];
                }

                // Copy over the Jacobians of the heat flux w.r.t. the additional DOF
                std::copy(std::begin(*hydra.getFlatdXdAdditionalDOF()) + (dim * dim + dim * dim + 1) * dof_size,
                          std::begin(*hydra.getFlatdXdAdditionalDOF()) + (dim * dim + dim * dim + 1 + dim) * dof_size,
                          jacobian_begin + getHeatFluxIndex() * dof_size);

                // Add the contributions from the deformation gradient
                for (unsigned int i = 0; i < dim; ++i) {
                    for (unsigned int j = 0; j < dim * dim; ++j) {
                        for (unsigned int k = 0; k < dim * dim; ++k) {
                            *(jacobian_begin + dof_size * (i + getHeatFluxIndex()) + k +
                              getDisplacementGradientIndex()) +=
                                (*hydra.getFlatdXdF())[dim * dim * (i + dim * dim + dim * dim + 1) + j] *
                                dFdGradU[dim * j + k];
                        }
                    }
                }

                // Add the contributions from the temperature
                for (unsigned int i = 0; i < 1; ++i) {
                    *(jacobian_begin + dof_size * (i + getHeatFluxIndex()) + getTemperatureIndex()) +=
                        (*hydra.getFlatdXdT())[i + dim * dim + dim * dim + 1];
                }

                // Copy over the Jacobians of the mass diffusion w.r.t. the additional DOF
                std::copy(std::begin(*hydra.getFlatdXdAdditionalDOF()) + (dim * dim + dim * dim + 1 + dim) * dof_size,
                          std::begin(*hydra.getFlatdXdAdditionalDOF()) +
                              (dim * dim + dim * dim + 1 + dim + dim) * dof_size,
                          jacobian_begin + getMassDiffusionIndex() * dof_size);

                // Add the contributions from the deformation gradient
                for (unsigned int i = 0; i < dim; ++i) {
                    for (unsigned int j = 0; j < dim * dim; ++j) {
                        for (unsigned int k = 0; k < dim * dim; ++k) {
                            *(jacobian_begin + dof_size * (i + getMassDiffusionIndex()) + k +
                              getDisplacementGradientIndex()) +=
                                (*hydra.getFlatdXdF())[dim * dim * (i + dim * dim + dim * dim + 1 + dim) + j] *
                                dFdGradU[dim * j + k];
                        }
                    }
                }

                // Add the contributions from the temperature
                for (unsigned int i = 0; i < 1; ++i) {
                    *(jacobian_begin + dof_size * (i + getMassDiffusionIndex()) + getTemperatureIndex()) +=
                        (*hydra.getFlatdXdT())[i + dim * dim + dim * dim + 1 + dim];
                }

                // Copy over the Jacobians of the mass change rate w.r.t. the additional DOF
                std::copy(std::begin(*hydra.getFlatdXdAdditionalDOF()) +
                              (dim * dim + dim * dim + 1 + 2 * dim) * dof_size,
                          std::begin(*hydra.getFlatdXdAdditionalDOF()) +
                              (dim * dim + dim * dim + 1 + 2 * dim + 1) * dof_size,
                          jacobian_begin + getMassChangeRateIndex() * dof_size);

                // Copy over the Jacobians of the internal heat generation rate w.r.t. the additional DOF
                std::copy(std::begin(*hydra.getFlatdXdAdditionalDOF()) +
                              (dim * dim + dim * dim + 1 + 2 * dim + 1) * dof_size,
                          std::begin(*hydra.getFlatdXdAdditionalDOF()) +
                              (dim * dim + dim * dim + 1 + 2 * dim + 2) * dof_size,
                          jacobian_begin + getInternalHeatGenerationRateIndex() * dof_size);

            } catch (tardigradeHydra::convergence_error &e) {
                return 1;

            } catch (std::exception &e) {
                tardigradeErrorTools::captureNestedExceptions(e, output_message);

                return 2;
            }

            return 0;
        }

    }  // namespace DefinedPlasticEvolution

}  // namespace tardigradeCMML
