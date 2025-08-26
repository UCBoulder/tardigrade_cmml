/*!
=====================================================================
|                tardigrade_DefinedPlasticEvolution.h               |
=====================================================================
| A material where the plastic deformation is defined by an         |
| external velocity gradient.                                       |
---------------------------------------------------------------------
| The output material vector is comprised of the following          |
| responses.                                                        |
|                                                                   |
| Stress response          : Linear, isotropic elasticity           |
|     permanent deformation from externally provided velocity       |
|     gradient                                                      |
| Predicted energy         : Linearly related to the temperature    |
|     via the specific heat.                                        |
| Mass generation rate     : Zero                                   |
| Body force               : Zero                                   |
| Interphasic force        : Zero                                   |
| Heat flux                : Linearly related to the temperature    |
| Internal heat generation : Zero                                   |
| Interphasic heat transfer: Zero                                   |
| trace mass change velocity gradient: Zero                         |
| mass diffusion : linear w.r.t. density gradient                   |
|                                                                   |
| The class can be used as the base class for more advanced         |
| behaviors                                                         |
=====================================================================
*/

#ifndef TARDIGRADE_DEFINEDPLASTICEVOLUTION_H
#define TARDIGRADE_DEFINEDPLASTICEVOLUTION_H

#include<vector>
#include "tardigrade_error_tools.h"
#define TARDIGRADE_CMML_SKIP_HEADER_ONLY_MATERIALS_INCLUDE
#include "tardigrade_cmml.h"
#undef TARDIGRADE_CMML_SKIP_HEADER_ONLY_MATERIALS_INCLUDE
#include "tardigrade_BasicSolid.h"
#include "tardigrade_hydraDOFVelocityGradientDeformation.h"

namespace tardigradeCMML{

    namespace DefinedPlasticEvolution{

        /*!
         * The basic solid material model hydra class
         */
        class definedPlasticEvolutionHydra : public tardigradeCMML::BasicSolid::basicSolidHydra{

            public:

                using tardigradeCMML::BasicSolid::basicSolidHydra::basicSolidHydra;
                using tardigradeCMML::BasicSolid::basicSolidHydra::setResidualClasses;

                template<
                    class parameter_iter
                >
                void setDefinedDeformationParameters( const parameter_iter &value_start, const parameter_iter &value_end ){
                    /*!
                     * Set the parameters for the defined deformation calculation
                     * 
                     * \param &value_start: The starting value of the parameter iterator
                     * \param &value_end: The stopping value of the parameter iterator
                     */

                    _defined_deformation_parameters = std::vector< double >( value_start, value_end );

                }

                template<
                    class parameter_iter
                >
                void setMassDiffusionParameters( const parameter_iter &value_start, const parameter_iter &value_end ){
                    /*!
                     * Set the parameters for the mass diffusion
                     *
                     * \param &value_start: The starting value of the parameter iterator
                     * \param &value_end: The stopping value of the parameter iterator
                     */

                    _mass_diffusion_parameters = std::vector< double >( value_start, value_end );

                }

                void setDefinedVelocityGradientIndex( const unsigned int &value ){
                    /*!
                     * Set the starting index in the additional degrees of freedom of the additional velocity gradient that defines the permanent deformation
                     *
                     * \param &value: The value of the velocity gradient index
                     */

                    _defined_velocity_gradient_index = value;
                    _defined_velocity_gradient_index_set = true;

                }

                void setDensityGradientIndex( const unsigned int &value ){
                    /*!
                     * Set the starting index in the additional degrees of freedom of the density gradient that defines the mass diffusion
                     *
                     * \param &value: The value of the density gradient index
                     */

                    _density_gradient_index = value;
                    _density_gradient_index_set = true;

                }

                const unsigned int getDefinedVelocityGradientIndex( ){
                    /*!
                     * Get the starting index of the velocity gradient that defines the permanent deformation in the additional degree of freedom vector
                     */

                    TARDIGRADE_ERROR_TOOLS_CHECK( _defined_velocity_gradient_index_set, "The defined velocity gradient index must be set before it is called" );

                    return _defined_velocity_gradient_index;

                }

                const unsigned int getDensityGradientIndex( ){
                    /*!
                     * Get the starting index of the density gradient that defines the mass diffusion
                     */

                    TARDIGRADE_ERROR_TOOLS_CHECK( _density_gradient_index_set, "The density gradient index must be set before it is called" );

                    return _density_gradient_index;

                }

                virtual void setResidualClasses( ){
                    /*!
                     * Define all of the residual classes
                     */

                    stress =
                        tardigradeHydra::linearElasticity::residual(
                            this, getStressSize( ), *getStressParameters( )
                        );

                    defined_deformation =
                        tardigradeHydra::dofVelocityGradientDeformation::residual(
                            this, getStressSize( ), 1, getDefinedVelocityGradientIndex( ), *getDefinedDeformationParameters( ), getIntegrationParameter( )
                        );

                    internal_energy =
                        tardigradeHydra::linearInternalEnergy::residual(
                            this, getInternalEnergySize( ), *getInternalEnergyParameters( ), 1,  18
                        );

                    heat_conduction =
                        tardigradeHydra::fourierHeatConduction::residual(
                            this, getHeatConductionSize( ), *getHeatConductionParameters( ), getTemperatureGradientIndex( ), 1, 19
                        );

                    mass_diffusion =
                        tardigradeHydra::fourierHeatConduction::residual(
                            this, getHeatConductionSize( ), *getMassDiffusionParameters( ), getDensityGradientIndex( ), 1, 22
                        );

                    std::vector< tardigradeHydra::residualBase* > residuals( 5 );

                    residuals[ 0 ] = &stress;

                    residuals[ 1 ] = &defined_deformation;

                    residuals[ 2 ] = &internal_energy;

                    residuals[ 3 ] = &heat_conduction;

                    residuals[ 4 ] = &mass_diffusion;

                    setResidualClasses( residuals );

                }

                const std::vector< double > *getDefinedDeformationParameters( ){ /*! Get the defined deformation parameters */ return &_defined_deformation_parameters; }

                const std::vector< double > *getMassDiffusionParameters( ){ /*! Get the mass diffusion parameters */ return &_mass_diffusion_parameters; }

                const double getIntegrationParameter( ){ /*! Get the integration parameter */ return _integration_parameter; }

            protected:

                tardigradeHydra::fourierHeatConduction::residual mass_diffusion; //!< The residual that defines the mass diffusion (mathematically the same as Fourier heat conduction with a density gradient)
                tardigradeHydra::dofVelocityGradientDeformation::residual defined_deformation; //!< The residual that defines the externally defined permanent deformation

            private:

                double                            _integration_parameter = 0.5; //!< The integration parameter for the defined deformation (0 is explicit, 1 is implicit)
                std::vector< double >          _defined_deformation_parameters; //!< The parameters associated with the calculation of the defined deformation
                unsigned int                  _defined_velocity_gradient_index; //!< The index in the additional DOF vector that the defined velocity gradient is defined in
                bool              _defined_velocity_gradient_index_set = false; //!< The flag for if _defined_velocity_gradient_index is set
                unsigned int                           _density_gradient_index; //!< The index in the additional DOF vector that the denensity gradient is defined in
                bool                       _density_gradient_index_set = false; //!< The flag for if _density_gradient_index is set
                std::vector< double >          _mass_diffusion_parameters; //!< The parameters associated with the mass diffusion

        };

        /*!
         * A material where a plastic configuration evolves due to an externally defined velocity gradient
         */
        class DefinedPlasticEvolution : public tardigradeCMML::BasicSolid::BasicSolid
        {

            public:

                DefinedPlasticEvolution( ){
                    /*!
                     * A material where a plastic configuration evolves 
                     * due to an externally defined velocity gradient
                     *
                     * Stress response          : Linear, isotropic elasticity
                     *     permanent deformation from externally provided velocity
                     *     gradient
                     * Predicted energy         : Linearly related to the temperature
                     *     via the specific heat.
                     * Mass generation rate     : Zero
                     * Body force               : Zero
                     * Interphasic force        : Zero
                     * Heat flux                : Linearly related to the temperature
                     * Internal heat generation : Zero
                     * Interphasic heat transfer: Zero
                     * trace mass change velocity gradient: Zero
                     * mass diffusion: Linearly related to density gradient
                     *
                     * The class can be used as the base class for more advanced
                     * behaviors
                     */

                    setName( "DefinedPlasticEvolution" );
                    setEvaluateModelResultSize( 26 );

                }

                virtual int evaluate_model(
                    const time_type         &current_time,       const time_type &dt,
                    const current_dof_type  *current_dof_begin,  const previous_dof_type *previous_dof_begin,
                    const unsigned int      dof_size,
                    const parameter_type    *parameters_begin,
                    const unsigned int      parameters_size,
                    sdvs_type               *sdvs_begin,
                    const unsigned int      sdvs_size,
                    result_type             *result_begin,
                    const unsigned int      result_size,
                    std::string             &output_message
                ) override;

                virtual int evaluate_model(
                    const time_type         &current_time,       const time_type &dt,
                    const current_dof_type  *current_dof_begin,  const previous_dof_type *previous_dof_begin,
                    const unsigned int      dof_size,
                    const parameter_type    *parameters_begin,
                    const unsigned int      parameters_size,
                    sdvs_type               *sdvs_begin,
                    const unsigned int      sdvs_size,
                    result_type             *result_begin,
                    const unsigned int      result_size,
                    jacobian_type           *jacobian_begin,
                    additional_type         *additional_begin,
                    const unsigned int      additional_size,
                    std::string             &output_message
                ) override;

                void setDefinedVelocityGradientIndex( unsigned int value ){
                    /*!
                     * Set the index the velocity gradient that defines the permanent deformation is located in the dof vector
                     * 
                     * \param value: The velocity gradient index
                     */

                    _defined_velocity_gradient_index = value;

                }

                void setDensityGradientIndex( unsigned int value ){
                    /*!
                     * Set the index the density gradient that defines the mass diffusion is located in the dof vector
                     * 
                     * \param value: The density gradient index
                     */

                    _density_gradient_index = value;

                }

                void setMassDiffusionCoefficient( double value ){
                    /*!
                     * Set the mass diffusion coefficient
                     */

                    _mass_diffusion_coefficient = value;

                }

                const unsigned int getDefinedVelocityGradientIndex( ){ /*! Get the defined velocity gradient index */ return _defined_velocity_gradient_index; }

                const unsigned int getDensityGradientIndex( ){ /*! Get the density gradient index */ return _density_gradient_index; }

                const double getMassDiffusionCoefficient( ){ /*! Get the mass diffusion coefficient */ return _mass_diffusion_coefficient; }

                const unsigned int getMassDiffusionIndex( ){
                    /*!
                     * Get the mass diffusion index in the result vector
                     */

                    return _mass_diffusion_index;
                }

            protected:

                virtual void extract_parameters(const  parameter_type *parameters_begin, const unsigned int parameters_size ){
                    /*!
                     * Extract the material parameters
                     *
                     * The parameter vector is assumed to be organized as
                     *
                     * defined velocity gradient index, density gradient index, displacement gradient index, temperature index, temperature gradient index, lambda, mu, specifc heat, conductivity,
                     * mass diffusion coefficient
                     *
                     * where lambda and mu are the linear elastic Lame parameters, specific heat is the linear
                     * relationship between temperature and the internal energy per unit mass, and the conductivity is the
                     * conductivity parameter for the heat flux.
                     *
                     * \param *parameters_begin: A pointer to the initial element of the parameter vector
                     * \param parameters_size: The size of the parameters vector
                     */

                    TARDIGRADE_ERROR_TOOLS_CHECK( parameters_size == 10, "The parameters vector must have ten values" )

                    setDefinedVelocityGradientIndex( ( unsigned int )( *( parameters_begin + 0 ) + 0.5 ) );
                    setDensityGradientIndex( ( unsigned int )( *( parameters_begin + 1 ) + 0.5 ) );
                    setMassDiffusionCoefficient( *( parameters_begin + parameters_size - 1 ) );

                    BasicSolid::extract_parameters( parameters_begin + 2, parameters_size - 3 );

                }

            private:

                unsigned int _defined_velocity_gradient_index; //!< The index of the velocity gradient that defines the deformation evolution

                unsigned int _density_gradient_index; //!< The index of the density gradient that defines the mass diffusion

                float _mass_diffusion_coefficient; //!< The diffusion coefficient of the mass

                unsigned int _mass_diffusion_index = 23; //!< The index of the mass diffusion in the result vector
        };

        //! Register the material in the library
        REGISTER_MATERIAL( tardigradeCMML::DefinedPlasticEvolution::DefinedPlasticEvolution );

    }

}

#include "tardigrade_DefinedPlasticEvolution.cpp"

#endif
