/*!
=====================================================================
|                 tardigrade_BasicReactingSolid.h                   |
=====================================================================
| A basic reacting solid.                                           |
---------------------------------------------------------------------
| The output material vector is comprised of the following          |
| responses.                                                        |
|                                                                   |
| Stress response          : Linear, isotropic elasticity           |
| Predicted energy         : Linearly related to the temperature    |
|     via the specific heat.                                        |
| Mass generation rate     : Zero                                   |
| Body force               : Zero                                   |
| Interphasic force        : Zero                                   |
| Heat flux                : Linearly related to the temperature    |
| Internal heat generation : Zero                                   |
| Interphasic heat transfer: Zero                                   |
| trace mass change velocity gradient: Zero                         |
|                                                                   |
| The class can be used as the base class for more advanced         |
| behaviors                                                         |
=====================================================================
*/

#ifndef TARDIGRADE_GENERALREACTINGSOLID_H
#define TARDIGRADE_GENERALREACTINGSOLID_H

#include<vector>
#include "tardigrade_error_tools.h"
#include "tardigrade_cmml.h"
#include "tardigrade_hydraLinearElasticity.h"
#include "tardigrade_hydraLinearInternalEnergy.h"
#include "tardigrade_hydraFourierHeatConduction.h"

namespace tardigradeCMML{

    class hydra : public tardigradeHydra::hydraBase{

        public:

            using tardigradeHydra::hydraBase::hydraBase;
            using tardigradeHydra::hydraBase::setResidualClasses;

            template<
                class parameter_iter
            >
            void setStressParameters( const parameter_iter &value_start, const parameter_iter &value_end ){
                /*!
                 * Set the parameters for the stress calculation
                 *
                 * \param &value_start: The starting value of the parameter iterator
                 * \param &value_end: The stopping value of the parameter itertor
                 */

                _stress_parameters = std::vector< double >( value_start, value_end );

            }

            template<
                class parameter_iter
            >
            void setInternalEnergyParameters( const parameter_iter &value_start, const parameter_iter &value_end ){
                /*!
                 * Set the parameters for the internal energy
                 *
                 * \param &value_start: The starting value of the parameter iterator
                 * \param &value_end: The stopping value of the parameter itertor
                 */

                _internal_energy_parameters = std::vector< double >( value_start, value_end );

            }

            template<
                class parameter_iter
            >
            void setHeatConductionParameters( const parameter_iter &value_start, const parameter_iter &value_end ){
                /*!
                 * Set the parameters for the heat conduction
                 *
                 * \param &value_start: The starting value of the parameter iterator
                 * \param &value_end: The stopping value of the parameter itertor
                 */

                _heat_conduction_parameters = std::vector< double >( value_start, value_end );

            }

            virtual void setResidualClasses( ){

                stress          = tardigradeHydra::linearElasticity::residual( this, getStressSize( ), *getStressParameters( ) );

                internal_energy = tardigradeHydra::linearInternalEnergy::residual( this, getInternalEnergySize( ), *getInternalEnergyParameters( ), 1,  9 );

                heat_conduction = tardigradeHydra::linearInternalEnergy::residual( this, getHeatConductionSize( ), *getHeatConductionParameters( ), *getTemperatureGradientIndex( ), 1, 10 );

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                residuals[ 0 ] = &stress;

                residuals[ 1 ] = &internal_energy;

                residuals[ 3 ] = &heat_conduction;

                setResidualClasses( residuals );

            }

            const std::vector< double > *getStressParameters( ){ /*! Get the stress parameters */ return &_stress_parameters; }

            const std::vector< double > *getInternalEnergyParameters( ){ /*! Get the internal energy parameters */ return &_internal_energy_parameters; }

            const std::vector< double > *getHeatConductionParameters( ){ /*! Get the heat conduction parameters */ return &_heat_conduction_parameters; }

            const unsigned int getStressSize( ){ /*! Get the number of nonlinear unknowns in the stress response */ return _stress_size; }

            const unsigned int getInternalEnergySize( ){ /*! Get the number of nonlinear unknowns in the internal energy response */ return _internal_energy_size; }

            const unsigned int getHeatConductionSize( ){ /*! Get the number of nonlinear unknowns in the heat conduction response */ return _heat_conduction_size; }

        protected:

            tardigradeHydra::linearElasticity::residual      stress;
            tardigradeHydra::linearInternalEnergy::residual  internal_energy;
            tardigradeHydra::fourierHeatConduction::residual heat_conduction;

        private:

            unsigned int _stress_size = 9; //!< The number of nonlinear unknowns in the stress response (must be at least 9)
            unsigned int _internal_energy_size = 1; //!< The number of nonlinear unknowns in the internal energy response (must be at least 1)
            unsigned int _heat_conduction_size = 3; //!< The number of nonlinear unknowns in the heat conduction response (must be at least 3)
            std::vector< double > _stress_parameters; //!< The parameters associated with the calculation of the stress
            std::vector< double > _internal_energy_parameters; //!< The parameters associated with the calculation of the internal energy
            std::vector< double > _heat_conduction_parameters; //!< The parameters associated with the calculation of the heat conduction

    };

    /*!
     * A general reacting continuum solid
     */
    template<
        typename time_type,
        class current_dof_iter,
        class previous_dof_iter,
        class parameter_iter,
        class sdvs_iter,
        class result_iter,
        class jacobian_iter,
        class additional_iter
    >
    class BasicReactingSolid : public tardigradeCMML::CMMLMaterial<
        time_type, current_dof_iter, previous_dof_iter, parameter_iter,
        sdvs_iter, result_iter, jacobian_iter, additional_iter
    >{

        public:

            using tardigradeCMML::CMMLMaterialBase::setName;

            BasicReactingSolid( ){
                /*!
                 * A general reacting solid.
                 * 
                 * Stress response          : Linear, isotropic elasticity
                 * Predicted energy         : Linearly related to the temperature
                 *     via the specific heat.
                 * Mass generation rate     : Zero
                 * Body force               : Zero
                 * Interphasic force        : Zero
                 * Heat flux                : Linearly related to the temperature
                 * Internal heat generation : Zero
                 * Interphasic heat transfer: Zero
                 * trace mass change velocity gradient: Zero
                 *
                 * The class can be used as the base class for more advanced
                 * behaviors
                 * 
                 */

                setName( "BasicReactingSolid" );

            }

            void setDisplacementIndex( unsigned int value ){
                /*!
                 * Set the index the displacement is located in the dof vector
                 * 
                 * \param value: The displacement index
                 */

                _displacement_index = value;

            }

            void setTemperatureIndex( unsigned int value ){
                /*!
                 * Set the index the temperature is located in the dof vector
                 * 
                 * \param value: The displacement index
                 */

                _temperature_index = value;

            }

            virtual int evaluate_model(
                const time_type         &current_time,       const time_type &dt,
                const current_dof_iter  &current_dof_begin,  const current_dof_iter  &current_dof_end,
                const previous_dof_iter &previous_dof_begin, const previous_dof_iter &previous_dof_end,
                const parameter_iter    &parameters_begin,   const parameter_iter    &parameters_end,
                sdvs_iter               sdvs_begin,          sdvs_iter               sdvs_end,
                result_iter             result_begin,        result_iter             result_end,
                std::string             &output_message
            );

            virtual int evaluate_model(
                const time_type         &current_time,       const time_type &dt,
                const current_dof_iter  &current_dof_begin,  const current_dof_iter  &current_dof_end,
                const previous_dof_iter &previous_dof_begin, const previous_dof_iter &previous_dof_end,
                const parameter_iter    &parameters_begin,   const parameter_iter    &parameters_end,
                sdvs_iter               sdvs_begin,          sdvs_iter               sdvs_end,
                result_iter             result_begin,        result_iter             result_end,
                jacobian_iter           jacobian_begin,      jacobian_iter           jacobian_end,
                additional_iter         additional_begin,    additional_iter         additional_end,
                std::string             &output_message
            );

        protected:

            virtual void extract_parameters(const  parameter_iter &parameters_begin, const parameter_iter &parameters_end ){
                /*!
                 * Extract the material parameters
                 *
                 * The parameter vector is assumed to be organized as
                 *
                 * displacement index, temperature index, lambda, mu, specifc heat, conductivity
                 *
                 * where lambda and mu are the linear elastic Lame parameters, specific heat is the linear
                 * relationship between temperature and the internal energy, and the conductivity is the
                 * conductivity parameter for the heat flux.
                 *
                 * \param &parameters_begin: The starting iterator of the parameter vector
                 * \param &parameters_end: The stopping iterator of the parameter vector
                 */

                TARDIGRADE_ERROR_TOOLS_CHECK( ( parameters_end - parameters_begin ) == 6, "The parameters vector must have six indices" )

                setDisplacementIndex( ( unsigned int )( *( parameters_begin + 0 ) + 0.5 ) );
                setTemperatureIndex(  ( unsigned int )( *( parameters_begin + 1 ) + 0.5 ) );

            }

        private:

            unsigned int _expected_material_size = 23;
            unsigned int _displacement_index;
            unsigned int _temperature_index;

    };

//    REGISTER_MATERIAL( tardigradeCMML::BasicReactingSolid );

}

#include "tardigrade_BasicReactingSolid.cpp"

#endif
