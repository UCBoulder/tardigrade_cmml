/*!
=====================================================================
|                         tardigrade_cmml.h                         |
=====================================================================
| (C)ollected (M)aterial (M)odel (L)ibrary                          |
| prounced Camel                                                    |
|                                                                   |
| A header file which defines a class which registers all of the    |
| available material models and their interface. The library is     |
| implemented in such a way that someone can register a new         |
| constitutive model and have it available for use merely by        |
| calling it by name. This allows us to re-use code as much as      |
| possible without needing to do time consuming rewrites of already |
| existing code.                                                    |
---------------------------------------------------------------------
| Note: Registration approach taken from stackoverflow question     |
|       compile time plugin system 2                                |
=====================================================================
*/

#ifndef TARDIGRADE_CMML_H
#define TARDIGRADE_CMML_H

#include<array>
#include<math.h>
#include<list>
#include<string>
#include<map>
#include<memory>
#include<iostream>

namespace tardigradeCMML{

    typedef double time_type; //!< Define the time variable type
    typedef double current_dof_type; //!< Define the type of the current degree of freedom
    typedef double previous_dof_type; //!< Define the type of the previous degree of freedom
    typedef double parameter_type; //!< Define the type of the parameter vector
    typedef double sdvs_type; //!< Define the type of the sdvs vector
    typedef double result_type; //!< The type of the result vector
    typedef double jacobian_type; //!< The type of the jacobian vector
    typedef double additional_type; //!< The type of the additional information returned by the model

    /*!
     * The base material model class. Used primarily so that the self registration macro works.
     */
    class CMMLMaterialBase{

        public:

            CMMLMaterialBase( ){
                /*!
                 * The constructor for the base class
                 */

            }

            const std::string* getName( ){
                /*!
                 * Return a reference to the material name
                 */

                return &_name;

            }

            const unsigned int getEvaluateModelResultSize( ){
                /*!
                 * Get the expected size of the output from evaluate model
                 */

                return _evaluate_model_result_size;

            }

        protected:

            void setName( std::string name ){
                /*!
                 * Set the material model's name
                 * 
                 * \param name: The new name for the material model
                 */

                 _name = name;

            }

            void setEvaluateModelResultSize( const unsigned int &result ){
                /*!
                 * Set the size of the result from evaluate model
                 * 
                 * \param &result: The new value
                 */

                _evaluate_model_result_size = result;

            }

        private:
            std::string _name;
            unsigned int _evaluate_model_result_size;

    };

    /*!
     * Class template which defines the general interface to constitutive models
     */
    class CMMLMaterial : public CMMLMaterialBase{

        public:
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
            ){
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
             * \param *current_dof_begin: A pointer to the starting element of the current values of the degrees of freedom
             * \param *previous_dof_begin: A pointer to the starting element of the previous values of the degrees of freedom
             * \param dof_size: The size of the dof arrays
             * \param *parameters_begin: A pointer to the starting element of the model parameters
             * \param parameters_size: The size of the parameter vector
             * \param *sdvs_begin: A pointer to the starting element of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param sdvs_size: The size of the SDVS vector
             * \param *result_begin: A pointer to the starting element of the result vector.
             * \param result_size: The size of the result
             * \param &output_message: An output string containing messages from the code.
             */

                output_message = "evaluate_model is not implemented";
                return 2;

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
                jacobian_type           *jacobian_begin,
                additional_type         *additional_begin,
                const unsigned int      additional_size,
                std::string             &output_message
            ){
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
             * \param *current_dof_begin: A pointer to the starting element of the current values of the degrees of freedom
             * \param *previous_dof_begin: A pointer to the starting element of the previous values of the degrees of freedom
             * \param dof_size: The size of the dof arrays
             * \param *parameters_begin: A pointer to the starting element of the model parameters
             * \param parameters_size: The size of the parameter vector
             * \param *sdvs_begin: A pointer to the starting element of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param sdvs_size: The size of the SDVS vector
             * \param *result_begin: A pointer to the starting element of the result vector.
             * \param result_size: The size of the result
             * \param *jacobian_begin: A pointer to the starting element of the jacobian of the result w.r.t. the dof. (must have a size of result_size * dof_size)
             * \param *additional_begin: A pointer to the starting element of the additional information returned by the model
             * \param additional_size: The size of the additional information array
             * \param &output_message: An output string containing messages from the code.
             */

                output_message = "evaluate_model with additional information is not implemented";
                return 2;

            }

    };

    /*!
     * The base class for the CMML material registrar
     */
    class CMMLMaterialRegistrar{

        public:

            //! Get a pointer to the requested material
            virtual std::unique_ptr< CMMLMaterial > GetMaterial( ) = 0;

    };

    /*!
     * The factory that provides the common interface to CMML materials.
     * Materials register themselves here and the factory can serve them on
     * demand.
     * It is a Singleton
     */

    class MaterialFactory {

        public:
            static MaterialFactory& Instance( ){
                /*! Return an instance of the factory */
                static MaterialFactory instance;
                return instance;
            }

            void Register( CMMLMaterialRegistrar* registrar, std::string name ){
                /*!
                 * Register the material model in the library
                 * 
                 * \param *registrar: The registrar object
                 * \param name: The name of the model
                 */

                 registry_[name] = registrar;
            }

            std::unique_ptr< CMMLMaterial > GetMaterial( std::string name ){
                /*!
                 * Get a reference to the material model.
                 * 
                 * Throws out_of_range if material is unknown
                 * 
                 * \param name: The name of the model to retrieve
                 */

                try{
                    CMMLMaterialRegistrar* registrar = registry_.at(name);
                    return registrar->GetMaterial( );
                }
                catch(...){
                    std::string message = "Exception when attempting to access: " + name + "\n";
                    message += "material models available are:\n";
                    for ( auto it = registry_.begin( ); it != registry_.end( ); it++ ){
                        message += it->first + "\n";
                    }
                    throw std::runtime_error( message );
                }
                return NULL;
            }

            void PrintMaterials( ){
                /*! Prints all of the materials registered in the library*/
                std::string message = "Materials available in the library:\n";
                for ( auto it = registry_.begin( ); it != registry_.end( ); it++ ){
                    message += "    " + it->first + "\n";
                }
                std::cerr << message;
            }

        private:
            /*! Holds pointers to material registrars */
            std::map< std::string, CMMLMaterialRegistrar* > registry_;

            /*! Make constructors private and forbid cloning */
            MaterialFactory( ): registry_( ) { };
            
            MaterialFactory( MaterialFactory const & ) = delete;

            void operator = ( MaterialFactory const & ) = delete;

    };

    /*!
     * Helper class that registers a material upon construction.
     * Actually, the registrar registers itself, and the proxied material is only
     * created on-demand. This mechanism can be shortened by directly
     * registering an instance of the material, but the assumption here is that
     * instanciating the material can be heavy and not necessary.
     */

    template<
        class TCMMLMaterial
    >
    class MaterialRegistrar : public CMMLMaterialRegistrar{

        public:
            MaterialRegistrar( std::string classname );

            std::unique_ptr< CMMLMaterial > GetMaterial( );

        private:
            /*! That is not really used here, but could be useful */
            std::string classname_;

    };

    /* template functions in header */

    template<
        class TCMMLMaterial
    >
    MaterialRegistrar< TCMMLMaterial >::MaterialRegistrar( std::string classname ): classname_( classname ){
        /*!
         * Build a registration to the provided class
         */
        MaterialFactory &factory = MaterialFactory::Instance( );
        factory.Register( this, classname );
    }

    template<
        class TCMMLMaterial
    >
    std::unique_ptr< CMMLMaterial >
    MaterialRegistrar<TCMMLMaterial>::GetMaterial( ){
        /*!
         * Get a reference to the requested material model
         */
        std::unique_ptr< CMMLMaterial > material( new TCMMLMaterial( ) );
        return material;
    }

}

/*
 * Here is the trick: upon creation of the global variable, the class created
 * out of the template will get instanciated once, and will register itself.
 * The template contains the information to create a material instance.
 * An unnamed namespace is used to enclose this later unused variable in the
 * compilation unit.
 */
#define REGISTER_MATERIAL(CLASSNAME) \
    namespace { \
        static tardigradeCMML::MaterialRegistrar<CLASSNAME> \
        _registrar( #CLASSNAME ); \
    }

#include "tardigrade_cmml.cpp"

#endif
