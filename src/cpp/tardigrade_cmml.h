/*!
=====================================================================
|                         tardigrade_cmml.h                         |
=====================================================================
| (C)ollected (M)aterial (M)odel (L)ibrary                          |
| prounced Cammel                                                   |
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

        protected:

            void setName( std::string name ){
                /*!
                 * Set the material model's name
                 * 
                 * \param name: The new name for the material model
                 */

                 _name = name;

            }

        private:
            std::string _name;

    };

    /*!
     * Class template which defines the general interface to constitutive models
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
    class CMMLMaterial : public CMMLMaterialBase{

        public:
            virtual int evaluate_model(
                const time_type         &current_time,       const time_type &dt,
                const current_dof_iter  &current_dof_begin,  const current_dof_iter  &current_dof_end,
                const previous_dof_iter &previous_dof_begin, const previous_dof_iter &previous_dof_end,
                const parameter_iter    &parameters_begin,   const parameter_iter    &parameters_end,
                sdvs_iter               sdvs_begin,          sdvs_iter               sdvs_end,
                result_iter             result_begin,        result_iter             result_end,
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
             * \param &current_dof_begin: The starting iterator of the current values of the degrees of freedom
             * \param &current_dof_end: The stopping iterator of the current values of the degrees of freedom
             * \param &previous_dof_begin: The starting iterator of the previous values of the degrees of freedom
             * \param &previous_dof_end: The stopping iterator of the previous values of the degrees of freedom
             * \param &parameters_begin: The starting iterator of the model parameters
             * \param &parameters_end: The stopping iterator of the model parameters
             * \param &sdvs_begin: The starting iterator of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param &sdvs_end: The stopping iterator of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param &result_begin: The starting iterator of the result vector.
             * \param &result_end: The stopping iterator of the result vector.
             * \param &output_message: An output string containing messages from the code.
             */

                output_message = "evaluate_model is not implemented";
                return 2;

            }

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
             * \param &current_dof_begin: The starting iterator of the current values of the degrees of freedom
             * \param &current_dof_end: The stopping iterator of the current values of the degrees of freedom
             * \param &previous_dof_begin: The starting iterator of the previous values of the degrees of freedom
             * \param &previous_dof_end: The stopping iterator of the previous values of the degrees of freedom
             * \param &parameters_begin: The starting iterator of the model parameters
             * \param &parameters_end: The stopping iterator of the model parameters
             * \param &sdvs_begin: The starting iterator of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param &sdvs_end: The stopping iterator of the model state variables. These are initialized to the
             *     previous converged values and should be updated by the model
             * \param &result_begin: The starting iterator of the result vector.
             * \param &result_end: The stopping iterator of the result vector.
             * \param &jacobian_begin: The starting iterator of the jacobian of the result w.r.t. the dof.
             * \param &jacobian_end: The stopping iterator of the jacobian of the result w.r.t. the dof.
             * \param &additional_begin: The starting iterator of the additional information returned by the model
             * \param &additional_end: The stopping iterator of the additional information returned by the model
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
            virtual std::unique_ptr< CMMLMaterialBase > GetMaterial( ) = 0;

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

            std::unique_ptr< CMMLMaterialBase > GetMaterial( std::string name ){
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

            std::unique_ptr< CMMLMaterialBase > GetMaterial( );

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
    std::unique_ptr< CMMLMaterialBase >
    MaterialRegistrar<TCMMLMaterial>::GetMaterial( ){
        /*!
         * Get a reference to the requested material model
         */
        std::unique_ptr< CMMLMaterialBase > material( new TCMMLMaterial( ) );
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

#endif
