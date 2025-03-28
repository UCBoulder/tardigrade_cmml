/**
  *****************************************************************************
  * \file tardigrade_CMMLMultiphaseMaterialBase.h
  *****************************************************************************
  * A definition of the base class for multiphasic material responses
  *****************************************************************************
  */

#ifndef TARDIGRADE_MULTIPHASE_MATERIAL_RESPONSE_H
#define TARDIGRADE_MULTIPHASE_MATERIAL_RESPONSE_H

#include "tardigrade_cmml.h"

namespace tardigradeCMML{

    namespace multiphase{

        /*! The base class for multiphasic continuum models */
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
        class CMMLMultiphaseMaterialBase : public tardigradeCMML::CMMLMaterial<
            time_type, current_dof_iter, previous_dof_iter, parameter_iter,
            sdvs_iter, result_iter, jacobian_iter, additional_iter
        >{

            public:

                CMMLMultiphaseMaterialBase( ){
                    /*!
                     * The default constructor for a multiphasic material model
                     */
                }

        };

    }

}

#include "tardigrade_CMMLMultiphaseMaterialBase.cpp"

#endif
