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

namespace tardigradeCMML {

    namespace multiphase {

        /*! The base class for multiphasic continuum models */
        class CMMLMultiphaseMaterialBase : public tardigradeCMML::CMMLMaterial {
           public:
            CMMLMultiphaseMaterialBase() {
                /*!
                 * The default constructor for a multiphasic material model
                 */
            }
        };

    }  // namespace multiphase

}  // namespace tardigradeCMML

#include "tardigrade_CMMLMultiphaseMaterialBase.cpp"

#endif
