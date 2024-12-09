/**
  ******************************************************************************
  * \file tardigrade_finite_element_utilities.h
  ******************************************************************************
  * The header file for utilities which can assist using the balance equations
  * in finite element codes. We here assume that unknown quantities are
  * computed at the evaluation point using interpolation functions (here called
  * interp) and are projected to the nodes using test functions (here called
  * test).
  ******************************************************************************
  */

#ifndef TARDIGRADE_FINITE_ELEMENT_UTILITIES_H
#define TARDIGRADE_FINITE_ELEMENT_UTILITIES_H

#include<array>

#include "tardigrade_error_tools.h"

namespace tardigradeBalanceEquations{

    namespace finiteElementUtilities{

        typedef unsigned int size_type; //!< Define unsigned int as the default size type

        constexpr unsigned int dim = 3; //!< Set the dimension as 3D by default

        constexpr unsigned int sot_dim = dim * dim; //!< Set the dimensions of a standard second order tensor

        typedef double floatType; //!< Define the float type as a double

        typedef std::array< floatType, dim > floatVector; //!< Define a standard vector

        typedef std::array< floatType, sot_dim > secondOrderTensor; //!< Define a standard second-order tensor

        template<typename grad_iterator, typename output_iterator>
        void computeGradientSpatialJacobian( const grad_iterator &grad_a_start, const unsigned int grad_a_size, 
                                             floatVector grad_interp, const unsigned int index, output_iterator dgrad_adui_start );

        //! A base class for a simple finite element formulation useful for testing
        template<int dim, class node_in, class local_node_in, class local_point_in, class shape_functions_out, class grad_shape_functions_out>
        class FiniteElementBase{

            public:

                FiniteElementBase( const node_in &_x_begin, const node_in &_x_end, const node_in &_X_begin, const node_in &_X_end,
                                   const local_node_in &_local_node_xi_begin, const local_node_in &_local_node_xi_end );

                virtual void GetShapeFunctions( const local_point_in &xi_begin, const local_point_in &xi_end, shape_functions_out N_begin, shape_functions_out N_end ){
                    /*!
                     * Get the values of the shape functions
                     * 
                     * \param &xi_begin: The starting iterator of the local point
                     * \param &xi_end: The stopping iterator of the local point
                     * \param N_begin: The starting iterator of the shape functions
                     * \param N_end: The stopping iterator of the shape functions
                     */

                    throw std::logic_error("Function not implemented");

                }

                virtual void GetLocalShapeFunctionGradients( const local_point_in &xi_begin, const local_point_in &xi_end, grad_shape_functions_out dNdxi_begin, grad_shape_functions_out dNdxi_end ){
                    /*!
                     * Get the values of the gradients of the shape functions w.r.t. the local coordinates
                     * 
                     * \param &xi_begin: The starting iterator of the local point
                     * \param &xi_end: The stopping iterator of the local point
                     * \param dNdxi_begin: The starting iterator of the shape function gradients
                     * \param dNdxi_end: The stopping iterator of the shape function gradients
                     */

                    throw std::logic_error("Function not implemented");

                }

            protected:

                const node_in x_begin; //! Starting iterator for the current position of the nodal coordinates
                const node_in x_end; //! Stopping iterator for the current position of the nodal coordinates

                const node_in X_begin; //! Starting iterator for the reference position of the nodal coordinates
                const node_in X_end; //! Stopping iterator for the reference position of the nodal coordinates

                const local_node_in local_node_xi_begin; //! Starting iterator for the local nodal coordinates
                const local_node_in local_node_xi_end; //! Stopping iterator for the local nodal coordinates

        };

        //! An implementation of a linear hexahedral element
        template<typename T, class node_in, class local_point_in, class shape_functions_out, class grad_shape_functions_out>
        class LinearHex : public FiniteElementBase<3, node_in, typename std::array< T, 3 * 8 >::const_iterator, local_point_in, shape_functions_out, grad_shape_functions_out>{

            public:

                constexpr static std::array< T, 3 * 8 > local_nodes = { -1, -1, -1,
                                                                         1, -1, -1,
                                                                         1,  1, -1,
                                                                        -1,  1, -1,
                                                                        -1, -1,  1,
                                                                         1, -1,  1,
                                                                         1,  1,  1,
                                                                        -1,  1,  1 };

                LinearHex( const node_in &_x_begin, const node_in &_x_end, const node_in &_X_begin, const node_in &_X_end );

                using FiniteElementBase<dim,node_in,typename std::array< T, 3 * 8 >::const_iterator,local_point_in,shape_functions_out,grad_shape_functions_out>::FiniteElementBase;

                virtual void GetShapeFunctions( const local_point_in &xi_begin, const local_point_in &xi_end, shape_functions_out N_begin, shape_functions_out N_end ) override;

                virtual void GetLocalShapeFunctionGradients( const local_point_in &xi_begin, const local_point_in &xi_end, grad_shape_functions_out dNdxi_begin, grad_shape_functions_out dNdxi_end ) override;

        };

    }

}

#include "tardigrade_finite_element_utilities.cpp"

#endif
