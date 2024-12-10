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
        template<int dim, int local_dim, int node_count, class node_in, class local_node_in, class local_point_in, class shape_functions_out, class grad_shape_functions_out>
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

                virtual void GetGlobalShapeFunctionGradients( const local_point_in &xi_begin, const local_point_in &xi_end,
                                                              const node_in &node_positions_begin, const node_in &node_positions_end,
                                                              grad_shape_functions_out value_begin, grad_shape_functions_out value_end ){
                    /*!
                     * Compute the global gradient of the shape functions
                     *
                     * \param &xi_begin: The starting iterator of the local point
                     * \param &xi_end: The stopping iterator of the local point
                     * \param &node_positions_begin: The starting iterator of the nodal positions (row major)
                     * \param &node_positions_end: The stopping iterator of the nodal positions (row major)
                     * \param &value_begin: The starting iterator of the shape function global gradient (row major)
                     * \param &value_end: The stopping iterator of the shape function global gradient (row major)
                     */

                    throw std::logic_error("Function not implemented");
                     
                }

                virtual void GetVolumeIntegralJacobianOfTransformation( const local_point_in &xi_begin, const local_point_in &xi_end,
                                                                        typename std::iterator_traits<node_in>::value_type &value, const bool configuration = 1 ){
                    /*!
                     * Compute the value of the Jacobian of transformation from the local coordinates to the configuration for volume integrals
                     * 
                     * \param &xi_begin: The starting iterator of the local coordinates
                     * \param &xi_end: The stopping iterator of the local coordinates
                     * \param value: The Jacobian of transformation going from the local coordinates to the indicated configuration
                     * \param configuration: Compute the gradient w.r.t. the current configuration ( true ) or reference configuration ( false )
                     */

                    throw std::logic_error("Function not implemented");
                     
                }

                template<class quantity_in, class quantity_out>
                void InterpolateQuantity( const local_point_in &xi_begin, const local_point_in &xi_end, const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                          quantity_out value_begin, quantity_out value_end ){
                    /*!
                     * Interpolate the provided quantity to the local point
                     * 
                     * \param &xi_begin: The starting iterator of the local coordinates
                     * \param &xi_end: The stopping iterator of the local coordinates
                     * \param &quantity_begin: The starting iterator of the quantity at the nodes (row-major)
                     * \param &quantity_end: The stopping iterator of the quantity at the nodes (row-major)
                     * \param value_begin: The starting iterator for the interpolated value
                     * \param value_end: The stopping iterator for the interpolated value
                     */

                    const size_type quantity_dim = ( size_type )( value_end - value_begin );

                    TARDIGRADE_ERROR_TOOLS_CHECK( quantity_dim * node_count == ( size_type )( quantity_end - quantity_begin ), "The returned value size (" + std::to_string( quantity_dim ) + ") and the quantity dimension (" + std::to_string( ( size_type )( quantity_end - quantity_begin ) ) + ") are inconsistent with the node count (" + std::to_string( node_count ) + ")" );

                    GetShapeFunctions( xi_begin, xi_end, std::begin( _shapefunctions ), std::end( _shapefunctions ) );

                    std::fill( value_begin, value_end, 0 );

                    for ( auto N = std::begin( _shapefunctions ); N != std::end( _shapefunctions ); ++N ){

                        unsigned int node = ( size_type )( N - std::begin( _shapefunctions ) );

                        for ( auto v = value_begin; v != value_end; ++v ){

                            *v += ( *N ) * ( *( quantity_begin + quantity_dim * node + ( size_type )( v - value_begin ) ) );

                        }

                    }

                }

                template<class quantity_in, class quantity_gradient_out>
                void GetLocalQuantityGradient( const local_point_in &xi_begin, const local_point_in &xi_end, const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                               quantity_gradient_out value_begin, quantity_gradient_out value_end ){
                    /*!
                     * Compute the gradient of the quantity at a local point
                     * 
                     * \param &xi_begin: The starting iterator of the local coordinates
                     * \param &xi_end: The stopping iterator of the local coordinates
                     * \param &quantity_begin: The starting iterator of the quantity at the nodes (row-major)
                     * \param &quantity_end: The stopping iterator of the quantity at the nodes (row-major)
                     * \param value_begin: The starting iterator for the computed local gradient of the quantity
                     * \param value_end: The stopping iterator for the computed local gradient of the quantity in row-major form
                     */

                    const size_type quantity_dim = ( size_type )( value_end - value_begin ) / local_dim;

                    TARDIGRADE_ERROR_TOOLS_CHECK( quantity_dim * node_count == ( size_type )( quantity_end - quantity_begin ), "The returned value size (" + std::to_string( quantity_dim ) + ") and the quantity dimension (" + std::to_string( ( size_type )( quantity_end - quantity_begin ) ) + ") are inconsistent with the node count (" + std::to_string( node_count ) + ")" );

                    TARDIGRADE_ERROR_TOOLS_CATCH( GetLocalShapeFunctionGradients( xi_begin, xi_end, std::begin( _local_gradshapefunctions ), std::end( _local_gradshapefunctions ) ) );

                    std::fill( value_begin, value_end, 0 );

                    for ( unsigned int node = 0; node < node_count; ++node ){

                        for ( unsigned int row = 0; row < quantity_dim; ++row ){

                            for ( unsigned int col = 0; col < local_dim; ++col ){

                                *( value_begin + local_dim * row + col ) += _local_gradshapefunctions[ local_dim * node + col ] * ( *( quantity_begin + quantity_dim * node + row ) );

                            }

                        }

                    }

                }

                template<class quantity_in, class quantity_gradient_out>
                void GetGlobalQuantityGradient( const local_point_in &xi_begin, const local_point_in &xi_end, const quantity_in &quantity_begin, const quantity_in &quantity_end,
                                               quantity_gradient_out value_begin, quantity_gradient_out value_end, const bool configuration = true ){
                    /*!
                     * Compute the global gradient of the quantity at a local point
                     * 
                     * \param &xi_begin: The starting iterator of the local coordinates
                     * \param &xi_end: The stopping iterator of the local coordinates
                     * \param &quantity_begin: The starting iterator of the quantity at the nodes (row-major)
                     * \param &quantity_end: The stopping iterator of the quantity at the nodes (row-major)
                     * \param value_begin: The starting iterator for the computed global gradient of the quantity
                     * \param value_end: The stopping iterator for the computed global gradient of the quantity in row-major form
                     * \param configuration: Compute the gradient w.r.t. the current configuration ( true ) or reference configuration ( false )
                     */

                    const size_type quantity_dim = ( size_type )( value_end - value_begin ) / dim;

                    TARDIGRADE_ERROR_TOOLS_CHECK( quantity_dim * node_count == ( size_type )( quantity_end - quantity_begin ), "The returned value size (" + std::to_string( quantity_dim ) + ") and the quantity dimension (" + std::to_string( ( size_type )( quantity_end - quantity_begin ) ) + ") are inconsistent with the node count (" + std::to_string( node_count ) + ")" );

                    if ( configuration ){

                        TARDIGRADE_ERROR_TOOLS_CATCH( GetGlobalShapeFunctionGradients( xi_begin, xi_end, x_begin, x_end, std::begin( _global_gradshapefunctions ), std::end( _global_gradshapefunctions ) ) );

                    }
                    else{

                        TARDIGRADE_ERROR_TOOLS_CATCH( GetGlobalShapeFunctionGradients( xi_begin, xi_end, X_begin, X_end, std::begin( _global_gradshapefunctions ), std::end( _global_gradshapefunctions ) ) );

                    }

                    std::fill( value_begin, value_end, 0 );

                    for ( unsigned int node = 0; node < node_count; ++node ){

                        for ( unsigned int row = 0; row < quantity_dim; ++row ){

                            for ( unsigned int col = 0; col < local_dim; ++col ){

                                *( value_begin + local_dim * row + col ) += _global_gradshapefunctions[ local_dim * node + col ] * ( *( quantity_begin + quantity_dim * node + row ) );

                            }

                        }

                    }

                }

            protected:

                const node_in x_begin; //!< Starting iterator for the current position of the nodal coordinates
                const node_in x_end; //!< Stopping iterator for the current position of the nodal coordinates

                const node_in X_begin; //!< Starting iterator for the reference position of the nodal coordinates
                const node_in X_end; //!< Stopping iterator for the reference position of the nodal coordinates

                const local_node_in local_node_xi_begin; //!< Starting iterator for the local nodal coordinates
                const local_node_in local_node_xi_end; //!< Stopping iterator for the local nodal coordinates

                std::array<typename std::iterator_traits<shape_functions_out>::value_type, node_count> _shapefunctions; //!< A temporary storage container for shapefunction values

                std::array<typename std::iterator_traits<shape_functions_out>::value_type, local_dim * node_count> _local_gradshapefunctions; //!< A temporary storage container for local grad shapefunction values

                std::array<typename std::iterator_traits<shape_functions_out>::value_type, dim * node_count> _global_gradshapefunctions; //!< A temporary storage container for global grad shapefunction values

        };

        //! An implementation of a linear hexahedral element
        template<typename T, class node_in, class local_point_in, class shape_functions_out, class grad_shape_functions_out>
        class LinearHex : public FiniteElementBase<3, 3, 8, node_in, typename std::array< T, 3 * 8 >::const_iterator, local_point_in, shape_functions_out, grad_shape_functions_out>{

            public:

                //! The local nodes for an isoparametric linear hex element
                constexpr static std::array< T, 3 * 8 > local_nodes = { -1, -1, -1,
                                                                         1, -1, -1,
                                                                         1,  1, -1,
                                                                        -1,  1, -1,
                                                                        -1, -1,  1,
                                                                         1, -1,  1,
                                                                         1,  1,  1,
                                                                        -1,  1,  1 };

                LinearHex( const node_in &_x_begin, const node_in &_x_end, const node_in &_X_begin, const node_in &_X_end );

                using FiniteElementBase<3,3,8,node_in,typename std::array< T, 3 * 8 >::const_iterator,local_point_in,shape_functions_out,grad_shape_functions_out>::FiniteElementBase;

                virtual void GetShapeFunctions( const local_point_in &xi_begin, const local_point_in &xi_end, shape_functions_out N_begin, shape_functions_out N_end ) override;

                virtual void GetLocalShapeFunctionGradients( const local_point_in &xi_begin, const local_point_in &xi_end, grad_shape_functions_out dNdxi_begin, grad_shape_functions_out dNdxi_end ) override;

                virtual void GetGlobalShapeFunctionGradients( const local_point_in &xi_begin, const local_point_in &xi_end,
                                                              const node_in &node_positions_begin, const node_in &node_positions_end,
                                                              grad_shape_functions_out value_begin, grad_shape_functions_out value_end ) override;

                virtual void GetVolumeIntegralJacobianOfTransformation( const local_point_in &xi_begin, const local_point_in &xi_end,
                                                                        typename std::iterator_traits<node_in>::value_type &value, const bool configuration = 1 ) override;

        };

    }

}

#include "tardigrade_finite_element_utilities.cpp"

#endif
