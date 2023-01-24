// Copyright (c) 2022-2023 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_ISOSURFACING_3_EXPLICIT_CURVED_CARTESIAN_GRID_DOMAIN_3_H
#define CGAL_ISOSURFACING_3_EXPLICIT_CURVED_CARTESIAN_GRID_DOMAIN_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Isosurfacing_domain_3.h>
#include <CGAL/Isosurfacing_3/internal/Grid_topology_3.h>
#include <CGAL/Isosurfacing_3/internal/Cartesian_grid_3_view.h>
#include <CGAL/Isosurfacing_3/Zero_gradient.h>

#include <CGAL/assertions.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domains_grp
 *
 * \cgalModels `IsosurfacingDomainWithGradient_3`
 *
 * \brief A domain that represents a %Cartesian grid that discretizes an implicit function.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 * \tparam ImplicitFunction the type of the implicit function. It must be a model of `CopyConstructible` and implement
 *                          `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \sa `CGAL::Isosurfacing::create_implicit_Cartesian_grid_domain()`
 */
#ifdef DOXYGEN_RUNNING // Otherwise it shows what is behind "using" in the doc...
template <typename GeomTraits,
          typename ImplicitFunction,
          typename Gradient = Zero_gradient>
using Implicit_Cartesian_grid_domain_3 = unspecified_type;
#else
template <typename GeomTraits,
          typename PositionGrid,
          typename FunctionGrid,
          typename Gradient = Zero_gradient,
          typename Topology = internal::Grid_topology_3,
          typename Geometry = internal::Explicit_Cartesian_grid_geometry_3<GeomTraits, PositionGrid>,
          typename Function = internal::Explicit_Cartesian_grid_function_3<GeomTraits, FunctionGrid>>
using Explicit_curved_Cartesian_grid_domain_3 =
  internal::Isosurfacing_domain_3<GeomTraits,
                                  Topology,
                                  Geometry,
                                  Function,
                                  Gradient>;
#endif

/**
 * \ingroup IS_Domains_grp
 *
 * \brief creates a domain from an implicit function that can be used as input for isosurfacing algorithms.
 *
 * \details The implicit function is evaluated on the grid vertices of the virtual grid
 * defined by the bounding box and the spacing value. By not storing any function values implicitly,
 * fewer memory accesses are required in comparison to an `Explicit_Cartesian_grid_domain_3`.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 * \tparam ImplicitFunction the type of the implicit function. It must be a model of `CopyConstructible` and implement
 *                          `GeomTraits::FT operator()(const GeomTraits::Point_3& point) const`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \param bbox a bounding box that specifies the dimensions of the implicit function's domain
 * \param spacing the distance between discretization points
 * \param point_function the implicit function
 * \param grad a function that describes the gradient of the data
 * \param gt an instance of geometric traits
 *
 * \return a new `Implicit_Cartesian_grid_domain_3`
 */
template <typename GeomTraits,
          typename PositionGrid,
          typename FunctionGrid,
          typename Gradient = Zero_gradient>
Explicit_curved_Cartesian_grid_domain_3<GeomTraits, PositionGrid, FunctionGrid, Gradient>
create_explicit_curved_Cartesian_grid_domain(const PositionGrid& vertex_positions,
                                             const FunctionGrid& function,
                                             const Gradient& grad = Gradient(),
                                             const GeomTraits& gt = GeomTraits())
{
  using Domain = Explicit_curved_Cartesian_grid_domain_3<GeomTraits, PositionGrid, FunctionGrid, Gradient>;

  using Topology = typename Domain::Topology;
  using Geometry = typename Domain::Geometry;
  using Function = typename Domain::Function;

  const std::size_t size_i = vertex_positions.xdim();
  const std::size_t size_j = vertex_positions.ydim();
  const std::size_t size_k = vertex_positions.zdim();

  const Topology topo { size_i, size_j, size_k };
  const Geometry geom { vertex_positions };
  const Function func { function };

  return Domain{ topo, geom, func, grad, gt };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_EXPLICIT_CURVED_CARTESIAN_GRID_DOMAIN_3_H
