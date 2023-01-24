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

#ifndef CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_DOMAIN_3_H
#define CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_DOMAIN_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Isosurfacing_domain_3.h>
#include <CGAL/Isosurfacing_3/internal/Cartesian_grid_3_view.h>
#include <CGAL/Isosurfacing_3/internal/Implicit_Cartesian_grid_geometry_3.h>
#include <CGAL/Isosurfacing_3/internal/Grid_topology_3.h>
#include <CGAL/Isosurfacing_3/Zero_gradient.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domains_grp
 *
 * \cgalModels `IsosurfacingDomainWithGradient_3`
 *
 * \brief A domain that represents an explicitly stored %Cartesian grid.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \sa `CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain()`
 */
#ifdef DOXYGEN_RUNNING // Allow more than a Cartesian_grid_3
template <template <typename GeomTraits> class CGAL::Isosurfacing::Cartesian_grid_3,
          typename Gradient = Zero_gradient>
using Explicit_Cartesian_grid_domain_3 = unspecified_type;
#else
template <typename GeomTraits,
          typename Grid,
          typename Gradient = Zero_gradient,
          typename Topology = internal::Grid_topology_3,
          typename Geometry = internal::Implicit_Cartesian_grid_geometry_3<GeomTraits>,
          typename Function = internal::Explicit_Cartesian_grid_function_3<GeomTraits, Grid>>
using Explicit_Cartesian_grid_domain_3 =
  internal::Isosurfacing_domain_3<GeomTraits,
                                  Topology,
                                  Geometry,
                                  Function,
                                  Gradient>;
#endif

/**
 * \ingroup IS_Domains_grp
 *
 * \brief Creates a domain that can be used as input for isosurfacing algorithms.
 *
 * \warning As the domain will keep a pointer to the `grid` object, users must ensure that
 *          the lifetime of the `grid` object exceeds that of the object returned by this function.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 * \tparam Gradient the type of the gradient functor. It must be a model of `CopyConstructible` and implement
 *                  `GeomTraits::Vector_3 operator()(const GeomTraits::Point_3& point) const`.
 *
 * \param grid the %Cartesian grid containing input data
 * \param grad a function that describes the gradient of the data
 *
 * \return a new `Explicit_Cartesian_grid_domain_3`
 */
#ifdef DOXYGEN_RUNNING // Allow more than Cartesian_grid_3
template <typename GeomTraits,
          typename Gradient = Zero_gradient>
Explicit_Cartesian_grid_domain_3<Grid, Gradient>
create_explicit_Cartesian_grid_domain(const CGAL::Isosurfacing::Cartesian_grid_3<GeomTraits>& grid,
                                      const Gradient& grad = Gradient())
#else
// Actual code enables passing more than just a Cartesian_grid_3
template <typename GeomTraits,
          typename Grid,
          typename Gradient = Zero_gradient>
Explicit_Cartesian_grid_domain_3<GeomTraits, Grid, Gradient>
create_explicit_Cartesian_grid_domain(const Bbox_3& bbox,
                                      const Grid& grid,
                                      const Gradient& grad = Gradient(),
                                      const GeomTraits& gt = GeomTraits())
#endif
{
  using Domain = Explicit_Cartesian_grid_domain_3<GeomTraits, Grid, Gradient>;

  using Topology = typename Domain::Topology;
  using Geometry = typename Domain::Geometry;
  using Function = typename Domain::Function;

  auto vector = gt.construct_vector_3_object();

  const std::size_t size_i = grid.xdim();
  const std::size_t size_j = grid.ydim();
  const std::size_t size_k = grid.zdim();

  const typename Geometry::Vector_3 offset = vector(bbox.xmin(), bbox.ymin(), bbox.zmin());

  // calculate grid spacing
  const FT d_x = FT{ bbox.x_span() } / (size_i - 1);
  const FT d_y = FT{ bbox.y_span() } / (size_j - 1);
  const FT d_z = FT{ bbox.z_span() } / (size_k - 1);
  const typename Geometry::Vector_3 spacing = vector(d_x, d_y, d_z);

  const Topology topo { size_i, size_j, size_k };
  const Geometry geom { offset, spacing };
  const Function func { grid };

  return Domain{ topo, geom, func, grad, gt };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_EXPLICIT_CARTESIAN_GRID_DOMAIN_3_H
