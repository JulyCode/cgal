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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_EXPLICIT_CARTESIAN_GRID_FUNCTION_H
#define CGAL_ISOSURFACING_3_INTERNAL_EXPLICIT_CARTESIAN_GRID_FUNCTION_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Grid_topology_3.h>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <typename Grid, typename T>
class Cartesian_grid_3_view
{
public:
  using value_type = T;

  using Vertex_descriptor = typename Grid_topology_3::Vertex_descriptor;

public:
  Cartesian_grid_3_view(const Grid& grid)
    : grid{grid}
  { }

  // gets the value at vertex `v`
  value_type operator()(const Vertex_descriptor& v) const
  {
    return grid(v);
  }

private:
  const Grid& grid;
};

template<typename GeomTraits, typename Grid>
using Explicit_Cartesian_grid_function_3 = Cartesian_grid_3_view<Grid, typename GeomTraits::FT>;

template<typename GeomTraits, typename Grid>
using Explicit_Cartesian_grid_gradient_3 = Cartesian_grid_3_view<Grid, typename GeomTraits::Vector_3>;

template<typename GeomTraits, typename Grid>
using Explicit_Cartesian_grid_geometry_3 = Cartesian_grid_3_view<Grid, typename GeomTraits::Point_3>;

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_EXPLICIT_CARTESIAN_GRID_FUNCTION_H
