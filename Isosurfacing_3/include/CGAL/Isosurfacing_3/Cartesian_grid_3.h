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

#ifndef CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H
#define CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Grid_topology_3.h>
#include <CGAL/assertions.h>
#include <CGAL/Image_3.h>

#include <array>
#include <type_traits>
#include <vector>

namespace CGAL {
namespace Isosurfacing {

/**
 * \ingroup IS_Domains_grp
 *
 * \brief stores scalar values and gradients at the vertices of a %Cartesian grid.
 *
 * \tparam GeomTraits must be a model of `IsosurfacingTraits_3`.
 */
template <typename T>
class Cartesian_grid_3
{
public:
  using value_type = T;

  using Vertex_descriptor = typename internal::Grid_topology_3::Vertex_descriptor;

private:
  std::array<std::size_t, 3> m_sizes;

  std::vector<value_type> m_values;

public:
  /**
   * \brief creates a grid with `xdim * ydim * zdim` grid vertices.
   *
   * The grid covers the space described by a bounding box.
   *
   * \param xdim the number of grid vertices in the `x` direction
   * \param ydim the number of grid vertices in the `y` direction
   * \param zdim the number of grid vertices in the `z` direction
   * \param bbox the bounding box
   * \param gt the geometric traits
   *
   * \pre `xdim`, `ydim`, and `zdim` are (strictly) positive.
   */
  Cartesian_grid_3(const std::size_t xdim,
                   const std::size_t ydim,
                   const std::size_t zdim)
    : m_sizes{xdim, ydim, zdim}
  {
    CGAL_precondition(xdim > 0);
    CGAL_precondition(ydim > 0);
    CGAL_precondition(zdim > 0);

    // pre-allocate memory
    const std::size_t nv = xdim * ydim * zdim;
    m_values.resize(nv);
  }

  /**
   * \brief creates a grid from a `CGAL::Image_3`.
   *
   * The dimensions are read from the image. The values stored
   * in the image must be of type `Geom_traits::FT` or implicitly convertible to it.
   *
   * \param image the image providing the data
   */
  static Cartesian_grid_3<T> from_image(const Image_3& image);

  /**
   * \brief Reads the bounding box from a `CGAL::Image_3`.
   *
   * \param image the image providing the data
   */
  static Bbox_3 bbox_from_image(const Image_3& image);

  /**
   * \brief converts a grid to a `CGAL::Image_3`.
   *
   * The dimensions and bounding box are written to the image.
   * The values stored in this grid must be scalars.
   *
   * \param bbox the bounding box of the image
   */
  Image_3 to_image(const Bbox_3& bbox) const;

public:

  /**
   * \return the number of grid vertices in the `x` direction
   */
  std::size_t xdim() const
  {
    return m_sizes[0];
  }

  /**
   * \return the number of grid vertices in the `y` direction
   */
  std::size_t ydim() const
  {
    return m_sizes[1];
  }

  /**
   * \return the number of grid vertices in the `z` direction
   */
  std::size_t zdim() const
  {
    return m_sizes[2];
  }

  /**
   * \brief gets the value stored at the grid vertex described by a set of indices.
   *
   * \param x the index in the `x` direction
   * \param y the index in the `y` direction
   * \param z the index in the `z` direction
   *
   * \return the stored value
   */
  value_type operator()(const std::size_t x,
                        const std::size_t y,
                        const std::size_t z) const
  {
      return m_values[linear_index(x, y, z)];
  }

  /**
   * \brief gets the value stored at the grid vertex.
   *
   * \param v the vertex descriptor
   *
   * \return the stored value
   */
  value_type operator()(const Vertex_descriptor& v) const
  {
    return m_values[linear_index(v[0], v[1], v[2])];
  }

  /**
   * \brief gets the value stored at the grid vertex described by a set of indices.
   *
   * \note This function can be used to set the value at a grid vertex.
   *
   * \param x the index in the `x` direction
   * \param y the index in the `y` direction
   * \param z the index in the `z` direction
   *
   * \return a reference to the stored value
   */
  value_type& operator()(const std::size_t x,
                         const std::size_t y,
                         const std::size_t z)
  {
      return m_values[linear_index(x, y, z)];
  }

  /**
   * \brief gets the value stored at the grid vertex.
   *
   * \note This function can be used to set the value at a grid vertex.
   *
   * \param v the vertex descriptor
   *
   * \return a reference to the stored value
   */
  value_type& operator()(const Vertex_descriptor& v)
  {
    return m_values[linear_index(v[0], v[1], v[2])];
  }

private:
  std::size_t linear_index(const std::size_t x,
                           const std::size_t y,
                           const std::size_t z) const
  {
    CGAL_precondition(x < xdim() && y < ydim() && z < zdim());

    // convert (x, y, z) into a linear index to access the scalar values / gradient vectors
    return (z * ydim() + y) * xdim() + x;
  }
};

template<typename T>
Cartesian_grid_3<T>
Cartesian_grid_3<T>::
from_image(const Image_3& image)
{
  // create grid
  Cartesian_grid_3<T> grid{ image.xdim(), image.ydim(), image.zdim() };

  const value_type* data = static_cast<const value_type*>(image.data());
  // copy values
  for (std::size_t x = 0; x < grid.xdim(); ++x)
  {
    for (std::size_t y = 0; y < grid.ydim(); ++y)
    {
      for (std::size_t z = 0; z < grid.zdim(); ++z)
      {
        const std::size_t lid = grid.linear_index(x, y, z);
        grid(x, y, z) = data[lid];
      }
    }
  }
  return grid;
}

template<typename T>
Bbox_3
Cartesian_grid_3<T>::bbox_from_image(const Image_3& image)
{
  // compute bounding box
  const FT max_x = image.tx() + (image.xdim() - 1) * image.vx();
  const FT max_y = image.ty() + (image.ydim() - 1) * image.vy();
  const FT max_z = image.tz() + (image.zdim() - 1) * image.vz();
  return Bbox_3{ image.tx(), image.ty(), image.tz(), max_x, max_y, max_z };
}

template <typename T>
Image_3
Cartesian_grid_3<T>::
to_image(const Bbox_3& bbox) const
{
  // select number type
  WORD_KIND wordkind;
  if(std::is_floating_point<value_type>::value)
    wordkind = WK_FLOAT;
  else
    wordkind = WK_FIXED;

  // select signed or unsigned
  SIGN sign;
  if (std::is_signed<value_type>::value)
    sign = SGN_SIGNED;
  else
    sign = SGN_UNSIGNED;

  // get spacing
  const double vx = bbox.x_span() / (xdim() - 1);
  const double vy = bbox.x_span() / (ydim() - 1);
  const double vz = bbox.x_span() / (zdim() - 1);

  // create image
  _image* im = _createImage(xdim(), ydim(), zdim(),
                            1,           // vectorial dimension
                            vx, vy, vz,  // voxel size
                            sizeof(value_type),  // image word size in bytes
                            wordkind,    // image word kind WK_FIXED, WK_FLOAT, WK_UNKNOWN
                            sign);       // image word sign

  // error handling
  if(im == nullptr || im->data == nullptr)
    throw std::bad_alloc();  // @todo idk?

  // set min coordinates
  im->tx = bbox.xmin();
  im->ty = bbox.ymin();
  im->tz = bbox.zmin();

  // copy data
  value_type* data = static_cast<value_type*>(im->data);
  for(std::size_t x=0; x<xdim(); ++x) {
    for(std::size_t y=0; y<ydim(); ++y) {
      for(std::size_t z=0; z<zdim(); ++z)
      {
        const std::size_t lid = linear_index(x, y, z);
        data[lid] = operator(x, y, z);
      }
    }
  }

  return Image_3{ im, Image_3::OWN_THE_DATA };
}

} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_CARTESIAN_GRID_3_H
