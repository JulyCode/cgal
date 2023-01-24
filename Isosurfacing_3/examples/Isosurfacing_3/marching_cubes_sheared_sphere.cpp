#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Explicit_curved_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>

#include <CGAL/Bbox_3.h>

#include <CGAL/boost/graph/IO/OFF.h>

#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;

using FunctionGrid = CGAL::Isosurfacing::Cartesian_grid_3<FT>;
using PositionGrid = CGAL::Isosurfacing::Cartesian_grid_3<Point>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

int main(int, char**)
{
  FunctionGrid function{ 50, 50, 50 };
  PositionGrid positions{ 50, 50, 50 };

  // create a domain with given bounding box and grid spacing
  auto domain = CGAL::Isosurfacing::create_explicit_curved_Cartesian_grid_domain<Kernel>(positions, function);
  using Domain = decltype(domain);

  // compute and store function values at all grid points
  auto init_grid = [&](const Domain::Vertex_descriptor& v)
  {
    const Point p{ (int) v[0] - 25, (int) v[1] - 25, (int) v[2] - 25 };
    function(v) = sqrt((p - CGAL::ORIGIN).squared_length());

    positions(v) = Point{ p[0] + p[1], p[1], p[2] }; // shear x direction
  };

  domain.iterate_vertices(init_grid);

  // prepare collections for the output indexed mesh
  Point_range points;
  Polygon_range polygons;

  // execute marching cubes with an isovalue of 20
  CGAL::Isosurfacing::marching_cubes(domain, 20, points, polygons);

  // save ouput indexed mesh to a file, in the OFF format
  CGAL::IO::write_OFF("result.off", points, polygons);

  return EXIT_SUCCESS;
}
