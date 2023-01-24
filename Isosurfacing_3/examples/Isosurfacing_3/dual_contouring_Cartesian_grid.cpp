#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Explicit_Cartesian_grid_gradient_3.h>
#include <CGAL/Isosurfacing_3/Explicit_Cartesian_grid_domain_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>

#include <CGAL/boost/graph/IO/OFF.h>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Vector = typename Kernel::Vector_3;

using FunctionGrid = CGAL::Isosurfacing::Cartesian_grid_3<FT>;
using GradientGrid = CGAL::Isosurfacing::Cartesian_grid_3<Vector>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

int main(int, char**)
{
  // create bounding box and grid
  const CGAL::Bbox_3 bbox{-1., -1., -1.,  1., 1., 1.};
  FunctionGrid func_grid{ 30, 30, 30 };
  GradientGrid grad_grid{ 30, 30, 30 };

  // gradient field
  CGAL::Isosurfacing::Explicit_Cartesian_grid_gradient_3<Kernel, GradientGrid> gradient(grad_grid, bbox);

  // create domain from scalar and gradient fields
  auto domain = CGAL::Isosurfacing::create_explicit_Cartesian_grid_domain<Kernel>(bbox, func_grid, gradient);
  using Domain = decltype(domain);

  // compute field values and gradients
  auto init_grid = [&](const Domain::Vertex_descriptor& v)
  {
    const Vector direction = domain.point(v) - CGAL::ORIGIN;
    const FT distance = CGAL::approximate_sqrt(direction.squared_length());

    func_grid(v) = distance;
    grad_grid(v) = direction / distance; // @todo check division / 0
  };

  domain.iterate_vertices(init_grid);

  Point_range points;
  Polygon_range polygons;

  // run dual contouring isosurfacing
  CGAL::Isosurfacing::dual_contouring(domain, 0.8, points, polygons);

  // write output indexed surface mesh to file, in OFF format
  CGAL::IO::write_OFF("result.off", points, polygons);

  return EXIT_SUCCESS;
}
