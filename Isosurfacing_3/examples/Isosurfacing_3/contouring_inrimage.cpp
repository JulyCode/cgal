#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_gradients_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>

#include <CGAL/Image_3.h>

#include <CGAL/Isosurfacing_3/IO/Image_3.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <vector>

#include <iostream>
#include <fstream>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Interpolated_discrete_values_3<Grid>;
using Gradients = CGAL::Isosurfacing::Interpolated_discrete_gradients_3<Grid>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

namespace IS = CGAL::Isosurfacing;

void run_marching_cubes(const Grid& grid,
                        const FT isovalue,
                        const Values& values)
{
  using Domain = IS::Marching_cubes_domain_3<Grid, Values, IS::Linear_interpolation_edge_intersection>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Marching Cubes with isovalue = " << isovalue << std::endl;

  // fill up values

  // create a domain from the grid
  Domain domain { grid, values };

  // prepare collections for the output indexed soup
  Point_range points;
  Polygon_range triangles;

  // execute marching cubes
  IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles);

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #triangles: " << triangles.size() << std::endl;

  // save output indexed mesh to a file, in the OFF format
  CGAL::IO::write_polygon_soup("marching_cubes_inrimage.off", points, triangles);
}

void run_dual_contouring(const Grid& grid,
                         const FT isovalue,
                         const Values& values)
{
  using Domain = IS::Dual_contouring_domain_3<Grid, Values, Gradients, IS::Linear_interpolation_edge_intersection>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Dual Contouring with isovalue = " << isovalue << std::endl;

  // fill up values and gradients
  const FT step = CGAL::approximate_sqrt(grid.spacing().squared_length()) * 0.01; // finite difference step
  // Gradients gradients { values, step };
  Gradients gradients { grid };
  gradients.compute_discrete_gradients_from_discrete_values(values);
  Domain domain { grid, values, gradients };

  Point_range debug_points;
  Polygon_range debug_faces;
  auto debug_grid_creator = [&](const Domain::cell_descriptor& c)
    {
      std::vector<std::size_t> cell_vertices;
      for (const auto& v : domain.cell_vertices(c)) {
        cell_vertices.push_back(debug_points.size());
        debug_points.push_back(domain.point(v));
      }
      debug_faces.push_back({cell_vertices[0], cell_vertices[1], cell_vertices[3], cell_vertices[2]});
      debug_faces.push_back({cell_vertices[5], cell_vertices[1], cell_vertices[3], cell_vertices[7]});
      debug_faces.push_back({cell_vertices[4], cell_vertices[5], cell_vertices[1], cell_vertices[0]});
      debug_faces.push_back({cell_vertices[4], cell_vertices[0], cell_vertices[2], cell_vertices[6]});
      debug_faces.push_back({cell_vertices[4], cell_vertices[5], cell_vertices[7], cell_vertices[6]});
      debug_faces.push_back({cell_vertices[6], cell_vertices[7], cell_vertices[3], cell_vertices[2]});
    };
    domain.template for_each_cell<CGAL::Sequential_tag>(debug_grid_creator);
  CGAL::IO::write_polygon_soup("dual_contouring_inrimage_debug_grid.off", debug_points, debug_faces);

  Point_range points;
  Polygon_range triangles;

  // run dual contouring isosurfacing
  IS::dual_contouring<CGAL::Sequential_tag>(domain, isovalue, points, triangles);

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #triangles: " << triangles.size() << std::endl;
  CGAL::IO::write_polygon_soup("dual_contouring_inrimage.off", points, triangles);
}

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : "../examples/Isosurfacing_3/FullHead.inr";//CGAL::data_file_path("images/skull_2.9.inr");//
  const FT isovalue = (argc > 2) ? std::stod(argv[2]) : 1139;//2.9;//

  // load volumetric image from a file
  CGAL::Image_3 image;
  if(!image.read(fname))
  {
    std::cerr << "Error: Cannot read image file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // convert image to a Cartesian grid
  Grid grid;
  Values values { grid }; // 'values' keeps a reference to the grid
  if(!IS::IO::convert_image_to_grid(image, grid, values))
  {
    std::cerr << "Error: Cannot convert image to Cartesian grid" << std::endl;
    return EXIT_FAILURE;
  }

  std::ofstream file("full_head_half_res.txt");
  file << grid.xdim()/2 << " " << grid.ydim()/2 << " " << grid.zdim()/2 << std::endl;
  file << grid.spacing()[0]*2 << " " << grid.spacing()[1]*2 << " " << grid.spacing()[2]*2 << std::endl;

  for (int x = 0; x < grid.xdim(); x+=2) {
    for (int y = 0; y < grid.ydim(); y+=2) {
      for (int z = grid.zdim() - 1; z >= 0; z-=2) {
      // for (int z = 0; z < grid.zdim(); z++) {
        auto tmp = values(x, y, z);
        values(x, y, z) = 2 * isovalue - values(x, y, z);
        auto tmp2 = values(x, y, z);
        auto tmp3 = values(x, y, z);

        file << values(x, y, z) - isovalue << "\n";
      }
    }
  }
  file.close();

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  run_marching_cubes(grid, isovalue, values);

  // run_dual_contouring(grid, isovalue, values);

  return EXIT_SUCCESS;
}
