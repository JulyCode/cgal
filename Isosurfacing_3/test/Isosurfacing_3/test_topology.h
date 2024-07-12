#ifndef CGAL_ISOSURFACING_TEST_TOPOLOGY_H
#define CGAL_ISOSURFACING_TEST_TOPOLOGY_H

#include "test_util.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>

#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/Random.h>

#include <cassert>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

namespace IS = CGAL::Isosurfacing;

template<typename T>
class Union_find
{
private:
    struct Entry
    {
        T ele;
        std::size_t rank;
    };

public:
  Union_find() : m_num_sets(0) {}

  void insert(const T& element)
  {
    if (m_parent.find(element) == m_parent.end())
    {
      m_parent[element] = {element, 0};
      m_num_sets++;
    }
  }

  const T& find(const T& element) const
  {
    T& parent = m_parent[element].ele;
    if (parent != element)
      parent = find(parent);  // path compression

    return parent;
  }

  T& find(const T& element)
  {
    T& parent = m_parent[element].ele;
    if (parent != element)
      parent = find(parent);  // path compression

    return parent;
  }

  void union_sets(const T& element1, const T& element2)
  {
    const T& root1 = find(element1);
    const T& root2 = find(element2);
    Entry& entry1 = m_parent[root1];
    Entry& entry2 = m_parent[root2];
    
    if (root1 != root2)
    {
      // Union by rank
      if (entry1.rank > entry2.rank)
        entry2.ele = root1;
      else if (entry1.rank < entry2.rank)
        entry1.ele = root2;
      else
      {
        entry2.ele = root1;
        entry1.rank++;
      }
    }
    m_num_sets--;
  }

  std::size_t num_sets() const
  {
    return m_num_sets;
  }

private:
  mutable std::unordered_map<T, Entry> m_parent;
  std::size_t m_num_sets;
};

template <typename PointRange, typename PolygonRange>
std::size_t count_objects(PointRange points, PolygonRange polygons)
{
  Union_find<std::size_t> connected_sets;

  for (std::size_t i = 0; i < points.size(); i++)
    connected_sets.insert(i);

  for (const auto& face : polygons)
  {
    for (std::size_t i = 0; i < face.size(); i++)
      connected_sets.union_sets(face[i], face[(i + 1) % face.size()]);
  }

  return connected_sets.num_sets();
}

template <typename PolygonMesh>
std::size_t euler_characteristic(const PolygonMesh& mesh)
{
  return mesh.number_of_vertices() - mesh.number_of_edges() + mesh.number_of_faces();
}

template <typename Grid>
bool check_closed_not_empty(const Grid& grid, const IS::Interpolated_discrete_values_3<Grid>& values, const typename Grid::Geom_traits::FT iso = 0)
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  FT boundary_min = std::numeric_limits<FT>::max();
  FT total_min = std::numeric_limits<FT>::max();
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        total_min = std::min(total_min, values(x, y, z));
        if (x == 0 || y == 0 || z == 0 || x == grid.xdim() - 1 || y == grid.ydim() - 1 || z == grid.zdim() - 1)
          boundary_min = std::min(boundary_min, values(x, y, z));
      }
    }
  }

  // empty cell
  if (total_min > iso)
  {
    return false;
  }

  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        const bool is_boundary = x == 0 || y == 0 || z == 0 || x == grid.xdim() - 1 || y == grid.ydim() - 1 || z == grid.zdim() - 1;

        // closed field
        if (!is_boundary && boundary_min - values(x, y, z) < 1e-6)
          return false;
      }
    }
  }
  return true;
}

template <typename Grid>
bool check_iso_vertices(const Grid& grid, const IS::Interpolated_discrete_values_3<Grid>& values, const typename Grid::Geom_traits::FT iso = 0) {
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        if (std::abs(values(x, y, z)) < iso + 1e-6)
          return true;
      }
    }
  }
  return false;
}

template <typename Grid>
bool check_equal_edges(const Grid& grid, const IS::Interpolated_discrete_values_3<Grid>& values) {
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        if (x < grid.xdim() - 1 && std::abs(values(x, y, z) - values(x + 1, y, z)) < 1e-6)
          return true;
        if (y < grid.ydim() - 1 && std::abs(values(x, y, z) - values(x, y + 1, z)) < 1e-6)
          return true;
        if (z < grid.zdim() - 1 && std::abs(values(x, y, z) - values(x, y, z + 1)) < 1e-6)
          return true;
      }
    }
  }
  return false;
}

template <typename Grid>
bool check_iso_edges(const Grid& grid, const IS::Interpolated_discrete_values_3<Grid>& values, const typename Grid::Geom_traits::FT iso = 0) {
  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        if (x < grid.xdim() - 1 && std::abs(values(x, y, z)) < iso + 1e-6 && std::abs(values(x + 1, y, z)) < iso + 1e-6)
          return true;
        if (y < grid.ydim() - 1 && std::abs(values(x, y, z)) < iso + 1e-6 && std::abs(values(x, y + 1, z)) < iso + 1e-6)
          return true;
        if (z < grid.zdim() - 1 && std::abs(values(x, y, z)) < iso + 1e-6 && std::abs(values(x, y, z + 1)) < iso + 1e-6)
          return true;
      }
    }
  }
  return false;
}

template <typename Grid>
void generate_random(const Grid& grid, IS::Interpolated_discrete_values_3<Grid>& values)
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  CGAL::Random rand;

  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        if (x == 0 || y == 0 || z == 0 || x == grid.xdim() - 1 || y == grid.ydim() - 1 || z == grid.zdim() - 1)
          values(x, y, z) = rand.uniform_real<FT>(1, 2);
        else
          values(x, y, z) = rand.uniform_real<FT>(-1, 1);
      }
    }
  }
}

template <typename Grid>
void generate_random_iso_vertices(const Grid& grid, IS::Interpolated_discrete_values_3<Grid>& values)
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  CGAL::Random rand;

  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        if (x == 0 || y == 0 || z == 0 || x == grid.xdim() - 1 || y == grid.ydim() - 1 || z == grid.zdim() - 1)
          values(x, y, z) = rand.uniform_real<FT>(1, 2);
        else
        {
          const int config = rand.uniform_int(0, 3);
          if (config == 0)
            values(x, y, z) = rand.uniform_real<FT>(-1, -1e-6);
          else if (config == 1)
            values(x, y, z) = 0;
          else
            values(x, y, z) = rand.uniform_real<FT>(1e-6, 1);
        }
      }
    }
  }
}

template <typename Grid>
void generate_random_equal_edges(const Grid& grid, IS::Interpolated_discrete_values_3<Grid>& values)
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  CGAL::Random rand;

  const FT chance = 0.3;

  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        if (x == 0 || y == 0 || z == 0 || x == grid.xdim() - 1 || y == grid.ydim() - 1 || z == grid.zdim() - 1)
          values(x, y, z) = rand.uniform_real<FT>(1, 2);
        else if (x == 1 && y == 1 && z == 1)
          values(x, y, z) = rand.uniform_real<FT>(-1, 1);
        else
        {
          if (x > 1 && rand.uniform_01<FT>() < chance)
            values(x, y, z) = values(x - 1, y, z);
          else if (y > 1 && rand.uniform_01<FT>() < chance)
            values(x, y, z) = values(x, y - 1, z);
          else if (z > 1 && rand.uniform_01<FT>() < chance)
            values(x, y, z) = values(x, y, z - 1);
          else
            values(x, y, z) = rand.uniform_real<FT>(-1, 1);
        }
      }
    }
  }
}

template <typename Grid>
void generate_random_iso_edges(const Grid& grid, IS::Interpolated_discrete_values_3<Grid>& values)
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  CGAL::Random rand;

  const FT chance = 0.2;

  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        if (x == 0 || y == 0 || z == 0 || x == grid.xdim() - 1 || y == grid.ydim() - 1 || z == grid.zdim() - 1)
          values(x, y, z) = rand.uniform_real<FT>(1, 2);
        else if (x == 1 && y == 1 && z == 1)
          values(x, y, z) = rand.uniform_real<FT>(-1, 1);
        else
        {
          if (x > 1 && rand.uniform_01<FT>() < chance)
          {
            values(x, y, z) = 0;
            values(x - 1, y, z) = 0;
          }
          else if (y > 1 && rand.uniform_01<FT>() < chance)
          {
            values(x, y, z) = 0;
            values(x - 1, y, z) = 0;
          }
          else if (z > 1 && rand.uniform_01<FT>() < chance)
          {
            values(x, y, z) = 0;
            values(x - 1, y, z) = 0;
          }
          else
            values(x, y, z) = rand.uniform_real<FT>(-1, 1);
        }
      }
    }
  }
}

template <typename Grid>
void generate_random_ambiguous_faces(const Grid& grid, IS::Interpolated_discrete_values_3<Grid>& values)
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  CGAL::Random rand;

  const FT chance = 0.2;

  std::array<bool, 6> ambiguous_faces;
  for (int i = 0; i < 6; i++)
    ambiguous_faces[i] = rand.uniform_01<FT>() < chance;

  for(std::size_t x=0; x<grid.xdim(); ++x) {
    for(std::size_t y=0; y<grid.ydim(); ++y) {
      for(std::size_t z=0; z<grid.zdim(); ++z)
      {
        if (x == 0 || y == 0 || z == 0 || x == grid.xdim() - 1 || y == grid.ydim() - 1 || z == grid.zdim() - 1)
          values(x, y, z) = rand.uniform_real<FT>(1, 2);
        else if (x == 1 && y == 1 && z == 1)
          values(x, y, z) = rand.uniform_real<FT>(-1, 1);
        else
        {
          if (x > 1 && rand.uniform_01<FT>() < chance)
          {
            values(x, y, z) = 0;
            values(x - 1, y, z) = 0;
          }
          else if (y > 1 && rand.uniform_01<FT>() < chance)
          {
            values(x, y, z) = 0;
            values(x - 1, y, z) = 0;
          }
          else if (z > 1 && rand.uniform_01<FT>() < chance)
          {
            values(x, y, z) = 0;
            values(x - 1, y, z) = 0;
          }
          else
            values(x, y, z) = rand.uniform_real<FT>(-1, 1);
        }
      }
    }
  }
}


enum class AmbiguousCase {
  CONTOUR_6_VTS,
  CONTOUR_7_VTS_CONFIG_A,
  CONTOUR_7_VTS_CONFIG_B,
  CONTOUR_8_VTS_CONFIG_A,
  CONTOUR_8_VTS_CONFIG_B,
  CONTOUR_9_VTS_CONFIG_A,
  // CONTOUR_9_VTS_CONFIG_B,
  MC_4_TUNNEL,
  MC_7_TUNNEL,
  MC_10_TUNNEL,
  MC_12_TUNNEL,
  MC_13_TUNNEL,
  MC_13_12_VTS_CONFIG_A,
  MC_13_12_VTS_CONFIG_B,
  MC_13_SIMPLE
};

template <typename Grid>
typename Grid::Geom_traits::FT generate_predefined_inner_ambiguity(const Grid& grid, IS::Interpolated_discrete_values_3<Grid>& values, const AmbiguousCase topology_case) {
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  CGAL::Random rand;

  assert(grid.xdim() == 4);
  assert(grid.ydim() == 4);
  assert(grid.zdim() == 4);

  std::array<FT, 8> vals;
  FT iso;

  switch (topology_case) {
    case AmbiguousCase::CONTOUR_6_VTS:
      vals = {-0.960492426392903, 0.793207329559554, 0.916735525189067, -0.422761282626275, -0.934993247757551, -0.850129305868777, -0.0367116785741896, -0.656740699156587};
      iso = 0.0387;
      break;
    case AmbiguousCase::CONTOUR_7_VTS_CONFIG_A:
      vals = {10.2967816247556, 9.45145192686147, 9.54753711271687, 10.6482067822841, 9.81494966341055, 9.31168538578250, 9.80950580411527, 10.7451536262220};
      iso = 9.8588;
      break;
    case AmbiguousCase::CONTOUR_7_VTS_CONFIG_B:
      vals = {9.9998593195995547, 9.9993381282115549, 9.9979160205452544, 9.9986053863704142, 9.9999374908631235, 9.999424800002032, 9.9983922749132219, 9.999579324965488};
      iso = 9.9994608191478135;
      break;
    case AmbiguousCase::CONTOUR_8_VTS_CONFIG_A:
      vals = {0.454797708726920, 0.801330575352402, 0.649991492712356, -0.973974554763863,	-0.134171007607162,	-0.0844698148589140, -0.826313795402046, 0.433391503783462};
      iso = 0.0097;
      break;
    case AmbiguousCase::CONTOUR_8_VTS_CONFIG_B:
      vals = {9.9985934885536665, 9.9998695572230147, 9.9999045831713928, 9.999316745478131, 9.9986117521866866, 9.9998754368055813, 9.9999031760062458, 9.9992041920402936};
      iso = 9.99946;
      break;
    case AmbiguousCase::CONTOUR_9_VTS_CONFIG_A:
      vals = {-15.6504952739285, 2.90290077342601, 24.5454566157887, -24.5274127623786, 21.6741877710053, -4.49696327433901, -19.7891575872492, -15.5588482753161};
      iso = -1.8061;
      break;
    // case AmbiguousCase::CONTOUR_9_VTS_CONFIG_B:
    //   vals = {1763.0000000000000, 1052.0000000000000, 1815.0000000000000, 1050.0000000000000, 960.00000000000000, 1325.0000000000000, 1150.0000000000000, 1260.0000000000000};
    //   iso = 1233.5999999999999;
    //   break;
    case AmbiguousCase::MC_4_TUNNEL:
      vals = {-7.70146936482581, -3.21868369245987, -5.44023748418735, 15.6051950593180, 12.7611835388515, -4.46952393442309, -11.7240576326183, -9.23038948829007};
      iso = -1.7660;
      break;
    case AmbiguousCase::MC_7_TUNNEL:
      vals = {-3.42744283804455, 0.621278122151001, 4.48110777981235, -1.95551129669134, 2.30448107596369, -1.04182240925489, -3.51087814405650, -6.44976786808517};
      iso = -0.6;
      break;
    case AmbiguousCase::MC_10_TUNNEL:
      vals = {-0.100000000000000, -6.11000000000000, 2, 10.2000000000000, 10.8000000000000, 1.80000000000000, -8.20000000000000, -0.180000000000000};
      iso = 1.10918;
      break;
    case AmbiguousCase::MC_12_TUNNEL:
      vals = {-3.37811990337124, 0.473258332744286, 2.54344310345736, 7.87658724379480, 4.38700713005133, -1.49950251870885, -4.21025867362045, -1.00233824192217};
      iso = 0.0708;
      break;
    case AmbiguousCase::MC_13_TUNNEL:
      vals = {2.74742516087490, -3.39187542578189, -12.5297639669456, 0.431517989649243, -6.92460546400188, 2.52228314017858, 14.6950568276448, -10.0732624062474};
      iso = -1.3064;
      break;
    case AmbiguousCase::MC_13_12_VTS_CONFIG_A:
      vals = {0.546912886195662, -0.421103532406922, -0.643375084081520, 0.855507421818445, -0.260686312588506, 0.206413666735986, 0.237274227130530, -0.183297728364877};
      iso = 0.0293;
      break;
    case AmbiguousCase::MC_13_12_VTS_CONFIG_B:
      vals = {1069, 843, 950, 1133, 958, 1029, 1198, 946};
      iso = 1007.4;
      break;
    case AmbiguousCase::MC_13_SIMPLE:
      vals = {0.520482995461163, -0.839814387388296, -0.467491517013617, 0.937814095887345, -0.825777099007084, 0.506695544835103, 0.345318915961394, -0.861107217966913};
      iso = 0.0293;
      break;
    default:
      assert(false);
  }

  FT min = std::numeric_limits<FT>::max();
  FT max = std::numeric_limits<FT>::lowest();
  for (const FT v : vals)
  {
    min = std::min(min, v);
    max = std::max(max, v);
  }
  max += 1e-3;

  for(std::size_t x=0; x<grid.xdim(); ++x)
    for(std::size_t y=0; y<grid.ydim(); ++y)
      for(std::size_t z=0; z<grid.zdim(); ++z)
        values(x, y, z) = rand.uniform_real<FT>(max, 2 * max - min);

  values(1, 1, 1) = vals[0];
  values(1, 1, 2) = vals[1];
  values(1, 2, 1) = vals[2];
  values(1, 2, 2) = vals[3];
  values(2, 1, 1) = vals[4];
  values(2, 1, 2) = vals[5];
  values(2, 2, 1) = vals[6];
  values(2, 2, 2) = vals[7];

  assert(check_closed_not_empty(grid, values, iso));

  return iso;
}

enum class SingularCase {
  CASE_1,
  CASE_2,
  CASE_3,
  CASE_4
};

template <typename Grid>
typename Grid::Geom_traits::FT generate_predefined_singular(const Grid& grid, IS::Interpolated_discrete_values_3<Grid>& values, const SingularCase topology_case) {
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  CGAL::Random rand;

  assert(grid.xdim() == 4);
  assert(grid.ydim() == 4);
  assert(grid.zdim() == 4);

  std::array<FT, 8> vals;
  FT iso;
  bool invert;

  switch (topology_case) {
    case SingularCase::CASE_1:
      vals = {899.00000000000000, 1339.0000000000000, 1210.0000000000000, 1220.0000000000000, 754.00000000000000, 998.00000000000000, 928.00000000000000, 1264.0000000000000};
      iso = 1233.6;
      invert = false;
      break;
    case SingularCase::CASE_2:
      vals = {1003.0000000000000, 1523.0000000000000, 1230.0000000000000, 1246.0000000000000, 1062.0000000000000, 1448.0000000000000, 1293.0000000000000, 1029.0000000000000};
      iso = 1233.6;
      invert = true;
      break;
    case SingularCase::CASE_3:
      vals = {1311.0000000000000, 1225.0000000000000, 1255.0000000000000, 1256.0000000000000, 1203.0000000000000, 1237.0000000000000, 1118.0000000000000, 1288.0000000000000};
      iso = 1233.6;
      invert = true;
      break;
    case SingularCase::CASE_4:
      vals = {1763.0000000000000, 1052.0000000000000, 1815.0000000000000, 1050.0000000000000, 960.00000000000000, 1325.0000000000000, 1150.0000000000000, 1260.0000000000000};
      iso = 1233.6;
      invert = false;
      break;
    default:
      assert(false);
  }

  FT min = std::numeric_limits<FT>::max();
  FT max = std::numeric_limits<FT>::lowest();
  for (FT& v : vals)
  {
    if (invert)
      v = 2 * iso - v;

    min = std::min(min, v);
    max = std::max(max, v);
  }
  max += 1e-3;

  for(std::size_t x=0; x<grid.xdim(); ++x)
    for(std::size_t y=0; y<grid.ydim(); ++y)
      for(std::size_t z=0; z<grid.zdim(); ++z)
        values(x, y, z) = rand.uniform_real<FT>(max, 2 * max - min);

  values(1, 1, 1) = vals[0];
  values(1, 1, 2) = vals[1];
  values(2, 1, 1) = vals[2];
  values(2, 1, 2) = vals[3];
  values(1, 2, 1) = vals[4];
  values(1, 2, 2) = vals[5];
  values(2, 2, 1) = vals[6];
  values(2, 2, 2) = vals[7];

  assert(check_closed_not_empty(grid, values, iso));

  return iso;
}



template <typename Grid>
void generate_random_easy(const Grid& grid, IS::Interpolated_discrete_values_3<Grid>& values)
{
  using Geom_traits = typename Grid::Geom_traits;
  using FT = typename Geom_traits::FT;

  while (true) {
    generate_random(grid, values);

    if (!check_closed_not_empty(grid, values))
      continue;

    if (check_iso_vertices(grid, values))
      continue;

    if (check_equal_edges(grid, values))
      continue;

    bool easy = true;

    for(std::size_t x=0; x<grid.xdim(); ++x) {
      for(std::size_t y=0; y<grid.ydim(); ++y) {
        for(std::size_t z=0; z<grid.zdim(); ++z)
        {
          // singular face diagonals
          if (x < grid.xdim() - 1 && y < grid.ydim() - 1 && std::abs(values(x, y, z) - values(x + 1, y + 1, z)) < 1e-6)
            easy = false;
          if (x < grid.xdim() - 1 && z < grid.zdim() - 1 && std::abs(values(x, y, z) - values(x + 1, y, z + 1)) < 1e-6)
            easy = false;
          if (y < grid.ydim() - 1 && z < grid.zdim() - 1 && std::abs(values(x, y, z) - values(x, y + 1, z + 1)) < 1e-6)
            easy = false;

          // singular cell diagonal
          if (x < grid.xdim() - 1 && y < grid.ydim() - 1 && z < grid.zdim() - 1 && std::abs(values(x, y, z) - values(x + 1, y + 1, z + 1)) < 1e-6)
            easy = false;

          // ambiguous faces
          if (x < grid.xdim() - 1 && y < grid.ydim() - 1 && values(x, y, z) * values(x + 1, y + 1, z) > 0 && values(x + 1, y, z) * values(x, y + 1, z) > 0 && values(x, y, z) * values(x + 1, y, z) < 0)
            easy = false;
          if (x < grid.xdim() - 1 && z < grid.zdim() - 1 && values(x, y, z) * values(x + 1, y, z + 1) > 0 && values(x + 1, y, z) * values(x, y, z + 1) > 0 && values(x, y, z) * values(x + 1, y, z) < 0)
            easy = false;
          if (y < grid.ydim() - 1 && z < grid.zdim() - 1 && values(x, y, z) * values(x, y + 1, z + 1) > 0 && values(x, y + 1, z) * values(x, y, z + 1) > 0 && values(x, y, z) * values(x, y + 1, z) < 0)
            easy = false;

          // interior ambiguity
          if (x < grid.xdim() - 1 && y < grid.ydim() - 1 && z < grid.zdim() - 1 && values(x, y, z) * values(x + 1, y + 1, z + 1) > 0
              && values(x, y, z) * values(x + 1, y, z) < 0 && values(x, y, z) * values(x, y + 1, z) < 0 && values(x, y, z) * values(x, y, z + 1) < 0
              && values(x + 1, y + 1, z + 1) * values(x + 1, y + 1, z) < 0 && values(x + 1, y + 1, z + 1) * values(x + 1, y, z + 1) < 0 && values(x + 1, y + 1, z + 1) * values(x, y + 1, z + 1) < 0)
            easy = false;
        }
      }
    }

    if (easy)
      break;
  }
}

void test_random()
{
  using K = CGAL::Simple_cartesian<double>;
  using Grid = IS::Cartesian_grid_3<K>;
  using Values = IS::Interpolated_discrete_values_3<Grid>;
  using Domain = IS::Marching_cubes_domain_3<Grid, Values>;
  using FT = typename K::FT;
  using Point = typename K::Point_3;
  using Mesh = CGAL::Surface_mesh<Point>;
  using Point_range = std::vector<Point>;
  using Triangle_range = std::vector<std::vector<std::size_t> >;

  Grid grid { Point{-1., -1., -1.}, Point{1., 1., 1.}, std::array<std::size_t, 3>{4, 4, 4} };
  Values values { grid };
  Domain domain { grid, values };

  Mesh debug_grid;
  auto debug_grid_creator = [&](const Domain::cell_descriptor& c)
    {
      std::vector<Mesh::Vertex_index> cell_vertices;
      for (const auto& v : domain.cell_vertices(c)) {
        cell_vertices.push_back(debug_grid.add_vertex(domain.point(v)));
      }
      debug_grid.add_face(cell_vertices[6], cell_vertices[2], cell_vertices[0], cell_vertices[4]);
      debug_grid.add_face(cell_vertices[1], cell_vertices[3], cell_vertices[7], cell_vertices[5]);
      debug_grid.add_face(cell_vertices[0], cell_vertices[1], cell_vertices[5], cell_vertices[4]);
      debug_grid.add_face(cell_vertices[6], cell_vertices[7], cell_vertices[3], cell_vertices[2]);
      debug_grid.add_face(cell_vertices[2], cell_vertices[3], cell_vertices[1], cell_vertices[0]);
      debug_grid.add_face(cell_vertices[4], cell_vertices[5], cell_vertices[7], cell_vertices[6]);
    };
    domain.template for_each_cell<CGAL::Sequential_tag>(debug_grid_creator);
  CGAL::IO::write_OFF("debug_grid.off", debug_grid);

  FT iso = 0;
  // generate_random_easy(grid, values);
  // iso = generate_predefined_inner_ambiguity(grid, values, AmbiguousCase::MC_13_TUNNEL);
  iso = generate_predefined_singular(grid, values, SingularCase::CASE_3);

  {
    Grid grid_high_res { Point{-1., -1., -1.}, Point{1., 1., 1.}, std::array<std::size_t, 3>{151, 151, 151} };
    IS::Value_function_3<Grid> values_high_res { values, grid_high_res };


    Grid grid_high_res_large { Point{-50., -50., -50.}, Point{50., 50., 50.}, std::array<std::size_t, 3>{301, 301, 301} };

    struct Extra {
      Grid& g;
      Values& values;

      FT operator()(const typename IS::partition_traits<Grid>::vertex_descriptor& v) {
        return this->operator()(IS::partition_traits<Grid>::point(v, g));
      }

      FT operator()(const K::Point_3& p) {
        typename K::Compute_x_3 x_coord = g.geom_traits().compute_x_3_object();
        typename K::Compute_y_3 y_coord = g.geom_traits().compute_y_3_object();
        typename K::Compute_z_3 z_coord = g.geom_traits().compute_z_3_object();
        typename K::Construct_vertex_3 vertex = g.geom_traits().construct_vertex_3_object();

        // trilinear interpolation of stored values
        const K::Iso_cuboid_3& span = g.span();
        const K::Vector_3& spacing = g.spacing();

        // TODO: only works for cartesian grid
        // calculate min index including border case
        const K::Point_3& min_p = vertex(span, 0);
        std::size_t i = 1;
        std::size_t j = 1;
        std::size_t k = 1;

        // calculate coordinates of min index
        const FT min_x = x_coord(min_p) + i * spacing[0];  // TODO: x_coord ...
        const FT min_y = y_coord(min_p) + j * spacing[1];
        const FT min_z = z_coord(min_p) + k * spacing[2];

        // interpolation factors between 0 and 1
        const FT f_i = (x_coord(p) - min_x) / spacing[0];
        const FT f_j = (y_coord(p) - min_y) / spacing[1];
        const FT f_k = (z_coord(p) - min_z) / spacing[2];

        // read the value at all 8 corner points
        const FT g000 = values(typename IS::partition_traits<Grid>::vertex_descriptor{i + 0, j + 0, k + 0});
        const FT g001 = values(typename IS::partition_traits<Grid>::vertex_descriptor{i + 0, j + 0, k + 1});
        const FT g010 = values(typename IS::partition_traits<Grid>::vertex_descriptor{i + 0, j + 1, k + 0});
        const FT g011 = values(typename IS::partition_traits<Grid>::vertex_descriptor{i + 0, j + 1, k + 1});
        const FT g100 = values(typename IS::partition_traits<Grid>::vertex_descriptor{i + 1, j + 0, k + 0});
        const FT g101 = values(typename IS::partition_traits<Grid>::vertex_descriptor{i + 1, j + 0, k + 1});
        const FT g110 = values(typename IS::partition_traits<Grid>::vertex_descriptor{i + 1, j + 1, k + 0});
        const FT g111 = values(typename IS::partition_traits<Grid>::vertex_descriptor{i + 1, j + 1, k + 1});

        // interpolate along all axes by weighting the corner points
        const FT lambda000 = (FT(1) - f_i) * (FT(1) - f_j) * (FT(1) - f_k);
        const FT lambda001 = (FT(1) - f_i) * (FT(1) - f_j) * f_k;
        const FT lambda010 = (FT(1) - f_i) * f_j * (FT(1) - f_k);
        const FT lambda011 = (FT(1) - f_i) * f_j * f_k;
        const FT lambda100 = f_i * (FT(1) - f_j) * (FT(1) - f_k);
        const FT lambda101 = f_i * (FT(1) - f_j) * f_k;
        const FT lambda110 = f_i * f_j * (FT(1) - f_k);
        const FT lambda111 = f_i * f_j * f_k;

        // add weighted corners
        return g000 * lambda000 + g001 * lambda001 +
              g010 * lambda010 + g011 * lambda011 +
              g100 * lambda100 + g101 * lambda101 +
              g110 * lambda110 + g111 * lambda111;
      }
    };

    IS::Value_function_3<Grid> values_high_res_large {Extra{grid, values}, grid_high_res };
    IS::Marching_cubes_domain_3<Grid, IS::Value_function_3<Grid>> domain_high_res { grid_high_res, values_high_res };

    Point_range points_high_res;
    Triangle_range triangles_high_res;
    IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain_high_res, iso, points_high_res, triangles_high_res, CGAL::parameters::use_topologically_correct_marching_cubes(false));

    #define CGAL_TESTUISTE_ISOSURFACING_OUTPUT
    #ifdef CGAL_TESTUISTE_ISOSURFACING_OUTPUT
    CGAL::IO::write_polygon_soup("test_random_high_res.off", points_high_res, triangles_high_res);
    #endif
  }

  {
    Point_range points_mc;
    Triangle_range triangles_mc;
    IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain, iso, points_mc, triangles_mc, CGAL::parameters::use_topologically_correct_marching_cubes(false));

    #define CGAL_TESTUISTE_ISOSURFACING_OUTPUT
    #ifdef CGAL_TESTUISTE_ISOSURFACING_OUTPUT
    CGAL::IO::write_polygon_soup("test_random_mc.off", points_mc, triangles_mc);
    #endif
  }

  Point_range points;
  Triangle_range triangles;
  IS::marching_cubes<CGAL::Sequential_tag>(domain, iso, points, triangles, CGAL::parameters::use_topologically_correct_marching_cubes(true));

  #define CGAL_TESTUISTE_ISOSURFACING_OUTPUT
  #ifdef CGAL_TESTUISTE_ISOSURFACING_OUTPUT
  CGAL::IO::write_polygon_soup("test_random.off", points, triangles);
  #endif
  
  assert(points.size() && triangles.size());
  assert(!has_duplicate_points(points, triangles));
  assert(!has_duplicate_polygons(points, triangles));
  assert(!has_isolated_vertices(points, triangles));

  assert(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles));
  Mesh m;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, triangles, m);

  assert(!has_degenerate_faces(m));
  assert(is_manifold(m));
  assert(is_watertight(m));
}

#endif // CGAL_ISOSURFACING_TEST_TOPOLOGY_H
