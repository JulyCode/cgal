#ifndef CGAL_ISOSURFACING_TEST_UTIL_H
#define CGAL_ISOSURFACING_TEST_UTIL_H

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <cassert>
#include <iostream>

template <typename PointRange, typename PolygonRange>
bool has_duplicate_points(PointRange points, PolygonRange polygons) // intentional copies
{
  return CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(points, polygons) != 0;
}

template <typename PointRange, typename PolygonRange>
bool has_duplicate_polygons(PointRange points, PolygonRange polygons)
{
  return CGAL::Polygon_mesh_processing::merge_duplicate_polygons_in_polygon_soup(points, polygons) != 0;
}

template <typename PointRange, typename PolygonRange>
bool has_isolated_vertices(PointRange points, PolygonRange polygons)
{
  return CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(points, polygons) != 0;
}

template <typename PolygonMesh>
bool is_manifold(PolygonMesh& pmesh)
{
  return CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(pmesh, CGAL::parameters::dry_run(true)) == 0;
}

template <typename PolygonMesh>
bool has_degenerate_faces(PolygonMesh& pmesh)
{
  std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor> dfs;
  CGAL::Polygon_mesh_processing::degenerate_faces(pmesh, std::inserter(dfs, dfs.begin()));
  return !dfs.empty();
}

// template <typename PointRange, typename PolygonRange>
// bool is_closed(PointRange points, PolygonRange polygons) {
//   struct Edge {
//     std::size_t v0, v1;
//   };
//   struct EdgeHash {
//     std::size_t operator()(const Edge& e) const {
//       std::size_t hash = 17;
//       hash = hash * 31 + e.v0;
//       hash = hash * 31 + e.v1;
//       return hash;
//     }
//   };
//   struct EdgeEqual {
//     std::size_t operator()(const Edge& e1, const Edge& e2) const {
//       return e1.v0 == e2.v0 && e1.v1 == e2.v1;
//     }
//   };

//   std::unordered_map<Edge, std::size_t, EdgeHash, EdgeEqual> edges;

//   for (const auto& face : polygons) {
//     for (std::size_t i = 0; i < face.size(); i++) {
//       const std::size_t v0 = face[i];
//       const std::size_t v1 = face[(i + 1) % face.size()];

//       Edge edge = {v0, v1};
//       if (v0 > v1) {
//         edge = {v1, v0};
//       }

//       if (edges.find(edge) == edges.end()) {
//         edges[edge] = 1;
//       } else {
//         edges[edge]++;
//       }
//     }
//   }

//   for (auto& [edge, occurrences] : edges) {
//     if (occurrences != 2) {
//       return false;
//     }
//   }
//   return true;
// }

// template <typename PolygonMesh>
// bool is_closed(const PolygonMesh& pmesh) {
//   for (const auto& edge : pmesh.edges()) {
//     if (pmesh.is_border(edge)) {
//       return false;
//     }
//   }
//   return true;
// }

template <typename PolygonMesh>
bool is_self_intersecting(const PolygonMesh& pmesh) {
  return CGAL::Polygon_mesh_processing::does_self_intersect(pmesh);
}

template <typename PolygonMesh>
bool is_watertight(PolygonMesh& pmesh) {
  const bool manifold = is_manifold(pmesh);
  const bool closed = is_closed(pmesh);
  const bool self_intersecting = is_self_intersecting(pmesh);
  return manifold && closed && !self_intersecting;
}


// template <typename PolygonMesh>
// bool check_mesh_distance(const PolygonMesh& m0, const PolygonMesh& m1)
// {
//   auto dist = CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance<CGAL::Sequential_tag>(
//     m0, m1, CGAL::parameters::number_of_points_per_area_unit(4000));
//   std::cout << dist << std::endl;
//   return true;
// }

#endif // CGAL_ISOSURFACING_TEST_UTIL_H
