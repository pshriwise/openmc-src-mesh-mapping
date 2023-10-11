
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "openmc/capi.h"
#include "openmc/mesh.h"
#include "openmc/source.h"

int main(int argc, char** argv) {

  std::cout << "Hello world" << std::endl;

  openmc_init(0, nullptr, nullptr);

  int total_sites = 100000;
  int max_sites_incr = 100;
  int sites_sampled {};

  // accessing mesh
  const auto& mesh = openmc::model::meshes[openmc::model::mesh_map[1]];

  std::cout << "Mesh ID: " << mesh->id() << std::endl;
  std::cout << "Mesh bins: " << mesh->n_bins() << std::endl;

  std::vector<size_t> bin_counts(mesh->n_bins(), 0);

  while (sites_sampled < total_sites) {

    int n_sites = std::min(total_sites - sites_sampled, max_sites_incr);

    // sample source sites
    std::vector<openmc::SourceSite> sites(n_sites);

    uint64_t seed = 1;

    for (auto& s : sites) {
      s = openmc::sample_external_source(&seed);
    }

    sites_sampled += n_sites;

    for (const auto& s : sites) {
      int bin = mesh->get_bin(s.r);
      bin_counts[bin] += 1;
    }
  }

  std::vector<double> rel_strengths(bin_counts.begin(), bin_counts.end());
  for (auto& rs : rel_strengths) { rs /= sites_sampled; }

  // Determine total source strength
  double total_strength = 0.0;
  for (auto& s : openmc::model::external_sources)
    total_strength += s->strength();

  std::cout << "Total external source strength: " << total_strength << std::endl;

  // TODO: extend to all mesh types
  openmc::UnstructuredMesh* umesh_ptr = dynamic_cast<openmc::UnstructuredMesh*>(mesh.get());
  if (!umesh_ptr) { std::cerr << "Non-unstructured mesh used" << std::endl; return 1; }

  // generate a text output as CSV
  // Each row contains: OpenMC bin #, centroid (xyz), volume, rel. src strength
  std::ofstream output("mesh_src_strengths.csv");
  output << "bin, cx, cy, cz, volume, rel. src" << std::endl;

  for (int i = 0; i < mesh->n_bins(); i++) {
    std::stringstream line {};
    line << i << ", ";
    auto centroid = umesh_ptr->centroid(i);
    line << centroid.x << ", " << centroid.y << ", " << centroid.z << ", ";
    line << mesh->volume(i) << ", ";
    line << rel_strengths[i] << std::endl;

    output << line.str();

  }

  output.close();

  openmc_finalize();

  return 0;

}
