
#include <iostream>
#include <vector>

#include "openmc/capi.h"
#include "openmc/mesh.h"
#include "openmc/source.h"

int main(int argc, char** argv) {

  std::cout << "Hello world" << std::endl;

  openmc_init(argc, argv, nullptr);

  int num_sites = 10;

  // sample source sites
  std::vector<openmc::SourceSite> sites(num_sites);

  uint64_t seed = 1;

  for (int i = 0; i < num_sites; i++) {
    sites[i] = openmc::sample_external_source(&seed);
    std::cout << "Source site xyz: " << sites[i].r.x << ", " << sites[i].r.y << ", " << sites[i].r.z << std::endl;
  }

  const auto& mesh = openmc::model::meshes[openmc::model::mesh_map[1]];

  std::cout << "Mesh ID: " << mesh->id() << std::endl;
  std::cout << "Mesh bins: " << mesh->n_bins() << std::endl;

  openmc_finalize();

  return 0;

}
