
#include <iostream>

#include "openmc/capi.h"
#include "openmc/mesh.h"

int main(int argc, char** argv) {

  std::cout << "Hello world" << std::endl;

  openmc_init(argc, argv, nullptr);

  const auto& mesh = openmc::model::meshes[openmc::model::mesh_map[1]];

  std::cout << "Mesh ID: " << mesh->id() << std::endl;
  std::cout << "Mesh bins: " << mesh->n_bins() << std::endl;

  openmc_finalize();

  return 0;

}
