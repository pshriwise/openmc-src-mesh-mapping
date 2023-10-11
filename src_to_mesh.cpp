
#include <iostream>

#include "openmc/capi.h"

int main(int argc, char** argv) {

  std::cout << "Hello world" << std::endl;

  openmc_init(argc, argv, nullptr);

  

  return 0;

}
