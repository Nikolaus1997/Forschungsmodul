
#include "computation/computation.h"
#include <iostream>
#include <cstdlib>
#ifdef DGFM_USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }
#ifdef DGFM_USE_MPI
  MPI_Init(&argc,&argv);
  int r=0,p=1; MPI_Comm_rank(MPI_COMM_WORLD,&r); MPI_Comm_size(MPI_COMM_WORLD,&p);
  if(r==0) std::printf("MPI up with %d ranks\n", p);
#endif


// std::cout <<R"(====================================================================================)"<<std::endl;
// std::cout <<R"(________    ________  .____    ________  ___________ ______________________________ )"<<std::endl; 
// std::cout <<R"(\______ \  /  _____/  |    |   \_____  \ \_   _____//   _____/\_   _____/\______   \)"<<std::endl;
// std::cout <<R"( |    |  \/   \  ___  |    |    /   |   \ |    __)_ \_____  \  |    __)_  |       _/)"<<std::endl;
// std::cout <<R"( |    `   \    \_\  \ |    |___/    |    \|        \/        \ |        \ |    |   \)"<<std::endl;
// std::cout <<R"(/_______  /\______  / |_______ \_______  /_______  /_______  //_______  / |____|_  /)"<<std::endl;
// std::cout <<R"(\/        \/          \/       \/        \/        \/         \/         \/         )"<<std::endl;
// std::cout <<R"(====================================================================================)"<<std::endl;

std::cout <<R"(============================================================================================)"<<std::endl;
std::cout <<R"( ________  ________          ___       ________  _______   ________  _______   ________     )"<<std::endl; 
std::cout <<R"(|\   ___ \|\   ____\        |\  \     |\   __  \|\  ___ \ |\   ____\|\  ___ \ |\   __  \    )"<<std::endl; 
std::cout <<R"(\ \  \_|\ \ \  \___|        \ \  \    \ \  \|\  \ \   __/ \ \  \___|\ \   __/ \ \  \|\  \   )"<<std::endl; 
std::cout <<R"( \ \  \ \\ \ \  \  ___       \ \  \    \ \  \\\  \ \  \   _\ \_____  \ \  \    \ \   _  _\   )"<<std::endl; 
std::cout <<R"(  \ \  \_\\ \ \  \|\  \       \ \  \____\ \  \\\  \ \  \_|\ \|____|\  \ \  \_\  \ \  \\  \| )"<<std::endl; 
std::cout <<R"(   \ \_______\ \_______\       \ \_______\ \_______\ \_______\____\_\  \ \_______\ \__\\ _\ )"<<std::endl; 
std::cout <<R"(    \|_______|\|_______|        \|_______|\|_______|\|_______|\_________\|_______|\|__|\|__|)"<<std::endl;
std::cout <<R"(                                                             \|_________|                   )"<<std::endl;
std::cout <<R"(============================================================================================)"<<std::endl;
                                                                                            
                                                                                            



  // read in the first argument

  std::string filename = argv[1];
  
  auto computation = Computation();
  computation.initialize(filename);
  computation.runSimulation();

  #ifdef DGFM_USE_MPI
  MPI_Finalize();
  #endif
  return EXIT_SUCCESS;
}