// g++ src/test_cosmo_distance.cpp -lcosmo -o test_cosmo_distance

#include <iostream>
#include <iomanip>
#include "cosmo/cosmo.h"

int main(int argc, char **argv) {

    double OmegaLambda(0.73), OmegaMatter(0);

    // set up cosmology
    cosmo::AbsHomogeneousUniversePtr cosmology;
    if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
    try {
        cosmology.reset(new cosmo::LambdaCdmUniverse(OmegaLambda, OmegaMatter));
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    double zmin(0.0), zmax(4.0), nz(101);
    double dz = (zmax - zmin) / (nz - 1);
    for(int i = 0; i < nz; ++i) {
        double z = i * dz;
        double distance = cosmology->getLineOfSightComovingDistance(z);
        std::cout << std::fixed << std::setprecision(8) << z << " "
                  << std::fixed << std::setprecision(8) << distance << std::endl;
    }

    return 0;
}