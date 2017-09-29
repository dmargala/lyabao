#include <iostream>
#include <iomanip>

#include <fstream>
#include <sstream>
#include <string>

#include <cmath>
#include "cosmo/cosmo.h"

const double pi = std::atan(1)*4;

double to_radians(double angle) {
    return angle * pi / 180.0;
}

/// A simple representation of an object in the sky.
struct SkyObject {
    double ra, dec, sin_dec, cos_dec, sin_ra, cos_ra;
    /// Create a new SkyObject
    /// @param ra_ Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
    /// @param dec_ Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
    SkyObject(double ra_, double dec_) :
    ra(ra_), dec(dec_) {
        sin_dec = std::sin(dec);
        cos_dec = std::cos(dec);
        sin_ra = std::sin(ra);
        cos_ra = std::cos(ra);
    }
    /// Copy constructor
    /// @param other The other SkyObject
    SkyObject(const SkyObject& other) {
        ra = other.ra; dec = other.dec;
        sin_dec = other.sin_dec; cos_dec = other.cos_dec;
        sin_ra = other.sin_ra; cos_ra = other.cos_ra;
    };
    SkyObject() {};
    /// Return cosine of the angular separation between this and the provided SkyObject
    /// @param other The other SkyObject
    double getCosSeparation(SkyObject const &other) const {
        return sin_dec*other.sin_dec + cos_dec*other.cos_dec*(sin_ra*other.sin_ra + cos_ra*other.cos_ra);
    }
};

int main(int argc, char **argv) {

    // set up cosmology
    double OmegaLambda(0.73), OmegaMatter(0);
    cosmo::AbsHomogeneousUniversePtr cosmology;
    if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
    try {
        cosmology.reset(new cosmo::LambdaCdmUniverse(OmegaLambda, OmegaMatter));
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }


    std::ifstream infile("test_trig3.txt");

    unsigned int num_lines(0);

    std::string line;
    for(std::string line; std::getline(infile, line); ) {

        if (line[0] == '#') continue;

        std::istringstream iss(line);

        double ra1, dec1, z1, r1, ra2, dec2, z2, r2,
            test_dr_para, test_dr_perp, delta_dr_para, delta_dr_perp;

        if (!(iss >> ra1 >> dec1 >> z1 >> r1 >> ra2 >> dec2 >> z2 >> r2
                  >> test_dr_para >> test_dr_perp >> delta_dr_para >> delta_dr_perp)) { break; } // error

        num_lines++;

        // double r1, r2;
        // r1 = cosmology->getLineOfSightComovingDistance(z1);
        // r2 = cosmology->getLineOfSightComovingDistance(z2);

        SkyObject los1(ra1, dec1);
        SkyObject los2(ra2, dec2);

        float cos_theta, r1_sq, two_r1_cos_theta, r_sq, dr_para, dr_perp;

        cos_theta = los1.getCosSeparation(los2);
        r1_sq = r1 * r1;
        two_r1_cos_theta = 2 * r1 * cos_theta;
        r_sq = r1_sq + (r2 - two_r1_cos_theta) * r2;

        dr_para = std::fabs(r1 - r2);
        dr_perp = std::sqrt(std::fabs(r_sq - dr_para*dr_para));

        std::cout << std::fixed << std::setprecision(8) << dr_para - test_dr_para << " "
                  << std::fixed << std::setprecision(8) << dr_perp - test_dr_perp << " "
                  << std::fixed << std::setprecision(8) << delta_dr_para << " "
                  << std::fixed << std::setprecision(8) << delta_dr_perp << std::endl;
    }

    std::cout << "num_lines: " << num_lines << std::endl;

    return 0;
}
