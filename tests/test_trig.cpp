#include <iostream>
#include <iomanip>
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

void trig(double los1, double los2, double r1, double r2) {

}


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

    double ra1, dec1, z1, r1;
    double ra2, dec2, z2, r2;

    double eps;

    double cos_theta, r1_sq, two_r1_cos_theta, r_sq, dr_para, dr_perp;


    // - ra2=ra1, dec2=dec1+eps, z1=z2: dr_par=0, dr_perp = dist * eps
    ra1 = to_radians(100);
    dec1 = to_radians(15);
    z1 = 2.33;

    eps = to_radians(1.0);

    ra2 = ra1;
    dec2 = dec1 + eps;
    z2 = z1;

    r1 = cosmology->getLineOfSightComovingDistance(z1);
    r2 = cosmology->getLineOfSightComovingDistance(z2);

    SkyObject los1(ra1, dec1);
    SkyObject los2(ra2, dec2);

    cos_theta = los1.getCosSeparation(los2);
    r1_sq = r1 * r1;
    two_r1_cos_theta = 2 * r1 * cos_theta;
    r_sq = r1_sq + (r2 - two_r1_cos_theta) * r2;
    dr_para = std::fabs(r1 - r2);
    dr_perp = std::sqrt(r_sq - dr_para*dr_para);

    std::cout << 0 << " " << r1 * eps << std::endl;
    std::cout << dr_para << " " << dr_perp << std::endl;
    std::cout << std::endl;

	// - ra2=ra1+eps, dec2=dec1, z1=z2: dr_par=0, dr_perp = dist * eps * cos(dec)
	ra2 = ra1 + eps;
	dec2 = dec1;

    los2 = SkyObject(ra2, dec2);

    cos_theta = los1.getCosSeparation(los2);
    r1_sq = r1 * r1;
    two_r1_cos_theta = 2 * r1 * cos_theta;
    r_sq = r1_sq + (r2 - two_r1_cos_theta) * r2;
    dr_para = std::fabs(r1 - r2);
    dr_perp = std::sqrt(r_sq - dr_para*dr_para);

    std::cout << 0 << " " << r1 * eps * std::cos(dec1) << std::endl;
    std::cout << dr_para << " " << dr_perp << std::endl;
    std::cout << std::endl;

	// - ra2=ra1, dec2=dec1, z2=z1+eps: dr_par=eps * d(dist)/dz, dr_perp = 0

	ra2 = ra1;
	dec2 = dec1;
	z2 = z1 + eps;
	r2 = cosmology->getLineOfSightComovingDistance(z2);

    los2 = SkyObject(ra2, dec2);

    cos_theta = los1.getCosSeparation(los2);

    r1_sq = r1 * r1;
    two_r1_cos_theta = 2 * r1 * cos_theta;
    r_sq = r1_sq + (r2 - two_r1_cos_theta) * r2;
    dr_para = std::fabs(r1 - r2);
    dr_perp = std::sqrt(std::fabs(r_sq - dr_para*dr_para));

    std::cout << eps * (r2 - r1) / (z2 - z1) << " " << 0 << std::endl;
    std::cout << dr_para << " " << dr_perp << std::endl;
    std::cout << std::endl;

    return 0;
}
