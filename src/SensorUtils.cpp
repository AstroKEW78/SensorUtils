#include "SensorUtils.h"

#include <cmath>> 

#include "float.h"

#include <vector>

#include <armadillo>

using namespace std;


 /**
   * Computes and returns phase angle, in radians, given the positions of the
   * observer and illuminator.
   *
   * Phase Angle: The angle between the vector from the surface intersection point to
   * the observer (usually the spacecraft) and the vector from the surface intersection
   * point to the illuminator (usually the sun).
   *
   * @author Kaj Williams
   *
   * @param observerBodyFixedPosition  Three dimensional position of the observer,
   *                                   in the coordinate system of the target body.
   * @param illuminatorBodyFixedPosition Three dimensional position for the illuminator,
   *                                     in the body-fixed coordinate system.
   * @param surfaceIntersection Three dimensional position for the ground (surface intersection) point,
   *                                     in the body-fixed coordinate system.
   * @return @b double Phase angle, in radians.
   */

double PhaseAngle(const std::vector<double> &observerBodyFixedPosition,
                                const std::vector<double> &illuminatorBodyFixedPosition,
				const std::vector<double> &surfaceIntersection) {

    //convert the vector to an arma::vec
    arma::vec observer = arma::zeros<arma::vec>(3);
    observer = arma::conv_to<arma::vec>::from(observerBodyFixedPosition);

    //convert the illuminatorBodyFixedPosition vector to an arma::vec
    arma::vec illuminator = arma::zeros<arma::vec>(3);
    illuminator = arma::conv_to<arma::vec>::from(illuminatorBodyFixedPosition);

    //convert the surfaceIntersection vector to an arma::vec
    arma::vec surface = arma::zeros<arma::vec>(3);
    surface = arma::conv_to<arma::vec>::from(surfaceIntersection);    

    // Get vector from surface point to observer and normalise it 
    arma::vec surfaceToObserver = arma::zeros<arma::vec>(3);
    arma::vec normSurfaceToObserver = arma::zeros<arma::vec>(3);
    surfaceToObserver = observer - surface; 
    normSurfaceToObserver = arma::normalise(surfaceToObserver);

    // Get vector from surface point to sun and normalise it
    arma::vec surfaceToSun = arma::zeros<arma::vec>(3);
    arma::vec normSurfaceToSun = arma::zeros<arma::vec>(3);
    surfaceToSun = illuminator - surface; 
    normSurfaceToSun = arma::normalise(surfaceToSun);

    double cos_angle=arma::dot(normSurfaceToObserver,normSurfaceToSun);
   
    if(cos_angle >= 1.0) return 0.0; 
    if(cos_angle <= -1.0) return M_PI;
	
    return acos(cos_angle);
    

}


/**
   * Computes and returns local incidence angle.
   *
   * @author Kaj Williams
   *
   * @param normal  local normal vector.
   *                                   
   * @param surfSpaceVect Three dimensional normalized surface spacecraft vector
   * @param surfaceIntersection Three dimensional position for the ground (surface intersection) point,
   *                                     in the body-fixed coordinate system.
   * @return @b double local incidence angle, in radians.
   */

double LocalIncidenceAngle(const std::vector<double> &normal,
                                const std::vector<double> &,
        const std::vector<double> &surfaceIntersection) {


    //convert the local normal vector to an arma::vec
    arma::vec localNormal = arma::zeros<arma::vec>(3);
    localNormal = arma::conv_to<arma::vec>::from(normal);


    // Check to make sure normal is valid
    SpiceDouble mag;
    unorm_c(normal,normal,&mag);
    if (mag == 0.) {
      success = false;
      return;
    }

    // get a normalized surface spacecraft vector
    //convert the local normal vector to an arma::vec
    arma::vec surfSpaceVector = arma::zeros<arma::vec>(3);
    surfSpaceVector = arma::conv_to<arma::vec>::from(surfSpaceVect);

    SpiceDouble surfSpaceVect[3], unitizedSurfSpaceVect[3], dist;
    std::vector<double> sB = bodyRotation()->ReferenceVector(
        instrumentPosition()->Coordinate());

    SpiceDouble pB[3];
    SurfacePoint surfacePoint = GetSurfacePoint();
    pB[0] = surfacePoint.GetX().kilometers();
    pB[1] = surfacePoint.GetY().kilometers();
    pB[2] = surfacePoint.GetZ().kilometers();

    vsub_c((SpiceDouble *) &sB[0], pB, surfSpaceVect);
    unorm_c(surfSpaceVect, unitizedSurfSpaceVect, &dist);

    // get a normalized surface sun vector
    SpiceDouble surfaceSunVect[3];
    vsub_c(m_uB, pB, surfaceSunVect);
    SpiceDouble unitizedSurfSunVect[3];
    unorm_c(surfaceSunVect, unitizedSurfSunVect, &dist);

    // use normalized surface spacecraft and surface sun vectors to calculate
    // the phase angle (in radians)
    phase = Angle(vsep_c(unitizedSurfSpaceVect, unitizedSurfSunVect),
        Angle::Radians);

    // use normalized surface spacecraft and local normal vectors to calculate
    // the emission angle (in radians)
    emission = Angle(vsep_c(unitizedSurfSpaceVect, normal),
        Angle::Radians);

    // use normalized surface sun and normal vectors to calculate the incidence
    // angle (in radians)
    incidence = Angle(vsep_c(unitizedSurfSunVect, normal),
        Angle::Radians);
  
    return incidence;
    
}
