/* vim: set ts=4: */
/**
 * copyright  (C) 2004
 * the icecube collaboration
 * $Id: EarthModelCalculator.h $
 *
 * @file EarthModelCalculator.h
 * @version $Revision: 1.16 $
 * @date $Date: 2012-03-13 10:40:52 -0500 (Tue, 13 Mar 2012) $
 * @author Kotoyo Hoshina <hoshina@icecube.wisc.edu>
 */
#pragma once
#ifndef EarthModelCalculator_h
#define EarthModelCalculator_h

#include <cmath>
#include <vector>
#include <cassert>

#include "LeptonInjector/math/Coordinates.h"
#include "LeptonInjector/utilities/Constants.h"

/**
 * @brief This is a namespace which provides a collection of stand-alone
 * functions that calculate various geometrical information
 * for particle propagation of the Earth.
 */

namespace LI {
namespace detector {
namespace EarthModelCalculator {
   enum LeptonRangeOption {DEFAULT, LEGACY, NUSIM};

  /**
   * calculate impact parameter with rewpect to origin
   * and scalar value t that fulfills p = p0 + t*d,
   * where p0 is start position of a track, d is
   * direction of the track, and p is most closest position
   * on a track from origin.
   * this function should work in any Cartecian coordinate
   *
   * @param[in] p0  particle start position
   *
   * @param[in] d  particle direction (unit vector)
   *
   * @param[out] t distance from p0 to p
   *
   * @param[out] p most closest position on a track from origin
   *
   * @return impact parameter, distance from origin to p
   *
   */
   double GetImpactParameter(const LeptonInjector::LI_Position &p0,
                             const LeptonInjector::LI_Direction &d,
                             double &t,
                             LeptonInjector::LI_Position &p);

  /**
   * This function returns intersection-positions between a track
   * and a sphere with radius r.
   * Note that the origin of track position and direction must be
   * at the center of the sphere. Return positions are in same
   * coordinate as input.
   * If there is only one intersection, it will be stored in both
   * output parameters.
   *
   * @param[in] pos track position
   * @param[in] dir track direction (unit vector)
   * @param[in] r   radius
   *
   * @param[out] startPos  track-entering-position to the sphere
   * @param[out] endPos    track-exitting-position from the sphere
   *
   * @return number of intersections
   */
   int GetIntersectionsWithSphere(
                           const LeptonInjector::LI_Position &pos,
                           const LeptonInjector::LI_Direction &dir,
                           double r,
                           LeptonInjector::LI_Position &startPos,
                           LeptonInjector::LI_Position &endPos);

  /**
   * wrapper function of GetIntersectionsWithSphere
   *
   * @param[in] pos track position
   * @param[in] dir track direction (unit vector)
   * @param[in] r   radius
   *
   * @param[out] enterdist distance from pos to track-entering-position
   *             negative value indicates behind the pos
   *
   * @param[out] exitdist distance from pos to track-exiting-position
   *             negative value indicates behind the pos
   *
   * @return number of intersections
   */
   int GetDistsToIntersectionsWithSphere(
                           const LeptonInjector::LI_Position &pos,
                           const LeptonInjector::LI_Direction &dir,
                           double r,
                           double & enterdist,
                           double & exitdist);

  /**
   * @brief Returns muon range in m.w.e.
   * If you need surviving length [m] of muon/tau, use
   * EarthModelService::GetLeptonRangeInMeter().
   *
   * @return range [m.w.e]
   *
   * Now MuonRange calculation offers three options:
   *
   * DEFAULT
   *    -> use Dima's fitting function and optimized parameter
   *       confirmed on Mar. 29, 2011
   *
   * LEGACY
   *    -> use legacy nugen equation and parameter
   *       obtained study with mmc, using ice medium
   *
   * NUSIM
   *    -> use Gary's fitting function (used in NUSIM)
   *
   * scale gives
   * a scaling factor for MuonRange in Meter.
   * See GetLeptonRangeInMeter().
   *
   * [Dima's equation]
   *
   * muon and tau: R = (1/b)*log[(bE/a)+1]
   * for the equation see arXiv:hep-ph/0407075, section 6.2 P16.
   *
   * muon: (legacy nugen parameter)
   * b~3.4*10^-4 [1/m]
   * a~2*10^-1 [GeV/m]
   * R_mu ~ 3000 * log(1.0 + 1.7*10^-3*E[GeV]) [mwe]
   *
   * tau:  (see: http://www.ice.phys.psu.edu/~toale/icecube/nutau/nutau-mc.html)
   * b~2.6*10^-5 [1/m]
   * a~1.5*10^3 [GeV/m]
   * R_tau ~ 38000 * log(1.0 + 1.8*10^-8*E[GeV]) [mwe]
   *
   *
   * [Gary's equation] (muon only)
   *
   * if (energy < 1e3) {
   *    R_mu = energy/0.002;
   * } else if (energy < 1e4) {
   *    R_mu = 1.0e5 * (3.5 + 9.0 * (TMath::Log10(energy) - 3.0));
   * } else {
   *    R_mu = 1.0e5 * (12.0 + 6.0 * (TMath::Log10(energy) - 4.0));
   * }
   *
   */
   double GetLeptonRange(double particle_energy,
                       bool   isTau = false,
                       LeptonRangeOption option = DEFAULT,
                       double scale = 1);

   /**
    * @brief unit conversion from g/cm2 to m.w.e
    */
   double ColumnDepthCGStoMWE(double cdep_CGS);

   /**
    * @brief unit conversion from m.w.e to g/cm2
    */
   double MWEtoColumnDepthCGS(double range_MWE);

}


} // namespace detector
} // namespace LI

#endif
