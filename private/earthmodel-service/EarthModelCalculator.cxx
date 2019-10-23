/**
    copyright  (C) 2004
    the icecube collaboration
    $Id: EarthModelCalculator.cxx 86266 2012-03-13 15:40:52Z claudio.kopper $

    @version $Revision: 1.1 $
    @date $Date: 2012-03-13 10:40:52 -0500 (Tue, 13 Mar 2012) $
    @author

    @todo

*/
#include "earthmodel-service/EarthModelCalculator.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Constants.h"

using namespace std;
//------------------------------------------------------
// Oct.22 2012 K.Hoshina
// Following params and functions are used for finding 
// intersections between a cylinder and a particle.
// They are techinical tolerances and internal functions,
// so that I define it here as static for internal linkage.
//------------------------------------------------------

namespace earthmodel {

static double TOLERANCE = 1e-6;

//----------------------------------------
// namespace EarthModelCalculator start
//----------------------------------------

//__________________________________________________________
double EarthModelCalculator::GetImpactParameter(const I3Position &p0,
                                      const I3Direction &d,
                                      double &t,
                                      I3Position &p)
{
//
// calculate impact parameter to origin and scalar value t
// that fulfills impact position p = p0 + t*d.
// this function should work in any Cartecian coordinate
//

  // p  : most closest position from origin to a track
  // p0 : particle position
  // d  : particle direction (unit vector)
  // Now, vector p should be perpendicular to direction d:
  // p = p0 + td;
  // p*d = 0 => (p0 + td) * d = 0 => t = -p0*d

  t = -1.0 * p0*d;
  p = p0 + d*t;
  double r = p.Magnitude();

  return r;

}

//_________________________________________________________
int EarthModelCalculator::GetIntersectionsWithSphere(
                             const I3Position &pos,
                             const I3Direction &dir,
                             double r,
                             I3Position &startPos,
                             I3Position &endPos)
{

  // calc impact position
  I3Position impact_pos;
  double t = 0;
  double impact_param = GetImpactParameter(pos, dir, t, impact_pos);

  if (fabs(impact_pos* dir) > TOLERANCE) {
     log_warn("impact_pos must be most closest position "
               "to origin, but your impact_pos vector and "
               "dir are not parpendicular. Check params. ");
     return 0;
  }

  if (impact_param > r) {
     /*
     printf("impact_param %e is larger than r %e\n",
               impact_param, r);
     */
     return 0;
  } else if (impact_param == r) {
     /*
     printf("impact_param %e is same as r %e,"
               "startpos=endpos=(%e, %e, %e)\n",
                impact_param, r, impact_pos.GetX(),
                impact_pos.GetY(), impact_pos.GetZ());
     */
     startPos = impact_pos;
     endPos = impact_pos;
     return 1;
  }

  //
  // we obtain parameter t with the simple Pythagorean Proposition
  // impact_param^2 + t^2 = radius^2
  // t^2 = (radius + impact_param)*(radius - impact_param)

  t = sqrt((r + impact_param)*(r - impact_param));

  startPos = impact_pos - dir*t;
  endPos = impact_pos + dir*t;

  /*
  printf("impact_param %e is smaller than r %e,"
            "startpos=(%e, %e, %e), endpos=(%e, %e, %e)\n",
                impact_param, r, startPos.GetX(),
                startPos.GetY(), startPos.GetZ(),
                endPos.GetX(), 
                endPos.GetY(), startPos.GetZ());
  */

  return 2;

}

//_________________________________________________________
int EarthModelCalculator::GetDistsToIntersectionsWithSphere(
                             const I3Position &pos,
                             const I3Direction &dir,
                             double r,
                             double &enterdist,
                             double &exitdist)
{

   I3Position startp, endp;
   int n = GetIntersectionsWithSphere(pos, dir, r, startp, endp);

   if (n == 0) return 0; 

   enterdist = (startp - pos)*dir;
   exitdist  = (endp - pos)*dir;
   return n;
}

//_____________________________________________________________
double EarthModelCalculator::GetLeptonRange(double particle_energy, 
                                   bool   isTau,
                                   LeptonRangeOption option,
                                   double scale)
{
  //printf("called GetLeptonRange with e=%f, istau=%d, option=%f\n",
  //        particle_energy, (int)isTau, option);

  double range=0.;

  // First, calculate muon range

  if (option == LEGACY) {
     // use legacy parameter
     // rough parametrization of muon range (mwe):

     range=3e3*log(1+particle_energy/600);
     if (range>1e5) range=1e5;  //100,000 =100kmwe max rage

  } else if (option == DEFAULT) {
     // 10 < option < 20
     // use original parameters from Dima's internal report (2011 Mar 29) 
     static const double new_dima_a = 0.212/1.2;
     static const double new_dima_b = 0.251e-3/1.2;
     range = log(1 + particle_energy * new_dima_b/new_dima_a) / new_dima_b;
     
  } else if (option == NUSIM) {
     // 20 < option < 30
     // use Gary's NuSim parameter
     // but converted from [cm.w.e] to [m.w.e]
     if (particle_energy < 1e3) {
        range = particle_energy/2.0e-5;
     } else if (particle_energy < 1e4) {
        range = 1.0e3 * (3.5 + 9.0 * (log10(particle_energy) - 3.0));
     } else {
        range = 1.0e3 * (12.0 + 6.0 * (log10(particle_energy) - 4.0));
     }

  } else {
     log_error("The LeptonRange option %d is not supported! ", option);
  }

  if (isTau) {
    // currently we have only Pat's parameter.
    // rough parametrization of tau range (mwe):
    //range=38e3*log(1+particle_energy/5.6e7);
    // rough parametrization of tau->muon range (mwe):

    range=38e3*log(1+particle_energy/5.6e7) + range;
    if (range>3e5) range=3e5;//300,000 =300kmwe max rage
    //log_error("Particle is a tau: range = %f", range);
  }

  return range * scale;  // [m.w.e]

}

//_____________________________________________________________
double EarthModelCalculator::ColumnDepthCGStoMWE(double cdep_CGS) 
{
   return cdep_CGS/100.;   
}

//_____________________________________________________________
double EarthModelCalculator::MWEtoColumnDepthCGS(double cdep_MWE) 
{
   return cdep_MWE*100.;   
}

}

