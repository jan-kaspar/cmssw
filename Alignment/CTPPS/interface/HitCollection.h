/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#ifndef _HitCollection_h_
#define _HitCollection_h_

#include <vector>

class RPRecoHit;

struct Hit
{
  unsigned int id;    ///< detector decimal id
  double position;    ///< position in mm
  double sigma;       ///< uncertainty in mm

  Hit(unsigned int i=0, double p=0., double s=0.) : id(i), position(p), sigma(s) {}
};

typedef std::vector<Hit> HitCollection;

#endif
