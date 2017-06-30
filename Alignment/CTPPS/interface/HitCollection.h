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
  /// sensor id
  unsigned int id;

  /// first measurement: position and uncertainty in mm
  double pos1, sig1;

  /// second measurement: position and uncertainty in mm
  double pos2, sig2;

  Hit(unsigned int _id=0, double _p1=0., double _s1=0., double _p2=0., double _s2=0.) : id(_id), pos1(_p1), sig1(_s1), pos2(_p2), sig2(_s2)
  {
  }
};

typedef std::vector<Hit> HitCollection;

#endif
