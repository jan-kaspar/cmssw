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

//----------------------------------------------------------------------------------------------------

struct Hit
{
  /// sensor id
  unsigned int id;

  /// index of read-out direction (valid are: 1 or 2)
  unsigned int dirIdx;

  /// measurement position; mm
  double position;

  /// measurement position; mm
  double sigma;

  // TODO: estimate of z global ??
  //  for strips and diamond: copy constant from geometry
  //  for pixels: evaluate global z from local vector (m1, m2, 0)

  Hit(unsigned int _id=0, unsigned int _dirIdx=0, double _pos=0, double _sig=0) : id(_id), dirIdx(_dirIdx), position(_pos), sigma(_sig)
  {
  }
};

//----------------------------------------------------------------------------------------------------

typedef std::vector<Hit> HitCollection;

//----------------------------------------------------------------------------------------------------

#endif
