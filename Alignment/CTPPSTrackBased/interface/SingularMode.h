/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#ifndef Alignment_CTPPSTrackBased_SingularMode_h
#define Alignment_CTPPSTrackBased_SingularMode_h

#include <TVectorD.h>

/**
 *\brief 
 **/
struct SingularMode
{
  /// eigen value
  double val;

  /// eigen vector
  TVectorD vec;

  /// index
  unsigned int idx;
};

#endif
