/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#ifndef _SingularMode_h_
#define _SingularMode_h_

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
