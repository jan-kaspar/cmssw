/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#ifndef Alignment_CTPPSTrackBased_Utilities_h
#define Alignment_CTPPSTrackBased_Utilities_h

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include <TMatrixD.h>

extern void PrintId(unsigned int id);

extern void Print(TMatrixD& m, const char *label = NULL, bool mathematicaFormat = false);

#endif
