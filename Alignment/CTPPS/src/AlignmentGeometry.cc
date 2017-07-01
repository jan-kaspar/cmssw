/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "Alignment/CTPPS/interface/AlignmentGeometry.h"

#include "Alignment/CTPPS/interface/utilities.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

void AlignmentGeometry::Insert(unsigned int id, const DetGeometry &g)
{
  insert(value_type(id, g));
  
  rps.insert(CTPPSDetId(id).getRPId());
}

//----------------------------------------------------------------------------------------------------

// TODO: remove
/*
unsigned int AlignmentGeometry::MatrixIndexToDetId(unsigned int mi) const
{
  const_iterator it = FindByMatrixIndex(mi);
  if (it != end())
  {
    return it->first;
  } else {
    LogProblem("AlignmentGeometry") << ">> AlignmentGeometry::MatrixIndexToDetId > No detector corresponds to matrix index "
      << mi << ".";
    return 0;
  }
}
*/

//----------------------------------------------------------------------------------------------------

AlignmentGeometry::const_iterator AlignmentGeometry::FindByMatrixIndex(unsigned int mi) const
{
  for (const_iterator it = begin(); it != end(); ++it)
  {
    if (it->second.matrixIndex == mi)
      return it;
  }
  return end();
}

//----------------------------------------------------------------------------------------------------

void AlignmentGeometry::Print() const
{
  for (const_iterator it = begin(); it != end(); ++it) 
  {
    PrintId(it->first);

    printf(" [%2u] z = %+10.4f mm │ shift: x = %+7.3f mm, y = %+7.3f mm │ ",
        it->second.matrixIndex, it->second.z,
        it->second.sx, it->second.sy);

    for (const auto &dit : it->second.directionData)
    {
      printf("dir%u: %+.3f, %+.3f │ ", dit.first, dit.second.dx, dit.second.dy);
    }

    if (CTPPSDetId(it->first).subdetId() == CTPPSDetId::sdTrackingStrip)
      printf("%s", (it->second.isU) ? "U-det" : "V-det");

    printf("\n");
  }
}
