/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "Alignment/CTPPSTrackBased/interface/AlignmentGeometry.h"
#include "Alignment/CTPPSTrackBased/interface/Utilities.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

const DetGeometry& AlignmentGeometry::Get(unsigned int id) const
{
  auto it = sensorGeometry.find(id);
  if (it == sensorGeometry.end())
    throw cms::Exception("AlignmentGeometry") << "No geometry available for sensor " << id << ".";

  return it->second;
}

//----------------------------------------------------------------------------------------------------

void AlignmentGeometry::Insert(unsigned int id, const DetGeometry &g)
{
  sensorGeometry[id] = g;
}

//----------------------------------------------------------------------------------------------------

void AlignmentGeometry::Print() const
{
  for (const auto &it : sensorGeometry)
  {
    PrintId(it.first);

    printf(" z = %+10.4f mm │ shift: x = %+7.3f mm, y = %+7.3f mm │ ",
        it.second.z,
        it.second.sx, it.second.sy);

    for (const auto &dit : it.second.directionData)
    {
      printf("dir%u: %+.3f, %+.3f, %+.3f │ ", dit.first, dit.second.dx, dit.second.dy, dit.second.dz);
    }

    if (CTPPSDetId(it.first).subdetId() == CTPPSDetId::sdTrackingStrip)
      printf("%s", (it.second.isU) ? "U-det" : "V-det");

    printf("\n");
  }
}
