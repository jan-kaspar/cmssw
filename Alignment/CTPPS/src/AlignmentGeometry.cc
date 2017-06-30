/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "Alignment/CTPPS/interface/AlignmentGeometry.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

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
    const DetGeometry &d = it->second;
    printf("%u\t%+E\t%+E\t%+E\t%+E\t%+E\t%+E\t%+E\n", it->first, d.z, d.sx, d.sy, d.d1x, d.d1y, d.d2x, d.d2y);
  }
}

//----------------------------------------------------------------------------------------------------

void AlignmentGeometry::LoadFromFile(const std::string &filename)
{
  clear();

  FILE *f = fopen(filename.c_str(), "r");
  if (!f)
    throw cms::Exception("AlignmentGeometry::LoadFromFile") << "File `" << filename << "' can not be opened." << endl;

  while (!feof(f))
  {
    unsigned int id;
    float x, y, z, d1x, d1y, d2x, d2y;

    int res = fscanf(f, "%u%E%E%E%E%E%E%E", &id, &x, &y, &z, &d1x, &d1y, &d2x, &d2y);

    if (res == 8)
    {
      bool isU = false;

      CTPPSDetId detId(id);
      if (detId.subdetId() == CTPPSDetId::sdTrackingStrip)
      {
        TotemRPDetId stripDetId(id);
        const unsigned int rpNum = stripDetId.rp();
        const unsigned int plNum = stripDetId.plane();
        isU = (plNum % 2 != 0);
        if (rpNum == 2 || rpNum == 3)
          isU = !isU;
      }

      Insert(id, DetGeometry(z, x, y, d1x, d1y, d2x, d2y, isU));      
    } else {
      if (!feof(f))
      {
        throw cms::Exception("AlignmentGeometry::LoadFromFile") << "Cannot parse file `" << filename
          << "'. The format is probably wrong." << endl;
      }
    }
  }

  fclose(f);
}
