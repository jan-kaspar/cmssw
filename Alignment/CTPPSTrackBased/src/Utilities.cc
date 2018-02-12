/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "Alignment/CTPPSTrackBased/interface/Utilities.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSPixelDetId.h"

#include <cmath>

//----------------------------------------------------------------------------------------------------

void PrintId(unsigned int id)
{
  CTPPSDetId detId(id);

  if (detId.subdetId() == CTPPSDetId::sdTrackingStrip)
  {
    TotemRPDetId stDetId(id);
    printf("strip %u (%3u.%u)", id, 100*stDetId.arm() + 10*stDetId.station() + stDetId.rp(), stDetId.plane());
  }

  if (detId.subdetId() == CTPPSDetId::sdTimingDiamond)
  {
    CTPPSDiamondDetId diDetId(id);
    printf("dimnd %u (%3u.%u)", id, 100*diDetId.arm() + 10*diDetId.station() + diDetId.rp(), diDetId.plane());
  }

  if (detId.subdetId() == CTPPSDetId::sdTrackingPixel)
  {
    CTPPSPixelDetId piDetId(id);
    printf("pixel %u (%3u.%u)", id, 100*piDetId.arm() + 10*piDetId.station() + piDetId.rp(), piDetId.plane());
  }
}

//----------------------------------------------------------------------------------------------------

void Print(TMatrixD& m, const char *label, bool mathematicaFormat)
{
  if (mathematicaFormat) {
    printf("{");
    for (int i = 0; i < m.GetNrows(); i++) {
      if (i > 0) printf(", ");
      printf("{");
      for (int j = 0; j < m.GetNcols(); j++) {
        if (j > 0) printf(", ");
        printf("%.3f", m[i][j]);
      }
      printf("}");
    }
    printf("}\n");
    return;
  }

  if (label)
    printf("\n%s\n", label);

  printf("    | ");
  for (int j = 0; j < m.GetNcols(); j++)
    printf(" %9i", j);
  printf("\n------");
  for (int j = 0; j < m.GetNcols(); j++)
    printf("----------");
  printf("\n");

  for (int i = 0; i < m.GetNrows(); i++) {
    printf("%3i | ", i);
    for (int j = 0; j < m.GetNcols(); j++) {
      double v = m[i][j];
      if
        (fabs(v) >= 1E4) printf(" %+9.2E", v);
      else
        if (fabs(v) > 1E-6)
          printf(" %+9.2E", v);
        else
          printf("         0");
    }
    printf("\n");
  }
}
