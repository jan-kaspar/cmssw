/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *	Jan Kašpar (jan.kaspar@gmail.com)
 *
 ****************************************************************************/

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include "DataFormats/CTPPSAlignment/interface/RPAlignmentCorrectionData.h"
#include "DataFormats/CTPPSAlignment/interface/RPAlignmentCorrectionsData.h"

#include "Alignment/CTPPSTrackBased/interface/CommonMethods.h"
#include "Alignment/CTPPSTrackBased/interface/AlignmentGeometry.h"

#include <set>

#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

/**
 * NOTE ON ERROR PROPAGATION
 *
 * It is not possible to split (and merge again) the experimental errors between the RP and sensor
 * contributions. To do so, one would need to keep the entire covariance matrix. Thus, it has been
 * decided to save:
 *   RP errors = the uncertainty of the common shift/rotation
 *   sensor error = the full experimental uncertainty
 * In consequence: RP and sensor errors SHOULD NEVER BE SUMMED!
 **/
void CommonMethods::FactorRPFromSensorCorrections(const RPAlignmentCorrectionsData &inputAlignments,
  RPAlignmentCorrectionsData &expandedAlignments, RPAlignmentCorrectionsData &factoredAlignments,
  const AlignmentGeometry &geometry, bool equalWeights, unsigned int verbosity)
{
  // clean first
  expandedAlignments.clear();
  factoredAlignments.clear();

  // save full sensor alignments
  map<unsigned int, RPAlignmentCorrectionData> fullAlignments;
  map<unsigned int, set<unsigned int> > sensorsPerRP;
  for (auto it : inputAlignments.getSensorMap())
  {
    const auto &sensorId = it.first;

    // RP errors are coming from the previous iteration and shall be discarded!
    fullAlignments[sensorId] = inputAlignments.getFullSensorCorrection(sensorId, false);

    sensorsPerRP[CTPPSDetId(sensorId).getRPId()].insert(sensorId);
  }

  // convert full alignments to expandedAlignments
  for (const auto it : fullAlignments)
  {
    expandedAlignments.setSensorCorrection(it.first, it.second);
  }

  // do the factorization RP per RP
  for (auto rpit : sensorsPerRP)
  {
    CTPPSDetId rpId(rpit.first);
    const set<unsigned int> &sensors = rpit.second;

    if (verbosity)
      printf("* processing RP %u\n", rpit.first);

    // process strip RPs
    if (rpId.subdetId() == CTPPSDetId::sdTrackingStrip)
    {
      // TODO: this approach ignores uncertainties/weights

      // extract mean shx, shy and rotz
      const unsigned int N = sensors.size();

      TMatrixD B(N, 2), Vi(N, N);
      TVectorD M(N);

      double S_rotz=0.;

      unsigned int idx = 0;
      for (const auto &senId : sensors)
      {
        auto d2 = geometry.Get(senId).GetDirectionData(2);

        B(idx, 0) = d2.dx;
        B(idx, 1) = d2.dy;
        
        M(idx) = d2.dx * fullAlignments[senId].getShX() + d2.dy * fullAlignments[senId].getShY();

        Vi(idx, idx) = 1.;
        
        S_rotz += fullAlignments[senId].getRotZ();

        ++idx;
      }

      TMatrixD BT(TMatrixD::kTransposed, B);
      TMatrixD BTViB(BT, TMatrixD::kMult, Vi * B);
      TMatrixD BTViBi(TMatrixD::kInverted, BTViB);
	  TVectorD th_B(2);
      th_B = BTViBi * BT * Vi * M;

      double m_shx = th_B[0];
      double m_shy = th_B[1];
      double m_rotz = S_rotz / N;

      printf("    m_shx = %.3f, m_shy = %.3f, m_rotz = %.3f\n", m_shx, m_shy, m_rotz);

      factoredAlignments.addRPCorrection(rpId, RPAlignmentCorrectionData(m_shx, m_shy, 0., 0., 0., m_rotz));

      // calculate residuals
      for (const auto &senId : sensors)
      {
        auto d2 = geometry.Get(senId).GetDirectionData(2);

        const double de_s = d2.dx * (fullAlignments[senId].getShX() - m_shx) + d2.dy * (fullAlignments[senId].getShY() - m_shy);

        factoredAlignments.addSensorCorrection(senId, RPAlignmentCorrectionData(
          d2.dx * de_s,
          d2.dy * de_s,
          0.,
          0.,
          0.,
          fullAlignments[senId].getRotZ() - m_rotz
        ));
      }
    }

    // process pixel RPs
    if (rpId.subdetId() == CTPPSDetId::sdTrackingPixel)
    {
      // TODO: this approach ignores uncertainties/weights

      // extract mean shx, shy and rotz
      double S_1=0., S_shx=0., S_shy=0., S_rotz=0.;
      for (const auto &senId : sensors)
      {
        S_1 += 1.;
        S_shx += fullAlignments[senId].getShX();
        S_shy += fullAlignments[senId].getShY();
        S_rotz += fullAlignments[senId].getRotZ();
      }
    
      double m_shx = S_shx / S_1;
      double m_shy = S_shy / S_1;
      double m_rotz = S_rotz / S_1;

      printf("    S_1 = %.3f, m_shx = %.3f, m_shy = %.3f, m_rotz = %.3f\n", S_1, m_shx, m_shy, m_rotz);

      factoredAlignments.addRPCorrection(rpId, RPAlignmentCorrectionData(m_shx, m_shy, 0., 0., 0., m_rotz));

      // calculate residuals
      for (const auto &senId : sensors)
      {
        factoredAlignments.addSensorCorrection(senId, RPAlignmentCorrectionData(
          fullAlignments[senId].getShX() - m_shx,
          fullAlignments[senId].getShY() - m_shy,
          0.,
          0.,
          0.,
          fullAlignments[senId].getRotZ() - m_rotz
        ));
      }
    }
  }
}
