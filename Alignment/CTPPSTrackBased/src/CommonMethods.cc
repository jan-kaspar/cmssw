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

      double S_sh_z=0.;
      double S_rot_x=0., S_rot_y=0., S_rot_z=0.;

      unsigned int idx = 0;
      for (const auto &senId : sensors)
      {
        auto d2 = geometry.Get(senId).GetDirectionData(2);

        B(idx, 0) = d2.dx;
        B(idx, 1) = d2.dy;
        
        M(idx) = d2.dx * fullAlignments[senId].getShX() + d2.dy * fullAlignments[senId].getShY();

        Vi(idx, idx) = 1.;

        S_sh_z += fullAlignments[senId].getShZ();
        
        S_rot_x += fullAlignments[senId].getRotX();
        S_rot_y += fullAlignments[senId].getRotY();
        S_rot_z += fullAlignments[senId].getRotZ();

        ++idx;
      }

      TMatrixD BT(TMatrixD::kTransposed, B);
      TMatrixD BTViB(BT, TMatrixD::kMult, Vi * B);
      TMatrixD BTViBi(TMatrixD::kInverted, BTViB);
	  TVectorD th_B(2);
      th_B = BTViBi * BT * Vi * M;

      const double m_sh_x = th_B[0];
      const double m_sh_y = th_B[1];
      const double m_sh_z = S_sh_z / N;

      const double m_rot_x = S_rot_x / N;
      const double m_rot_y = S_rot_y / N;
      const double m_rot_z = S_rot_z / N;

      printf("    m_sh_x = %.3f, m_sh_y = %.3f, m_sh_z = %.3f, m_rot_x = %.3f, m_rot_y = %.3f, m_rot_z = %.3f\n",
        m_sh_x, m_sh_y, m_sh_z, m_rot_x, m_rot_y, m_rot_z);

      factoredAlignments.addRPCorrection(rpId, RPAlignmentCorrectionData(m_sh_x, m_sh_y, m_sh_z, m_rot_x, m_rot_y, m_rot_z));

      // calculate residuals
      for (const auto &senId : sensors)
      {
        auto d2 = geometry.Get(senId).GetDirectionData(2);

        const double de_s = d2.dx * (fullAlignments[senId].getShX() - m_sh_x) + d2.dy * (fullAlignments[senId].getShY() - m_sh_y);

        factoredAlignments.addSensorCorrection(senId, RPAlignmentCorrectionData(
          d2.dx * de_s,
          d2.dy * de_s,
          fullAlignments[senId].getShZ() - m_sh_z,
          fullAlignments[senId].getRotX() - m_rot_x,
          fullAlignments[senId].getRotY() - m_rot_y,
          fullAlignments[senId].getRotZ() - m_rot_z
        ));
      }
    }

    // process pixel RPs
    if (rpId.subdetId() == CTPPSDetId::sdTrackingPixel)
    {
      // TODO: this approach ignores uncertainties/weights

      // extract mean shx, shy and rotz
      // TODO: add other components
      double S_1=0.;
      double S_sh_x=0., S_sh_y=0., S_sh_z=0.;
      double S_rot_x=0., S_rot_y=0., S_rot_z=0.;
      for (const auto &senId : sensors)
      {
        const auto &a = fullAlignments[senId];

        S_1 += 1.;

        S_sh_x += a.getShX();
        S_sh_y += a.getShY();
        S_sh_z += a.getShZ();

        S_rot_x += a.getRotX();
        S_rot_y += a.getRotY();
        S_rot_z += a.getRotZ();
      }
    
      const double m_sh_x = S_sh_x / S_1;
      const double m_sh_y = S_sh_y / S_1;
      const double m_sh_z = S_sh_y / S_1;

      const double m_rot_x = S_rot_x / S_1;
      const double m_rot_y = S_rot_y / S_1;
      const double m_rot_z = S_rot_z / S_1;

      printf("    m_sh_x = %.3f, m_sh_y = %.3f, m_sh_z = %.3f, m_rot_x = %.3f, m_rot_y = %.3f, m_rot_z = %.3f\n",
        m_sh_x, m_sh_y, m_sh_z, m_rot_x, m_rot_y, m_rot_z);

      factoredAlignments.addRPCorrection(rpId, RPAlignmentCorrectionData(m_sh_x, m_sh_y, m_sh_z, m_rot_x, m_rot_y, m_rot_z));

      // calculate residuals
      for (const auto &senId : sensors)
      {
        factoredAlignments.addSensorCorrection(senId, RPAlignmentCorrectionData(
          fullAlignments[senId].getShX() - m_sh_x,
          fullAlignments[senId].getShY() - m_sh_y,
          fullAlignments[senId].getShZ() - m_sh_z,
          fullAlignments[senId].getRotX() - m_rot_x,
          fullAlignments[senId].getRotY() - m_rot_y,
          fullAlignments[senId].getRotZ() - m_rot_z
        ));
      }
    }
  }
}
