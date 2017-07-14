/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *	Jan Kašpar (jan.kaspar@gmail.com)
 *
 ****************************************************************************/

#include "Alignment/CTPPS/interface/CommonMethods.h"

// TODO: clean
/*
#include "FWCore/Utilities/interface/typelookup.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
*/

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include "DataFormats/CTPPSAlignment/interface/RPAlignmentCorrectionData.h"

#include "Alignment/CTPPS/interface/AlignmentGeometry.h"


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
void CommonMethods::FactorRPFromSensorCorrections(const RPAlignmentCorrectionsData &input,
  RPAlignmentCorrectionsData &expanded, RPAlignmentCorrectionsData &factored,
  const AlignmentGeometry &geometry, bool equalWeights, unsigned int verbosity)
{
  // TODO
#if 0
  // clean first
  expanded.Clear();
  factored.Clear();

  const auto &sensors = input.GetSensorMap();
  const auto &rps = input.GetRPMap();

  // save full sensor alignments
  map<unsigned int, RPAlignmentCorrectionData> origAlignments;
  map<unsigned int, set<unsigned int> > sensorsPerRP;
  for (auto it = sensors.begin(); it != sensors.end(); ++it)
  {
    const unsigned int sensorId = it->first;

    auto git = geometry.find(sensorId);
    if (git == geometry.end())
      continue;

    // RP errors are coming from the previous iteration and shall be discarded!
    origAlignments[sensorId] = input.GetFullSensorCorrection(sensorId, false);

    const DetGeometry &g = git->second;
    auto d1 = g.GetDirectionData(1);
    auto d2 = g.GetDirectionData(2);

    origAlignments[sensorId].xyTranslationToReadout(d1.dx, d1.dy, d2.dx, d2.dy);

    sensorsPerRP[CTPPSDetId(sensorId).getRPId()].insert(sensorId);
  }

  // do the factorization RP per RP
  for (auto it = sensorsPerRP.begin(); it != sensorsPerRP.end(); ++it)
  {
    unsigned int rpId = it->first;
    const set<unsigned int> &sensors = it->second;

    if (verbosity)
      printf("* processing RP %u\n", rpId);

    // get z0: typical z of the RP
    unsigned int N = 0;
    double z0 = 0;
    for (auto dit = sensors.begin(); dit != sensors.end(); ++dit)
    {
      auto git = geometry.find(*dit);
      const DetGeometry &d = git->second;
      N++;
      z0 += d.z;
    }
    z0 /= N;

    if (verbosity > 1)
      printf("\tN=%u, z0 = %E\n", N, z0);

    // skip RPs not listed in the geometry
    if (N == 0)
      continue;

    // shift fit variables
    TMatrixD A(N, 4), B(N, 2), V(N, N), Vi(N, N);
    TVectorD m(N);

    // rotation fit variables
    double Sr = 0., S1 = 0., Sss = 0.;

    // fit the shifts and rotations
    unsigned int idx = 0;
    for (auto dit = sensors.begin(); dit != sensors.end(); ++dit)
    {
      auto git = geometry.find(*dit);
      const DetGeometry &d = git->second;
      const RPAlignmentCorrectionData &oa = origAlignments[*dit];

      // shifts part
      double sh_r = oa.sh_r();
      double sh_r_e = oa.sh_r_e();

      // TODO: improve comments
      if (sh_r_e <= 0.)
        sh_r_e = 1E-8; // in mm
                        // 1E-8 seems to be a good value. It is significantly smaller
                        // than usual errors, but doesn't cause numerical problems like
                        // values below 1E-11

      double zeff = d.z - z0;

      A(idx, 0) = d.dx*zeff;
      A(idx, 1) = d.dx;
      A(idx, 2) = d.dy*zeff;
      A(idx, 3) = d.dy;

      B(idx, 0) = d.dx;
      B(idx, 1) = d.dy;

      V(idx, idx) = sh_r_e*sh_r_e;
      Vi(idx, idx) = (equalWeights) ? 1. : 1./sh_r_e/sh_r_e;
      m(idx) = sh_r;

      // rotations part
      double rot_z = oa.rot_z();
      double rot_z_e = oa.rot_z_e();
      if (rot_z_e <= 0.)
        rot_z_e = 1E-8; // rad

      double w = (equalWeights) ? 1. : 1. / rot_z_e / rot_z_e;
      Sr += rot_z * w;
      S1 += 1. * w;
      Sss += rot_z_e * rot_z_e;

      //printf("%u %u | %.3f +- %.3f | %.3f +- %.3f\n", *dit, idx, sh_r*1E3, sh_r_e*1E3, rot_z*1E3, rot_z_e*1E3);

      idx++;
    }

    // linear shift fit
    TMatrixD AT(TMatrixD::kTransposed, A);
    TMatrixD VRi(TMatrixD::kInverted, V);
    TMatrixD ATVRiA(AT, TMatrixD::kMult, VRi * A);
    TMatrixD ATVRiAi(ATVRiA);
    try {
      ATVRiAi = ATVRiA.Invert();
    }
    catch (...) {
      printf("ERROR in RPAlignmentCorrections::FactorRPFromSensorCorrections > AT A matrix is singular, skipping RP %u.\n", rpId);
      continue;
    }

    TVectorD th(4);
    th = ATVRiAi * AT * VRi * m;

    // g: intercepts (mm), h: slopes (rad), with errors
    double hx = th[0], hx_error = sqrt(ATVRiAi(0, 0));
    double gx = th[1], gx_error = sqrt(ATVRiAi(1, 1));
    double hy = th[2], hy_error = sqrt(ATVRiAi(2, 2));
    double gy = th[3], gy_error = sqrt(ATVRiAi(3, 3));

    // constant shift fit
    TMatrixD BT(TMatrixD::kTransposed, B);
    TMatrixD BTViB(BT, TMatrixD::kMult, Vi * B);
    TMatrixD BTViBi(TMatrixD::kInverted, BTViB);

    TMatrixD V_th_B_eW(BTViBi * BT * V * B * BTViBi);
    TMatrixD &V_th_B = (equalWeights) ? V_th_B_eW : BTViBi;

    TVectorD th_B(2);
    th_B = BTViBi * BT * Vi * m;
    double g0x = th_B[0], g0x_error = sqrt(V_th_B(0, 0));
    double g0y = th_B[1], g0y_error = sqrt(V_th_B(1, 1));

    // const rotation fit
    double rot_z_mean = Sr / S1;
    double rot_z_mean_error = (equalWeights) ? sqrt(Sss)/S1 : sqrt(1. / S1);

    // shift corrections
    TVectorD sc(B * th_B);

    // corrected/internal shift error matrix
    TMatrixD VR(V);
    VR -= B * BTViBi * BT;

    if (verbosity)
    {
      printf("\tshift fit\n");
      printf("\t\tconstant: gx=%.2E +- %.2E um, gy=%.2E +- %.2E um\n",
        g0x*1E3, g0x_error*1E3, g0y*1E3, g0y_error*1E3);
      printf("\t\tlinear  : gx=%.2E +- %.2E um, gy=%.2E +- %.2E um, hx=%.2E +- %.2E mrad, hy=%.2E +- %.2E mrad\n",
        gx*1E3, gx_error*1E3, gy*1E3, gy_error*1E3, hx*1E3, hx_error*1E3, hy*1E3, hy_error*1E3);
      printf("\trot_z fit\n");
      printf("\t\tconstant: mean = %.2E +- %.2E mrad\n", rot_z_mean*1E3, rot_z_mean_error*1E3);
    }

    // store factored values
    //  sh_r,  sh_r_e,  sh_x,  sh_x_e,  sh_y,  sh_y_e,  sh_z,  sh_z_e,  rot_z,  rot_z_e);
    factored.SetRPCorrection(rpId, RPAlignmentCorrectionData(0., 0., g0x, g0x_error, g0y, g0y_error, 0., 0., rot_z_mean, rot_z_mean_error));

    // calculate and store residuals for sensors
    idx = 0;
    for (auto dit = sensors.begin(); dit != sensors.end(); ++dit, ++idx)
    {
      AlignmentGeometry::const_iterator git = geometry.find(*dit);
      const DetGeometry &d = git->second;
      const RPAlignmentCorrectionData &oa = origAlignments[*dit];

      double s = oa.sh_r() - sc[idx];
      double s_e_full = oa.sh_r_e(); // keep the full error
      double s_e_res = sqrt(VR(idx, idx));

      double zeff = d.z - z0;
      double sp = s - d.dx*(hx*zeff+gx) - d.dy*(hy*zeff+gy);

      double rot_z_res = oa.rot_z() - rot_z_mean;
      double rot_z_e_full = oa.rot_z_e(); // keep the full error
      double rot_z_e_res = sqrt(rot_z_e_full*rot_z_e_full - rot_z_mean_error*rot_z_mean_error);

      if (verbosity > 1)
        printf("\t%u [%u] | sh=%.3f, sh_e_full=%.3f, sh_e_res=%.3f | sh_lin_res=%.3f | rot=%.3f, rot_e_full=%.3f, rot_e_res=%.3f\n",
          *dit, idx,
          s*1E3, s_e_full*1E3, s_e_res*1E3,
          sp,
          rot_z_res*1E3, rot_z_e_full*1E3, rot_z_e_res*1E3);

      RPAlignmentCorrectionData ac(
        s, s_e_full,
        s*d.dx, s_e_full*d.dx, s*d.dy, s_e_full*d.dy,   // sigma(sh_x) = sigma(sh_r) * dx
        oa.sh_z(), oa.sh_z_e(),
        rot_z_res, rot_z_e_full
      );
      factored.SetSensorCorrection(*dit, ac);
    }
  }
#endif
}
