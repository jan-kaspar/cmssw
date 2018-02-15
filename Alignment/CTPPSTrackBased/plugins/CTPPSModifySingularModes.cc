/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/CTPPSAlignment/interface/RPAlignmentCorrectionData.h"
#include "DataFormats/CTPPSAlignment/interface/RPAlignmentCorrectionsData.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/RPAlignmentCorrectionsMethods.h"

#include "Alignment/CTPPSTrackBased/interface/AlignmentTask.h"
#include "Alignment/CTPPSTrackBased/interface/CommonMethods.h"

/**
 *\brief Modifies the alignment modes unconstrained by the track-based alignment.
 **/
class CTPPSModifySingularModes : public edm::EDAnalyzer
{
  public:
    CTPPSModifySingularModes(const edm::ParameterSet &ps); 

  private:
    edm::ParameterSet ps;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

    virtual void analyze(const edm::Event &e, const edm::EventSetup &es) override
    {
    }
};

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

CTPPSModifySingularModes::CTPPSModifySingularModes(const ParameterSet &_ps) : ps(_ps)
{
}

//----------------------------------------------------------------------------------------------------

void CTPPSModifySingularModes::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  // get config parameters
  double z1 = ps.getUntrackedParameter<double>("z1");
  double z2 = ps.getUntrackedParameter<double>("z2");
  double de_x1 = ps.getUntrackedParameter<double>("de_x1");
  double de_x2 = ps.getUntrackedParameter<double>("de_x2");
  double de_y1 = ps.getUntrackedParameter<double>("de_y1");
  double de_y2 = ps.getUntrackedParameter<double>("de_y2");
  double de_rho1 = ps.getUntrackedParameter<double>("de_rho1");
  double de_rho2 = ps.getUntrackedParameter<double>("de_rho2");

  FileInPath inputFileInPath(ps.getUntrackedParameter<string>("inputFile"));
  string inputFile = inputFileInPath.fullPath();
  string outputFile = ps.getUntrackedParameter<string>("outputFile");

  // validate config parameters
  if (z1 == z2)
    throw cms::Exception("CTPPSModifySingularModes") << "z1 equals z2";

  // temporarily
  if (de_rho1 != 0. || de_rho2 != 0.)
    throw cms::Exception("CTPPSModifySingularModes") << "rotations are not yet supported";

  // calculate slopes and intercepts
  double a_x = (de_x2 - de_x1) / (z2 - z1), b_x = de_x1 - a_x * z1;
  double a_y = (de_y2 - de_y1) / (z2 - z1), b_y = de_y1 - a_y * z1;
  //double a_rho = (de_rho2 - de_rho1) / (z2 - z1), b_rho = de_rho1 - a_rho * z1;

  // get geometry
  ESHandle<CTPPSGeometry> geometry;
  es.get<VeryForwardRealGeometryRecord>().get(geometry);

  // get input alignments
  RPAlignmentCorrectionsDataSequence inputSequence = RPAlignmentCorrectionsMethods::loadFromXML(inputFile);
  const auto &input = inputSequence.begin()->second;

  // modify the singular modes
  RPAlignmentCorrectionsData output = input;

  for (auto &it : input.getSensorMap())
  {
    const auto &sensorId = it.first;

    const double z = geometry->getSensorTranslation(sensorId).z(); // TODO: this is crude approximation for pixels

    RPAlignmentCorrectionData d = it.second;
    d.setShX(d.getShX() + a_x*z + b_x);
    d.setShY(d.getShY() + a_y*z + b_y);

    output.setSensorCorrection(sensorId, d);
  }

  // build list of RPs
  vector<unsigned int> rps;
  unsigned int last_rp = 123456;
  for (auto &it : input.getSensorMap())
  {
      CTPPSDetId senId(it.first);
      unsigned int rpDecId = senId.arm()*100 + senId.station()*10 + senId.rp();

      if (last_rp != rpDecId)
      {
        rps.push_back(rpDecId);
        last_rp = rpDecId;
      }
  }

  // build alignment geometry
  AlignmentGeometry alignmentGeometry;
  vector<unsigned int> excludePlanes;
  AlignmentTask::BuildGeometry(rps, excludePlanes, geometry.product(), 0., alignmentGeometry);

  // factorise output
  RPAlignmentCorrectionsData outputExpanded;
  RPAlignmentCorrectionsData outputFactored;

  const bool equalWeights = false;
  CommonMethods::FactorRPFromSensorCorrections(output, outputExpanded, outputFactored,
    alignmentGeometry, equalWeights, 1);

  // save output
  RPAlignmentCorrectionsMethods::writeToXML(outputFactored, outputFile, false, false,
    true, true, true, true);

// TODO: use or remove
#if 0

  // prepare output - expand input to sensor level
  RPAlignmentCorrectionsData output;
  for (RPAlignmentCorrectionsData::mapType::const_iterator it = input.GetSensorMap().begin();
      it != input.GetSensorMap().end(); ++it)
  { 
    // TODO
    //unsigned int rawId = TotRPDetId::DecToRawId(it->first);
    unsigned int rawId = it->first;
    CLHEP::Hep3Vector d = geom->localToGlobalDirection(rawId, CLHEP::Hep3Vector(0., 1., 0.));

    RPAlignmentCorrectionData ac = input.GetFullSensorCorrection(it->first);
    ac.xyTranslationToReadout(d.x(), d.y());
    output.SetSensorCorrection(it->first, ac);
  }

  // apply singular-mode change
  printf("\tID      shift in x    shift in y    rotation about z\n");
  for (RPAlignmentCorrectionsData::mapType::const_iterator it = output.GetSensorMap().begin();
      it != output.GetSensorMap().end(); ++it)
  { 
    // TODO
    //unsigned int rawId = TotRPDetId::DecToRawId(it->first);
    unsigned int rawId = it->first;
    CLHEP::Hep3Vector d = geom->localToGlobalDirection(rawId, CLHEP::Hep3Vector(0., 1., 0.));
    double dx = d.x(), dy = d.y();
    CLHEP::Hep3Vector c = geom->GetDetTranslation(rawId);
    double cx = c.x(), cy = c.y(), z = c.z();

    double de_x = a_x * z + b_x;
    double de_y = a_y * z + b_y;
    double de_rho = a_rho * z + b_rho;

    printf("\t%u %+10.1f um %+10.1f um %+10.1f mrad\n", it->first, de_x*1E3, de_y*1E3, de_rho*1E3);
    //printf("\t\tcx=%e, cy=%E | dx=%E, dy=%E\n", cx, cy, dx, dy);

    double inc_s = +(dx*de_x + dy*de_y) - de_rho * (-dy*(cx + de_x) + dx*(cy + de_y));
    double inc_rho = de_rho;
    //printf("\t\t %E, %E\n", inc_s, inc_rho);

    RPAlignmentCorrectionData &ac = output.GetSensorCorrection(it->first);
    ac.setTranslationR1(ac.sh_r1() + inc_s, ac.sh_r1_e());
    // TODO: what about translation 2
    ac.setRotationZ(ac.rot_z() + inc_rho, ac.rot_z_e());
    ac.readoutTranslationToXY(dx, dy);
  }

  // factorize alignments and write output
  vector<unsigned int> rps;
  unsigned int last_rp = 123456;
  for (RPAlignmentCorrectionsData::mapType::const_iterator it = input.GetSensorMap().begin();
      it != input.GetSensorMap().end(); ++it)
  {
      unsigned int rp = it->first/10;
      if (last_rp != rp)
      {
        rps.push_back(rp);
        last_rp = rp;
      }
  }

  AlignmentGeometry alGeom;
  vector<unsigned int> excludePlanes;
  AlignmentTask::BuildGeometry(rps, excludePlanes, geom.product(), 0., alGeom);

  // TODO
  /*
  RPAlignmentCorrectionsData expanded, factored;
  output.FactorRPFromSensorCorrections(expanded, factored, alGeom);
  factored.WriteXMLFile(ps.getUntrackedParameter<string>("outputFile"));
  */
#endif
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(CTPPSModifySingularModes);
