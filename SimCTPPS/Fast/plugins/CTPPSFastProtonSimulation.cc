/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"

#include "Geometry/VeryForwardRPTopology/interface/RPTopology.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/Records/interface/VeryForwardMisalignedGeometryRecord.h"

#include "RecoCTPPS/OpticsParametrization/interface/LHCOpticsApproximator.h"

//----------------------------------------------------------------------------------------------------

class CTPPSFastProtonSimulation : public edm::EDProducer
{
  public:
    CTPPSFastProtonSimulation(const edm::ParameterSet &);
    virtual ~CTPPSFastProtonSimulation();

  protected:
    unsigned int verbosity;

    edm::EDGetTokenT<edm::HepMCProduct> tokenHepMC;

    /// angular limit to distinguish between forward, central and backward particles
    double thetaLimit;

    /// whether the beam position should be subracted from hit positions
    bool produceHitsRelativeToBeam;
    
    /// whether measurement values shall be rounded to the nearest strip
    bool roundToPitch;

    /// strip pitch in mm
    double pitch;

    /// size of insensitive margin at sensor's edge facing the beam, in mm
    double insensitiveMargin;
    
    /// internal variable: v position of strip 0, in mm
    double stripZeroPosition;

    LHCOpticsApproximator *opticsApproximator_45, *opticsApproximator_56;
    double opticsZ0_45, opticsZ0_56;

    double vtx0_y_45, vtx0_y_56;
    double half_crossing_angle_45, half_crossing_angle_56;

    std::map<unsigned int, LHCOpticsApproximator *> mRPOptics;

    virtual void produce(edm::Event &, const edm::EventSetup&) override;
};

//----------------------------------------------------------------------------------------------------

using namespace edm;
using namespace std;
using namespace HepMC;
using namespace CLHEP;

//----------------------------------------------------------------------------------------------------

CTPPSFastProtonSimulation::CTPPSFastProtonSimulation(const ParameterSet &ps) :
    verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),

    tokenHepMC(consumes<HepMCProduct>(ps.getParameter<InputTag>("tagHepMC"))),

    thetaLimit(ps.getParameter<double>("thetaLimit")),
    produceHitsRelativeToBeam(ps.getParameter<bool>("produceHitsRelativeToBeam")),
    roundToPitch(ps.getParameter<bool>("roundToPitch")),
    pitch(ps.getParameter<double>("pitch")),

    insensitiveMargin(ps.getParameter<double>("insensitiveMargin"))
{
  // v position of strip 0
  stripZeroPosition = RPTopology::last_strip_to_border_dist_ + (RPTopology::no_of_strips_-1)*RPTopology::pitch_
    - RPTopology::y_width_/2.;

  produces<DetSetVector<TotemRPRecHit>>();

  // load optics and beam parameters
  string opticsFile_45 = ps.getParameter<string>("opticsFile_45");
  string opticsObject_45 = ps.getParameter<string>("opticsObject_45");
  opticsZ0_45 = ps.getParameter<double>("opticsZ0_45");
  vtx0_y_45 = ps.getParameter<double>("vtx0_y_45");
  half_crossing_angle_45 = ps.getParameter<double>("half_crossing_angle_45");

  string opticsFile_56 = ps.getParameter<string>("opticsFile_56");
  string opticsObject_56 = ps.getParameter<string>("opticsObject_56");
  opticsZ0_56 = ps.getParameter<double>("opticsZ0_56");
  vtx0_y_56 = ps.getParameter<double>("vtx0_y_56");
  half_crossing_angle_56 = ps.getParameter<double>("half_crossing_angle_56");

  TFile *f_in_45 = TFile::Open(opticsFile_45.c_str());
  if (!f_in_45)
    throw cms::Exception("CTPPSFastProtonSimulation") << "Can't open file '" + opticsFile_45 + "'.";

  opticsApproximator_45 = (LHCOpticsApproximator *) f_in_45->Get(opticsObject_45.c_str());
  if (!opticsApproximator_45)
    throw cms::Exception("CTPPSFastProtonSimulation") << "Can't load object '" << opticsObject_45 << "'.";

  opticsApproximator_45 = new LHCOpticsApproximator(*opticsApproximator_45);
  
  TFile *f_in_56 = TFile::Open(opticsFile_56.c_str());
  if (!f_in_56)
    throw cms::Exception("CTPPSFastProtonSimulation") << "Can't open file '" + opticsFile_56 + "'.";

  opticsApproximator_56 = (LHCOpticsApproximator *) f_in_56->Get(opticsObject_56.c_str());
  if (!opticsApproximator_56)
    throw cms::Exception("CTPPSFastProtonSimulation") << "Can't load object '" << opticsObject_56 << "'.";

  opticsApproximator_56 = new LHCOpticsApproximator(*opticsApproximator_56);

  // TODO: debug
  {
    //mRPOptics[0] = (LHCOpticsApproximator *) f_in_45->Get("ip5_to_station_150_v_1_lhcb2");
    mRPOptics[2] = (LHCOpticsApproximator *) f_in_45->Get("ip5_to_station_150_h_1_lhcb2");
    mRPOptics[3] = (LHCOpticsApproximator *) f_in_45->Get("ip5_to_station_150_h_2_lhcb2");
    //mRPOptics[5] = (LHCOpticsApproximator *) f_in_45->Get("ip5_to_station_150_v_2_lhcb2");
    //mRPOptics[100] = (LHCOpticsApproximator *) f_in_56->Get("ip5_to_station_150_v_1_lhcb1");
    mRPOptics[102] = (LHCOpticsApproximator *) f_in_56->Get("ip5_to_station_150_h_1_lhcb1");
    mRPOptics[103] = (LHCOpticsApproximator *) f_in_56->Get("ip5_to_station_150_h_2_lhcb1");
    //mRPOptics[105] = (LHCOpticsApproximator *) f_in_56->Get("ip5_to_station_150_v_2_lhcb1");
  }

  // TODO: uncomment
  delete f_in_45;
  delete f_in_56;
}

//----------------------------------------------------------------------------------------------------

CTPPSFastProtonSimulation::~CTPPSFastProtonSimulation()
{
}

//----------------------------------------------------------------------------------------------------

void CTPPSFastProtonSimulation::produce(edm::Event &event, const EventSetup &es)
{
  if (verbosity > 5)
  {
    printf("\n--------------------------------------------------------------------------------------\n");
    printf(">> CTPPSFastProtonSimulation::produce > event %llu\n", event.id().event());
  }

  // get input
  Handle<HepMCProduct> hepMCProduct;
  event.getByToken(tokenHepMC, hepMCProduct);

  // get geometry
  ESHandle<TotemRPGeometry> geometry;
  es.get<VeryForwardMisalignedGeometryRecord>().get(geometry);

  // prepare output
  unique_ptr<DetSetVector<TotemRPRecHit>> stripHitColl(new DetSetVector<TotemRPRecHit>());

  // loop over vertices
  const GenEvent *mcEvent = hepMCProduct->GetEvent();
  for (GenEvent::vertex_const_iterator vit = mcEvent->vertices_begin(); vit != mcEvent->vertices_end(); ++vit)
  {
    // TODO: revert
    //const FourVector &vertex = (*vit)->position(); // in mm
    FourVector vertex = FourVector(0E-3, 100E-3, 0., 0.); // in mm

    printf("vertex: x=%f, y=%f, z=%f\n", vertex.x(), vertex.y(), vertex.z());
    
    // loop over outgoing particles
    for (GenVertex::particles_out_const_iterator pit = (*vit)->particles_out_const_begin(); pit != (*vit)->particles_out_const_end(); ++pit)
    {
      // keep only protons
      if ((*pit)->pdg_id() != 2212)
        continue;

      // TODO: revert
      //const FourVector &momentum = (*pit)->momentum();
      FourVector momentum = (*pit)->momentum();

      // TODO
      {
        double th_x = 0E-6;
        double th_y = 0E-6;
        double xi = 0.;
        double p0 = 6500.;

        double p = (1.-xi) * p0;
        double px = th_x * p;
        double py = th_y * p;
        double pz = (momentum.z() > 0) ? p : -p;

        momentum = FourVector(px, py, pz, p);
      }

      printf("momentum: px=%f, py=%f, pz=%f, E=%f\n", momentum.x(), momentum.y(), momentum.z(), momentum.e());

      // determine the arm of action 
      double th = momentum.theta();
      unsigned int arm = 123;
      if (th < thetaLimit) arm = 1;
      if (th > (M_PI - thetaLimit)) arm = 0;

      // skip central particles
      if (arm > 1)
        continue;

      printf("    arm = %u\n", arm);

      // transport proton
      LHCOpticsApproximator *opticsApproximator = (arm == 1) ? opticsApproximator_56 : opticsApproximator_45;
      const double opticsZ0 = (arm == 1) ? opticsZ0_56 : opticsZ0_45;
      const double z_sign = (arm == 1) ? +1. : -1.;
      const double vtx0_y = (arm == 1) ? vtx0_y_56 : vtx0_y_45;
      const double half_crossing_angle = (arm == 1) ? half_crossing_angle_56 : half_crossing_angle_45;

      const double p0 = opticsApproximator->GetBeamMomentum();
      const double p = momentum.rho();
      const double xi = 1. - p / p0;
      const double th_x_phys = momentum.x() / p;
      const double th_y_phys = momentum.y() / p;

      const bool check_appertures = false;
      const bool invert_beam_coord_sytems = true;

      double kin_tr_in[5] = { vertex.x()*1E-3, (th_x_phys + half_crossing_angle) * (1.-xi), vertex.y()*1E-3 + vtx0_y, th_y_phys * (1.-xi), -xi };
      double kin_tr_out[5];
      bool proton_transported = opticsApproximator->Transport(kin_tr_in, kin_tr_out, check_appertures, invert_beam_coord_sytems);
      const double b_x_tr = kin_tr_out[0], b_y_tr = kin_tr_out[2];
      const double a_x_tr = kin_tr_out[1]/(1.-xi), a_y_tr = kin_tr_out[3]/(1.-xi);

      double kin_be_in[5] = { 0., half_crossing_angle, vtx0_y, 0., 0. };
      double kin_be_out[5];
      opticsApproximator->Transport(kin_be_in, kin_be_out, check_appertures, invert_beam_coord_sytems);
      const double b_x_be = kin_be_out[0], b_y_be = kin_be_out[2];
      const double a_x_be = kin_be_out[1], a_y_be = kin_be_out[3];

      printf("    track: ax=%f, bx=%f, ay=%f, by=%f\n", a_x_tr, b_x_tr, a_y_tr, b_y_tr);
      printf("    beam: ax=%f, bx=%f, ay=%f, by=%f\n", a_x_be, b_x_be, a_y_be, b_y_be);

      // TODO: remove
      {
        for (const auto &it : mRPOptics)
        {
          if ((it.first / 100) != arm)
            continue;

          double out[5];
          it.second->Transport(kin_tr_in, out, check_appertures, invert_beam_coord_sytems);
          printf("    RP: %u, x = %f, y = %f\n", it.first, out[0]*1e3, out[2]*1e3);
        }
      }

      // stop if proton not transported
      if (!proton_transported)
        continue;

      // check all sensors known to geometry
      for (TotemRPGeometry::mapType::const_iterator it = geometry->beginDet(); it != geometry->endDet(); ++it)
      {
        // so far only works for strips
        CTPPSDetId detId(it->first);
        if (detId.subdetId() != CTPPSDetId::sdTrackingStrip)
          continue;

        // is in the correct arm?
        if (detId.arm() != arm)
          continue;

        // get geometry
        Hep3Vector gl_o = geometry->LocalToGlobal(detId, Hep3Vector(0, 0, 0));

        // evaluate positions (in mm) of track and beam
        const double de_z = (gl_o.z() - opticsZ0) * z_sign;

        const double x_tr = a_x_tr * de_z + b_x_tr * 1E3;
        const double y_tr = a_y_tr * de_z + b_y_tr * 1E3;

        const double x_be = a_x_be * de_z + b_x_be * 1E3;
        const double y_be = a_y_be * de_z + b_y_be * 1E3;

        /*
        cout << TotemRPDetId(it->first) << ", z = " << gl_o.z() << ", de z = " << (gl_o.z() - opticsZ0) <<
          " | track: x=" << x_tr << ", y=" << y_tr <<
          " | beam: x=" << x_be << ", y=" << y_be <<
          endl;
        */

        // global hit in coordinates "aligned to beam" (as in the RP alignment)
        Hep3Vector h_glo = (produceHitsRelativeToBeam) ? Hep3Vector(x_tr - x_be, y_tr - y_be, gl_o.z()) :
          Hep3Vector(x_tr, y_tr, gl_o.z());

        // transform hit global to local coordinates
        Hep3Vector h_loc = geometry->GlobalToLocal(detId, h_glo);

        double u = h_loc.x();
        double v = h_loc.y();

        // is it within detector?
        if (!RPTopology::IsHit(u, v, insensitiveMargin))
          continue;

        // round the measurement
        if (roundToPitch)
        {
          double m = stripZeroPosition - v;
          signed int strip = (int) floor(m / pitch + 0.5);

          v = stripZeroPosition - pitch * strip;
        }

        const double sigma = pitch / sqrt(12.);

        DetSet<TotemRPRecHit> &hits = stripHitColl->find_or_insert(detId);
        hits.push_back(TotemRPRecHit(v, sigma));
      }
    }
  }

  // save output
  event.put(move(stripHitColl));
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(CTPPSFastProtonSimulation);
