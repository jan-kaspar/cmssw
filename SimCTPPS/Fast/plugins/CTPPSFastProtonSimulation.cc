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
    printf("\n>> CTPPSFastProtonSimulation::produce > event %llu\n", event.id().event());

  // get input
  Handle<HepMCProduct> hepMCProduct;
  event.getByToken(tokenHepMC, hepMCProduct);

  // get geometry
  printf("dupa1\n");
  ESHandle<TotemRPGeometry> geometry;
  es.get<VeryForwardMisalignedGeometryRecord>().get(geometry);
  printf("dupa2\n");

  // prepare output
  unique_ptr<DetSetVector<TotemRPRecHit>> stripHitColl(new DetSetVector<TotemRPRecHit>());

  // loop over vertices
  const GenEvent *mcEvent = hepMCProduct->GetEvent();
  for (GenEvent::vertex_const_iterator vit = mcEvent->vertices_begin(); vit != mcEvent->vertices_end(); ++vit)
  {
    const FourVector &vertex = (*vit)->position(); // in mm

    printf("vertex: x=%f, y=%f, z=%f\n", vertex.x(), vertex.y(), vertex.z());
    
    // loop over outgoing particles
    for (GenVertex::particles_out_const_iterator pit = (*vit)->particles_out_const_begin(); pit != (*vit)->particles_out_const_end(); ++pit)
    {
      // keep only protons
      if ((*pit)->pdg_id() != 2212)
        continue;

      const FourVector &momentum = (*pit)->momentum();

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
      const double vtx0_y = (arm == 1) ? vtx0_y_56 : vtx0_y_45;
      const double half_crossing_angle = (arm == 1) ? half_crossing_angle_56 : half_crossing_angle_45;

      const double p0 = opticsApproximator->GetBeamMomentum();
      const double p = momentum.rho();
      const double xi = 1. - p / p0;
      const double th_x_phys = momentum.x() / p;
      const double th_y_phys = momentum.y() / p;

      double kin_tr_in[5] = { vertex.x(), (th_x_phys + half_crossing_angle) * (1.-xi), vertex.y() + vtx0_y, th_y_phys * (1.-xi), -xi };
      double kin_be_in[5] = { 0., half_crossing_angle, vtx0_y, 0., 0. };

      const bool check_appertures = false;
      const bool invert_beam_coord_sytems = true;

      double kin_tr_out[5];
      double kin_be_out[5];

      bool proton_transported = opticsApproximator->Transport(kin_tr_in, kin_tr_out, check_appertures, invert_beam_coord_sytems);
      opticsApproximator->Transport(kin_be_in, kin_be_out, check_appertures, invert_beam_coord_sytems);

      // stop if proton not transported
      if (!proton_transported)
        continue;

      // check all sensors known to geometry
      for (TotemRPGeometry::mapType::const_iterator it = geometry->beginDet(); it != geometry->endDet(); ++it)
      {
        // so far only works for strips
        CTPPSDetId detId(it->first);
        if (detId.subdetId() != CTPPSDetId::sdTrackingPixel)
          continue;

        // get geometry
        Hep3Vector gl_o = geometry->LocalToGlobal(detId, Hep3Vector(0, 0, 0));

        cout << TotemRPDetId(it->first) << "z = " << gl_o.z() << endl;
      }
    
#if 0
      // loop over all stations of the arm
      const set<unsigned int> stations = geom->StationsInArm(arm);
      for (set<unsigned int>::const_iterator st = stations.begin(); st != stations.end(); ++st)
      {
        //printf("station %u\n", *st);
        // beam misalignment for the station
        map<unsigned int, BeamMisalignment>::iterator bm = bms.find(*st);
  
        const set<unsigned int> RPs = geom->RPsInStation(*st);
        for (set<unsigned int>::const_iterator rp = RPs.begin(); rp != RPs.end(); ++rp)
        {
            //printf("RP %u\n", *rp);
            // center of the RP box, in mm
            double z0 = geom->GetRPDevice(*rp)->translation().z();
        
            // optical functions for the RP
            LHCOpticsApproximator *of = optFun->GetFunction(*rp);
            if (!of)
              throw cms::Exception("CTPPSFastProtonSimulation::produce") << "No optical function found for RP " << *rp;
        
            // get position and direction of the proton at the center of the RP box
            double xi = p.e() / optPar->GetBeamEnergy() - 1.;
            double parIn[5] = {vertex.x() * 1E-3, p.x()/p.rho(), vertex.y() * 1E-3, p.y()/p.rho(), xi};
            double parOut[5]; // in m, rad and xi dimensionless
            bool transportable = of->Transport(parIn, parOut, true);
        
            if (verbosity > 5)
              printf("RP%u\tz0 = %E\t p.z = %E\t transportable = %i, beam type = %i\n", *rp, z0, p.z(), transportable, of->GetBeamType());
  
            if (verbosity > 8)
            {
              printf("param input:  %E\t%E\t%E\t%E\t%E\n", parIn[0], parIn[1], parIn[2], parIn[3], parIn[4]);
              printf("param output: %E\t%E\t%E\t%E\t%E\n", parOut[0], parOut[1], parOut[2], parOut[3], parOut[4]);
            }
        
            // do not continue if the proton is stopped
            if (!transportable)
              continue;
        
            // collection of hits
            vector<RPRecoHit> hits;
  
            // u and v hit counter
            unsigned uHits = 0, vHits = 0;
        
            // loop over all detectors within the RP
            set<unsigned int> dets = geom->DetsInRP(*rp);
            for (set<unsigned int>::iterator dit = dets.begin(); dit != dets.end(); ++dit)
            {
              // calculate hit position in global coordinates
              unsigned int rawId = TotRPDetId::DecToRawId(*dit);
              double z = geom->GetDetector(rawId)->translation().z(); // in mm
              double x = parOut[0]*1E3  + parOut[1] * (z - z0); // in mm
              double y = parOut[2]*1E3  + parOut[3] * (z - z0); // in mm
  
              // apply beam misalignment
              if (bm != bms.end())
              {
                double dx = 0., dy = 0.; // in mm
                bm->second.CalculateShift(z, dx, dy);
                x += dx;
                y += dy;
              }
        
              if (verbosity > 5)
                printf("\t%u\t%u\t%E\t%e\t%E\n", *dit, rawId, z, x, y);
        
              // convert to local coordinates
              CLHEP::Hep3Vector hit(x, y, z);
              hit = geom->GlobalToLocal(rawId, hit);
              double u = hit.x(); // in mm
              double v = hit.y(); // in mm
        
              // is it within detector?
              if (!RPTopology::IsHit(u, v, insensitiveMargin))
                continue;
  
              if (TotRPDetId::IsStripsCoordinateUDirection(*dit))
                uHits++;
              else
                vHits++;
        
              // round the measurement
              if (roundToPitch)
              {
                double m = stripZeroPosition - v;
                signed int strip = (int) floor(m / pitch);
                double offset = m - strip*pitch;
                double dstrip = 0.;
        
                if (offset < dscrStart)
                {
                  dstrip = 0.;
                } else {
                  if (offset < dscrEnd)
                    dstrip = 0.5;
                  else
                    dstrip = 1.;
                }
        
                v = stripZeroPosition - pitch * (strip + dstrip);
                if (verbosity > 5)
                  printf(" | stripZeroPosition = %+8.4f, strip = %+6.1f", stripZeroPosition, strip+dstrip);
              }
        
              if (verbosity > 5)
                printf("\t\tv = %E\n", v);
        
              hits.push_back(RPRecoHit(rawId, v, pitch / sqrt(12.)));
            }
        
            if (verbosity > 5)
              printf("\tRP %i has %lu hits\n", *rp, hits.size());
        
            // insert to output collection
            bool fittable = (uHits > minUVSensorsPerRP && vHits > minUVSensorsPerRP);
            (*output)[*rp] = RPTrackCandidate(hits, fittable);
        }
      }
#endif

    }
  }

  // save output
  event.put(move(stripHitColl));
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(CTPPSFastProtonSimulation);
