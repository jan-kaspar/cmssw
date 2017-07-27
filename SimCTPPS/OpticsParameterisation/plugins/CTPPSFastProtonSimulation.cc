/****************************************************************************
 *
 * This is a part of CTPPS offline software
 * Authors:
 *   Jan Kašpar
 *   Laurent Forthomme
 *
 ****************************************************************************/

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "CondFormats/CTPPSOpticsObjects/interface/LHCOpticsApproximator.h"

#include "Geometry/VeryForwardRPTopology/interface/RPTopology.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/Records/interface/VeryForwardMisalignedGeometryRecord.h"

#include "CLHEP/Random/RandGauss.h"

#include <unordered_map>
#include <array>

class CTPPSFastProtonSimulation : public edm::stream::EDProducer<>
{
  public:
    explicit CTPPSFastProtonSimulation( const edm::ParameterSet& );
    ~CTPPSFastProtonSimulation();

    static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );

  private:
    struct CTPPSPotInfo {
      CTPPSPotInfo() : detid( 0 ), resolution( 0.0 ), z_position( 0.0 ), approximator( 0 ) {}
      CTPPSPotInfo( const CTPPSDetId& det_id, double resol, double z_position, LHCOpticsApproximator* approx ) :
        detid( det_id ), resolution( resol ), z_position( z_position ), approximator( approx ) {}

      CTPPSDetId detid;
      double resolution;
      double z_position;
      LHCOpticsApproximator* approximator;
    };

    virtual void beginRun( const edm::Run&, const edm::EventSetup& ) override;

    virtual void produce( edm::Event&, const edm::EventSetup& ) override;

    void transportProton(const HepMC::GenVertex* in_vtx, const HepMC::GenParticle* in_trk,
      edm::DetSetVector<TotemRPRecHit>& out_hits) const;

    bool produceRecHit( const CLHEP::Hep3Vector& coord_global, const CTPPSDetId& detid, TotemRPRecHit& rechit ) const;

    edm::EDGetTokenT<edm::HepMCProduct> hepMCToken_;

    edm::ParameterSet beamConditions_;
    double sqrtS_;
    double halfCrossingAngleSector45_, halfCrossingAngleSector56_;
    double yOffsetSector45_, yOffsetSector56_;

    bool roundToPitch_;

    bool checkApertures_;
    bool produceHitsRelativeToBeam_;

    /// strip pitch in mm
    double pitch_;

    /// size of insensitive margin at sensor's edge facing the beam, in mm
    double insensitiveMargin_;

    edm::FileInPath opticsFileBeam1_, opticsFileBeam2_;

    std::vector<edm::ParameterSet> detectorPackages_;

    std::vector<CTPPSPotInfo> pots_;

    /// map: RP id -> vector of sensor ids
    std::unordered_map<unsigned int, std::vector<CTPPSDetId> > strips_list_;

    /// internal variable: v position of strip 0, in mm
    double stripZeroPosition_;

    // TODO: doesn't need to be class member
    edm::ESHandle<TotemRPGeometry> geometry_;

    static const bool invertBeamCoordinatesSystem_;
};

//----------------------------------------------------------------------------------------------------

const bool CTPPSFastProtonSimulation::invertBeamCoordinatesSystem_ = true;

//----------------------------------------------------------------------------------------------------

CTPPSFastProtonSimulation::CTPPSFastProtonSimulation( const edm::ParameterSet& iConfig ) :
  hepMCToken_( consumes<edm::HepMCProduct>( iConfig.getParameter<edm::InputTag>( "beamParticlesTag" ) ) ),
  beamConditions_             ( iConfig.getParameter<edm::ParameterSet>( "beamConditions" ) ),
  sqrtS_                      ( beamConditions_.getParameter<double>( "sqrtS" ) ),
  halfCrossingAngleSector45_  ( beamConditions_.getParameter<double>( "halfCrossingAngleSector45" ) ),
  halfCrossingAngleSector56_  ( beamConditions_.getParameter<double>( "halfCrossingAngleSector56" ) ),
  yOffsetSector45_            ( beamConditions_.getParameter<double>( "yOffsetSector45" ) ),
  yOffsetSector56_            ( beamConditions_.getParameter<double>( "yOffsetSector56" ) ),
  roundToPitch_               ( iConfig.getParameter<bool>( "roundToPitch" ) ),
  checkApertures_             ( iConfig.getParameter<bool>( "checkApertures" ) ),
  produceHitsRelativeToBeam_  ( iConfig.getParameter<bool>( "produceHitsRelativeToBeam" ) ),
  pitch_                      ( iConfig.getParameter<double>( "pitch" ) ),
  insensitiveMargin_          ( iConfig.getParameter<double>( "insensitiveMargin" ) ),
  opticsFileBeam1_            ( iConfig.getParameter<edm::FileInPath>( "opticsFileBeam1" ) ),
  opticsFileBeam2_            ( iConfig.getParameter<edm::FileInPath>( "opticsFileBeam2" ) ),
  detectorPackages_           ( iConfig.getParameter< std::vector<edm::ParameterSet> >( "detectorPackages" ) )
{
  produces< edm::DetSetVector<TotemRPRecHit> >();

  // v position of strip 0
  stripZeroPosition_ = RPTopology::last_strip_to_border_dist_ + (RPTopology::no_of_strips_-1)*RPTopology::pitch_ - RPTopology::y_width_/2.;

  auto f_in_optics_beam1 = std::make_unique<TFile>( opticsFileBeam1_.fullPath().c_str() ),
       f_in_optics_beam2 = std::make_unique<TFile>( opticsFileBeam2_.fullPath().c_str() );

  // load optics and interpolators
  for ( const auto& rp : detectorPackages_ )
  {
    const std::string interp_name = rp.getParameter<std::string>( "interpolatorName" );
    const unsigned int raw_detid = rp.getParameter<unsigned int>( "potId" );
    const double det_resol = rp.getParameter<double>( "resolution" );
    const double z_position = rp.getParameter<double>( "zPosition" );
    CTPPSDetId detid( raw_detid );

    if ( detid.arm()==0 ) // sector 45 -- beam 2
      pots_.emplace_back( detid, det_resol, z_position, dynamic_cast<LHCOpticsApproximator*>( f_in_optics_beam2->Get( interp_name.c_str() ) ) );
    if ( detid.arm()==1 ) // sector 56 -- beam 1
      pots_.emplace_back( detid, det_resol, z_position, dynamic_cast<LHCOpticsApproximator*>( f_in_optics_beam1->Get( interp_name.c_str() ) ) );
  }
}

//----------------------------------------------------------------------------------------------------

CTPPSFastProtonSimulation::~CTPPSFastProtonSimulation()
{}


//----------------------------------------------------------------------------------------------------

void CTPPSFastProtonSimulation::beginRun( const edm::Run&, const edm::EventSetup& iSetup )
{
  // get geometry
  iSetup.get<VeryForwardMisalignedGeometryRecord>().get( geometry_ );

  std::ostringstream os;
  for ( const auto& rp : pots_ )
  {
    std::vector<CTPPSDetId>& list = strips_list_[rp.detid.rawId()];
    os << "\npot " << rp.detid << ":";
    for ( TotemRPGeometry::mapType::const_iterator it=geometry_->beginDet(); it!=geometry_->endDet(); ++it )
    {
      const CTPPSDetId detid( it->first );
      if ( detid.getRPId()!=rp.detid )
        continue;
      list.push_back( detid );
      os << "\n* " << TotemRPDetId( detid );
    }
  }

  edm::LogWarning("CTPPSFastProtonSimulation::beginRun") << "Hierarchy of DetIds for each pot" << os.str();
}

//----------------------------------------------------------------------------------------------------

void
CTPPSFastProtonSimulation::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  // get input
  edm::Handle<edm::HepMCProduct> hepmc_prod;
  iEvent.getByToken( hepMCToken_, hepmc_prod );

  // prepare outputs
  std::unique_ptr< edm::DetSetVector<TotemRPRecHit> > pRecHits( new edm::DetSetVector<TotemRPRecHit>() );

  // loop over event vertices
  auto evt = new HepMC::GenEvent( *hepmc_prod->GetEvent() );
  for ( auto it_vtx = evt->vertices_begin(); it_vtx != evt->vertices_end(); ++it_vtx )
  {
    auto vtx = *( it_vtx );

    // loop over outgoing particles
    for ( auto it_part = vtx->particles_out_const_begin(); it_part != vtx->particles_out_const_end(); ++it_part )
    {
      auto part = *( it_part );

      // accept only stable protons
      if ( part->pdg_id()!=2212 )
        continue;
      if ( part->status()!=1 && part->status()<83 )
        continue;

      transportProton(vtx, part, *pRecHits);
    }
  }

  iEvent.put( std::move( pRecHits ) );
}

//----------------------------------------------------------------------------------------------------

void CTPPSFastProtonSimulation::transportProton(const HepMC::GenVertex* in_vtx, const HepMC::GenParticle* in_trk,
  edm::DetSetVector<TotemRPRecHit>& out_hits) const
{
  /// xi is positive for diffractive protons, thus proton momentum p = (1-xi) * p_nom
  /// horizontal component of proton momentum: p_x = th_x * (1-xi) * p_nom
  
  // TODO: the code below uses the TOTEM/LHC convention for momentum coordinates. It is to be changed to the CMS one.

  const HepMC::FourVector vtx = in_vtx->position();
  const HepMC::FourVector mom = in_trk->momentum();

  // determine the LHC arm and related parameters
  unsigned int arm = 3;
  double half_cr_angle = 0.0, vtx_y_offset = 0.0;
  double z_sign;

  if (mom.z() < 0)  // sector 45
  {
    arm = 0;
    z_sign = -1;
    vtx_y_offset = yOffsetSector45_;
    half_cr_angle = halfCrossingAngleSector45_;
  } else {  // sector 56
    arm = 1;
    z_sign = +1;
    vtx_y_offset = yOffsetSector56_;
    half_cr_angle = halfCrossingAngleSector56_;
  }
  
  // transport the proton into each pot
  for ( const auto& rp : pots_ )
  {
    // first check the arm
    if ( rp.detid.arm() != arm)
      continue;

    // so far only works for strips
    if ( rp.detid.subdetId() != CTPPSDetId::sdTrackingStrip )
      continue;

    // calculate kinematics for optics parametrisation
    const double p0 = rp.approximator->GetBeamMomentum();
    const double p = mom.rho();
    const double xi = 1. - p / p0;
    const double th_x_phys = mom.x() / p;
    const double th_y_phys = mom.y() / p;

    // transport proton
    double kin_tr_in[5] = { vtx.x()*1E-3, (th_x_phys + half_cr_angle) * (1.-xi), vtx.y()*1E-3 + vtx_y_offset, th_y_phys * (1.-xi), -xi };
    double kin_tr_out[5];
    bool proton_transported = rp.approximator->Transport(kin_tr_in, kin_tr_out, checkApertures_, invertBeamCoordinatesSystem_);
    const double b_x_tr = kin_tr_out[0], b_y_tr = kin_tr_out[2];
    const double a_x_tr = kin_tr_out[1]/(1.-xi), a_y_tr = kin_tr_out[3]/(1.-xi);

    // determine beam position
    double kin_be_in[5] = { 0., half_cr_angle, vtx_y_offset, 0., 0. };
    double kin_be_out[5];
    rp.approximator->Transport(kin_be_in, kin_be_out, false, invertBeamCoordinatesSystem_);
    const double b_x_be = kin_be_out[0], b_y_be = kin_be_out[2];
    const double a_x_be = kin_be_out[1], a_y_be = kin_be_out[3];

    // stop if proton not transported
    if (!proton_transported)
      continue;

    // get z position of the approximator scoring plane
    const double approximator_z = rp.z_position * 1E3; // in mm

    // TODO: replace by a method from TotemRPGeometry ?
    // retrieve the sensor from geometry
    for ( const auto& detid : strips_list_.at( rp.detid.rawId() ) )
    {
      // get sensor geometry
      const double gl_o_z = geometry_->LocalToGlobal( detid, CLHEP::Hep3Vector() ).z(); // in mm

      // evaluate positions (in mm) of track and beam
      const double de_z = (gl_o_z - approximator_z) * z_sign;

      const double x_tr = a_x_tr * de_z + b_x_tr * 1E3;
      const double y_tr = a_y_tr * de_z + b_y_tr * 1E3;

      const double x_be = a_x_be * de_z + b_x_be * 1E3;
      const double y_be = a_y_be * de_z + b_y_be * 1E3;

      // determine hit global coordinates
      CLHEP::Hep3Vector h_glo(x_tr, y_tr, gl_o_z);
      if (produceHitsRelativeToBeam_)
          h_glo -= CLHEP::Hep3Vector(x_be, y_be, 0.);

      // make hit (if within acceptance)
      TotemRPRecHit hit; // all coordinates in mm
      if ( produceRecHit( h_glo, detid, hit ) )
      {
        edm::DetSet<TotemRPRecHit>& hits = out_hits.find_or_insert( detid );
        hits.push_back( hit );
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------

bool CTPPSFastProtonSimulation::produceRecHit( const CLHEP::Hep3Vector& coord_global, const CTPPSDetId& detid, TotemRPRecHit& rechit ) const
{
  // transform hit global to local coordinates
  const CLHEP::Hep3Vector h_loc = geometry_->GlobalToLocal( detid, coord_global );

  double u = h_loc.x();
  double v = h_loc.y();

  // is it within detector?
  if ( !RPTopology::IsHit( u, v, insensitiveMargin_ ) )
    return false;

  // round the measurement
  if ( roundToPitch_ )
  {
    double m = stripZeroPosition_ - v;
    int strip = static_cast<int>( floor( m/pitch_ + 0.5 ) );
    v = stripZeroPosition_ - pitch_ * strip;
  }

  const double sigma = pitch_ / sqrt( 12. );

  rechit = TotemRPRecHit( v, sigma );

  return true;
}

//----------------------------------------------------------------------------------------------------

void
CTPPSFastProtonSimulation::fillDescriptions( edm::ConfigurationDescriptions& descriptions )
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault( desc );
}

// define this as a plug-in
DEFINE_FWK_MODULE( CTPPSFastProtonSimulation );
