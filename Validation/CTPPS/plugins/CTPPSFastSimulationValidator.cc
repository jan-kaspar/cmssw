/****************************************************************************
 *
 * This is a part of CTPPS validation software
 * Authors:
 *   Jan Kašpar
 *   Laurent Forthomme
 *
 ****************************************************************************/


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include <map>

//----------------------------------------------------------------------------------------------------

class CTPPSFastSimulationValidator : public edm::one::EDAnalyzer<>
{
  public:
    explicit CTPPSFastSimulationValidator( const edm::ParameterSet& );
    ~CTPPSFastSimulationValidator();

  private:
    virtual void beginJob() override;

    virtual void analyze( const edm::Event&, const edm::EventSetup& ) override;

    virtual void endJob() override;

    edm::EDGetTokenT< std::vector<CTPPSLocalTrackLite> > simuTracksToken_;
    edm::EDGetTokenT< std::vector<CTPPSLocalTrackLite> > recoTracksToken_;

    std::string outputFile;

    struct RPPlots
    {
      TH2D *h2_xr_vs_xs=NULL, *h2_yr_vs_ys=NULL;
      TH1D *h_de_x, *h_de_y;

      void init()
      {
        h2_xr_vs_xs = new TH2D("", "", 100, -10., +10., 100, -10, +10.);
        h2_yr_vs_ys = new TH2D("", "", 100, -10., +10., 100, -10, +10.);
        h_de_x = new TH1D("", "", 100, -0., +0.);
        h_de_y = new TH1D("", "", 100, -0., +0.);
      }

      void fill(double simu_x, double simu_y, double reco_x, double reco_y)
      {
        if (h2_xr_vs_xs == NULL)
          init();

        h2_xr_vs_xs->Fill(simu_x, reco_x);
        h2_yr_vs_ys->Fill(simu_y, reco_y);

        h_de_x->Fill(reco_x - simu_x);
        h_de_y->Fill(reco_y - simu_y);
      }

      void write() const
      {
        h2_xr_vs_xs->Write("h2_xr_vs_xs");
        h2_yr_vs_ys->Write("h2_yr_vs_ys");
        h_de_x->Write("h_de_x");
        h_de_y->Write("h_de_y");
      }
    };

    std::map<unsigned int, RPPlots> rpPlots;
};

//----------------------------------------------------------------------------------------------------

CTPPSFastSimulationValidator::CTPPSFastSimulationValidator( const edm::ParameterSet& iConfig ) :
  simuTracksToken_( consumes< std::vector<CTPPSLocalTrackLite> >( iConfig.getParameter<edm::InputTag>( "simuTracksTag" ) ) ),
  recoTracksToken_( consumes< std::vector<CTPPSLocalTrackLite> >( iConfig.getParameter<edm::InputTag>( "recoTracksTag" ) ) ),
  outputFile( iConfig.getParameter<std::string>("outputFile") )
{
}

//----------------------------------------------------------------------------------------------------

CTPPSFastSimulationValidator::~CTPPSFastSimulationValidator()
{
}

//----------------------------------------------------------------------------------------------------

void CTPPSFastSimulationValidator::analyze( const edm::Event& iEvent, const edm::EventSetup& )
{
  // get input
  edm::Handle< std::vector<CTPPSLocalTrackLite> > simuTracks;
  iEvent.getByToken( simuTracksToken_, simuTracks );

  edm::Handle< std::vector<CTPPSLocalTrackLite> > recoTracks;
  iEvent.getByToken( recoTracksToken_, recoTracks );

  // process tracks
  for (const auto& simuTrack : *simuTracks)
  {

    for (const auto& recoTrack : *recoTracks)
    {
      if (simuTrack.getRPId() == recoTrack.getRPId())
      {
        CTPPSDetId rpId(simuTrack.getRPId());
        unsigned int rpDecId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();

        rpPlots[rpDecId].fill(simuTrack.getX(), simuTrack.getY(), recoTrack.getX(), recoTrack.getY());
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------

void CTPPSFastSimulationValidator::beginJob()
{
}

//----------------------------------------------------------------------------------------------------

void CTPPSFastSimulationValidator::endJob()
{
  TFile *f_out = TFile::Open(outputFile.c_str(), "recreate");
  
  for (const auto it : rpPlots)
  {
    gDirectory = f_out->mkdir(Form("RP %u", it.first));
    it.second.write();
  }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE( CTPPSFastSimulationValidator );
