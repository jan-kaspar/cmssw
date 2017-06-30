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
#include "FWCore/Framework/interface/ESWatcher.h"

#include "Alignment/CTPPS/interface/StraightTrackAlignment.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"

/**
 *\brief An EDAnalyzer that runs StraightTrackAlignment.
 **/
class CTPPSStraightTrackAligner : public edm::EDAnalyzer
{
  public:
    CTPPSStraightTrackAligner(const edm::ParameterSet &ps); 
    ~CTPPSStraightTrackAligner() {}

  private:
    unsigned int verbosity;

    bool worker_initialized;
    StraightTrackAlignment worker;

    edm::ESWatcher<VeryForwardRealGeometryRecord> geometryWatcher;

    virtual void beginJob() {}
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event &e, const edm::EventSetup &es);
    virtual void endJob();
};

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

CTPPSStraightTrackAligner::CTPPSStraightTrackAligner(const ParameterSet &ps) : 
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
  worker_initialized(false),
  worker(ps)
{
}

//----------------------------------------------------------------------------------------------------

void CTPPSStraightTrackAligner::beginRun(edm::Run const&, edm::EventSetup const& es)
{
}

//----------------------------------------------------------------------------------------------------

void CTPPSStraightTrackAligner::analyze(const edm::Event &e, const edm::EventSetup &es)
{
  if (geometryWatcher.check(es))
  {
    if (worker_initialized)
      throw cms::Exception("CTPPSStraightTrackAligner") <<
        "CTPPSStraightTrackAligner can't cope with changing geometry - change in event " << e.id() << endl;
  }

  if (!worker_initialized)
  {
    worker.Begin(es);
    worker_initialized = true;
  }

  worker.ProcessEvent(e, es);
}

//----------------------------------------------------------------------------------------------------

void CTPPSStraightTrackAligner::endJob()
{
  worker.Finish();
}

DEFINE_FWK_MODULE(CTPPSStraightTrackAligner);
