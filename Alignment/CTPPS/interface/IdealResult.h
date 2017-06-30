/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"

#include "DataFormats/CTPPSAlignment/interface/RPAlignmentCorrectionData.h"
#include "DataFormats/CTPPSAlignment/interface/RPAlignmentCorrectionsData.h"

#include "Alignment/CTPPS/interface/AlignmentAlgorithm.h"

#include <vector>

class AlignmentTask;


/**
 *\brief Calculates the ideal result of the StraightTrackAlignment.
 **/
class IdealResult : public AlignmentAlgorithm
{
  protected:
    edm::ESHandle<TotemRPGeometry> gReal, gMisaligned;

    bool useExtendedConstraints;

  public:
    /// dummy constructor (not to be used)
    IdealResult() {}

    /// normal constructor
    IdealResult(const edm::ParameterSet& ps, AlignmentTask *_t);

    virtual ~IdealResult() {}

    virtual std::string GetName() override
      { return "Ideal"; }

    virtual bool HasErrorEstimate() override
      { return false; }

    virtual void Begin(const edm::EventSetup&) override;

    virtual void Feed(const HitCollection&, const LocalTrackFit&) override {}

    virtual void SaveDiagnostics(TDirectory *) override {}

    virtual std::vector<SingularMode> Analyze() override;

    virtual unsigned int Solve(const std::vector<AlignmentConstraint>&, RPAlignmentCorrectionsData &result, TDirectory *dir) override;

    virtual void End() override {}
};
