/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#ifndef Alignment_CTPPSTrackBased_AlignmentAlgorithm_h
#define Alignment_CTPPSTrackBased_AlignmentAlgorithm_h

#include "FWCore/Framework/interface/EventSetup.h"

#include "Alignment/CTPPSTrackBased/interface/LocalTrackFit.h"
#include "Alignment/CTPPSTrackBased/interface/AlignmentGeometry.h"
#include "Alignment/CTPPSTrackBased/interface/HitCollection.h"
#include "Alignment/CTPPSTrackBased/interface/SingularMode.h"
#include "Alignment/CTPPSTrackBased/interface/AlignmentConstraint.h"
#include "Alignment/CTPPSTrackBased/interface/AlignmentResult.h"

#include <string>
#include <map>

class AlignmentTask;
class TDirectory;

namespace edm {
  class ParameterSet;
}

/**
 *\brief Abstract parent for all (track-based) alignment algorithms
 **/
class AlignmentAlgorithm
{
  protected:
    unsigned int verbosity;

    /// the tasked to be completed
    AlignmentTask *task;

    /// eigenvalues in (-singularLimit, singularLimit) are treated as singular
    double singularLimit;

  public:
    /// dummy constructor (not to be used)
    AlignmentAlgorithm() {}
    
    /// normal constructor
    AlignmentAlgorithm(const edm::ParameterSet& ps, AlignmentTask *_t);

    virtual ~AlignmentAlgorithm() {}

    virtual std::string GetName()
      { return "Base"; }

    /// returns whether this algorithm is capable of estimating result uncertainties
    virtual bool HasErrorEstimate() = 0;

    /// prepare for processing
    virtual void Begin(const edm::EventSetup&) = 0;

    /// process one track
    virtual void Feed(const HitCollection&, const LocalTrackFit&) = 0;

    /// saves diagnostic histograms/plots
    virtual void SaveDiagnostics(TDirectory *) = 0;
    
    /// analyzes the data collected, returns a list of singular modes
    virtual std::vector<SingularMode> Analyze() = 0;

    /// solves the alignment problem with the given constraints
    /// \param dir a directory (in StraightTrackAlignment::taskDataFileName) where
    /// intermediate results can be stored
    virtual unsigned int Solve(const std::vector<AlignmentConstraint>&,
      std::map<unsigned int, AlignmentResult> &results, TDirectory *dir = NULL) = 0;

    /// cleans up after processing
    virtual void End() = 0;
};

#endif
