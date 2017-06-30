/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#ifndef _LocalTrackFitter_h_
#define _LocalTrackFitter_h_

#include "DataFormats/CTPPSAlignment/interface/LocalTrackFit.h"

#include "Alignment/CTPPS/interface/AlignmentGeometry.h"
#include "Alignment/CTPPS/interface/HitCollection.h"

namespace edm {
  class ParameterSet;
}

/**
 *\brief Performs straight-line fit and outlier rejection.
 **/
class LocalTrackFitter
{
  public:
    /// dummy constructor (not to be used)
    LocalTrackFitter() {}
    
    /// normal constructor
    LocalTrackFitter(const edm::ParameterSet&);
    
    virtual ~LocalTrackFitter() {}
    
    /// runs the fit and outlier-removal loop
    /// returns true in case of success
    // TODO: add const
    bool Fit(HitCollection&, const AlignmentGeometry&, LocalTrackFit&) const;

  protected:
    /// verbosity level
    unsigned int verbosity;
    
    /// minimum of hits to accept data from a RP
    unsigned int minimumHitsPerProjectionPerRP;
    
    /// hits with higher ratio residual/sigma will be dropped
    double maxResidualToSigma;
    
    /// fits the collection of hits and removes hits with too high residual/sigma ratio
    /// \param failed whether the fit has failed
    /// \param selectionChanged whether some hits have been removed
    // TODO: add const
    void FitAndRemoveOutliers(HitCollection&, const AlignmentGeometry&, LocalTrackFit&, 
      bool &failed, bool &selectionChanged) const;
    
    /// removes the hits of pots with too few planes active
    // TODO: add const
    void RemoveInsufficientPots(HitCollection&, bool &selectionChanged) const;
};

#endif
