/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Alignment/CTPPS/interface/LocalTrackFitter.h"

#include <set>

#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

LocalTrackFitter::LocalTrackFitter(const edm::ParameterSet &ps) :
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
  minimumHitsPerProjectionPerRP(ps.getParameter<unsigned int>("minimumHitsPerProjectionPerRP")),
  maxResidualToSigma(ps.getParameter<double>("maxResidualToSigma"))
{
}

//----------------------------------------------------------------------------------------------------

bool LocalTrackFitter::Fit(HitCollection &selection, const AlignmentGeometry &geometry, LocalTrackFit &trackFit) const
{
  if (verbosity > 5)
    printf(">> LocalTrackFitter::Fit\n");

  bool selectionChanged = true;
  unsigned int loopCounter = 0;
  while (selectionChanged)
  {
    // fit/outlier-removal loop
    while (selectionChanged)
    {
      if (verbosity > 5)
        printf("* fit loop %u\n", loopCounter++);

      bool fitFailed = false;
      FitAndRemoveOutliers(selection, geometry, trackFit, fitFailed, selectionChanged);

      if (fitFailed)
      {
        if (verbosity > 5)
          printf("\tFIT FAILED\n");
        return false;
      }
    }

    // remove pots with too few active planes
    if (verbosity > 5)
      printf("* removing insufficient pots\n");
    RemoveInsufficientPots(selection, selectionChanged);
  }

  return true;
}

//----------------------------------------------------------------------------------------------------

void LocalTrackFitter::FitAndRemoveOutliers(HitCollection &selection, const AlignmentGeometry &geometry, 
  LocalTrackFit &trackFit, bool &failed, bool &selectionChanged) const
{
  if (verbosity > 5)
    printf(">> LocalTrackFitter::FitAndRemoveOutliers\n");

  if (selection.empty())
  {
    failed = true;
    return;
  }
  
  // build matrices and vectors
  TMatrixD A(selection.size(), 4);
  TMatrixD Vi(selection.size(), selection.size());
  TVectorD measVec(selection.size());
  unsigned int j = 0;
  for (HitCollection::iterator it = selection.begin(); it != selection.end(); ++it, ++j)
  {
	const unsigned int &detId = it->id;

    // make sure detector is in geometry
    AlignmentGeometry::const_iterator dit = geometry.find(detId);
    if (dit == geometry.end())
    {
      printf("ERROR in LocalTrackFitter::FitAndRemoveOutliers > detector %u not in geometry.\n", detId);
      continue;
    }

    const DetGeometry &d = dit->second;
    const auto &dirData = d.GetDirectionData(it->dirIdx);

    // TODO: use z from hit
    A(j, 0) = d.z * dirData.dx;
    A(j, 1) = dirData.dx;
    A(j, 2) = d.z * dirData.dy;
    A(j, 3) = dirData.dy;

    measVec(j) = it->position + dirData.s;  // in mm

    Vi(j, j) = 1. / it->sigma / it->sigma;
  }
  //Print(A, "alpha");

  // evaluate local track parameter estimates (h stands for hat)
  TMatrixD AT(4, selection.size());
  AT.Transpose(A);
  TMatrixD ATViA(4, 4);
  ATViA = AT * Vi * A;
  TMatrixD ATViAI(ATViA);
  try
  {
    ATViAI = ATViA.Invert();
  }
  catch (...)
  {
    failed = true;
    return;
  }
  TVectorD theta(4);
  theta = ATViAI * AT * Vi * measVec;

  // residuals
  TVectorD R(measVec);
  R -= A * theta;

  // save results to trackFit
  trackFit.ax = theta(0);
  trackFit.bx = theta(1);
  trackFit.ay = theta(2);
  trackFit.by = theta(3);
  trackFit.z0 = geometry.z0;
  trackFit.ndf = selection.size() - 4;
  trackFit.chi_sq = 0;
  for (int i = 0; i < R.GetNrows(); i++)
    trackFit.chi_sq += R(i)*R(i)*Vi(i, i);
  
  if (verbosity > 5)
  {
    printf("    ax = %.3f mrad, bx = %.4f mm, ay = %.3f mrad, by = %.4f mm, z0 = %.3f mm\n",
      trackFit.ax*1E3, trackFit.bx, trackFit.ay*1E3, trackFit.by, trackFit.z0);
    printf("    ndof = %i, chi^2/ndof/si^2 = %.3f\n", trackFit.ndf, trackFit.chi_sq / trackFit.ndf);
  }

  // check residuals
  selectionChanged = false;
  TVectorD interpolation(A * theta);
  j = 0;
  for (HitCollection::iterator it = selection.begin(); it != selection.end(); ++j)
  {
    if (verbosity > 5)
      printf("        %2u, %4u: interpolation = %+8.1f um, R = %+6.1f um, R / sigma = %+6.2f\n", j, 
        it->id, interpolation[j]*1E3, R[j]*1E3, R[j]/it->sigma);

    double resToSigma = R[j] / it->sigma;
    if (fabs(resToSigma) > maxResidualToSigma)
    {
      selection.erase(it);
      selectionChanged = true;
      if (verbosity > 5)
        printf("            Removed\n");
    } else
      ++it;
  } 
}

//----------------------------------------------------------------------------------------------------

void LocalTrackFitter::RemoveInsufficientPots(HitCollection &selection, bool &selectionChanged) const
{
  selectionChanged = false;

// TODO
#if 0
  // map: rp id -> (active u planes, active v planes)
  map<unsigned int, pair< set<unsigned int>, set<unsigned int> > > planeMap;
  for (HitCollection::iterator it = selection.begin(); it != selection.end(); ++it)
  {
    unsigned int decId = it->id;
    unsigned int rpId = decId / 10;
    if (TotRPDetId::IsStripsCoordinateUDirection(decId))
      planeMap[rpId].first.insert(decId);
    else
      planeMap[rpId].second.insert(decId);
  }

  selectionChanged = false;
  for (map<unsigned int, pair< set<unsigned int>, set<unsigned int> > >::iterator it = planeMap.begin();
    it != planeMap.end(); ++it)
  {
    if (it->second.first.size() < minimumHitsPerProjectionPerRP
        || it->second.second.size() < minimumHitsPerProjectionPerRP)
    {
      if (verbosity > 5)
        printf("\tRP %u: u=%lu, v=%lu\n", it->first, it->second.first.size(), it->second.second.size());
      
      // remove all hits from that RP
      for (HitCollection::iterator dit = selection.begin(); dit != selection.end();)
      {
        if (it->first == dit->id/10)
        {
          if (verbosity > 5)
            printf("\t\tremoving %u\n", dit->id);
          selection.erase(dit);
          selectionChanged = true;
        } else
          dit++;
      }
    }
  }
#endif
}
