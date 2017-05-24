#ifndef SimRomanPot_CTPPSOpticsParameterisation_ProtonReconstructionAlgorithm_h
#define SimRomanPot_CTPPSOpticsParameterisation_ProtonReconstructionAlgorithm_h

#include "TFile.h"
#include "TSpline.h"
#include "Fit/Fitter.h"

#include "LHCOpticsApproximator.h"

#include <map>
#include <string>
#include <cmath>

class ProtonReconstructionAlgorithm
{
  public:
    ProtonReconstructionAlgorithm() {}
    ~ProtonReconstructionAlgorithm();

    int Init(const std::string &optics_file_beam1, const std::string &optics_file_beam2);
    ProtonData Reconstruct(LHCSector sector, const TrackDataCollection &tracks) const;

  private:
    /// optics data associated with 1 RP
    struct RPOpticsData {
      LHCOpticsApproximator *optics;
      TSpline3 *s_xi_vs_x, *s_y0_vs_xi, *s_v_y_vs_xi, *s_L_y_vs_xi;
    };

    /// map: RP id --> optics data
    std::map<unsigned int, RPOpticsData> m_rp_optics;

    /// class for calculation of chi^2
    class ChiSquareCalculator {
      public:
        ChiSquareCalculator() {}

        const TrackDataCollection *tracks;
        const std::map<unsigned int, RPOpticsData> *m_rp_optics;
        double operator() (const double *parameters) const;
    };

    ChiSquareCalculator *chiSquareCalculator = NULL;

    // fitter object
    ROOT::Fit::Fitter *fitter = NULL;
};

#endif
