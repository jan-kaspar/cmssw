/****************************************************************************
*
* This is a part of CTPPS offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef _proton_reconstruction_h_
#define _proton_reconstruction_h_

#include "TFile.h"
#include "TSpline.h"
#include "Fit/Fitter.h"

#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

#include "RecoCTPPS/OpticsParametrization/interface/LHCOpticsApproximator.h"
#include "RecoCTPPS/ProtonReconstruction/interface/BeamConditions.h"
#include "RecoCTPPS/ProtonReconstruction/interface/ProtonData.h"

#include <map>
#include <string>
#include <cmath>

//----------------------------------------------------------------------------------------------------

class ProtonReconstruction
{
	public:

		ProtonReconstruction() {}

		int Init(const std::string &optics_file_beam1, const std::string &optics_file_beam2);

		ProtonData Reconstruct(LHCSector sector, const std::vector<CTPPSLocalTrackLite> &tracks) const;

		~ProtonReconstruction()
		{
			for (auto &it : m_rp_optics)
			{
				delete it.second.optics;
				delete it.second.s_xi_vs_x;
				delete it.second.s_y0_vs_xi;
				delete it.second.s_v_y_vs_xi;
				delete it.second.s_L_y_vs_xi;
			}

			if (chiSquareCalculator)
				delete chiSquareCalculator;

			if (fitter)
				delete fitter;
		}

	private:
		/// optics data associated with 1 RP
		struct RPOpticsData
		{
			LHCOpticsApproximator *optics;
			TSpline3 *s_xi_vs_x, *s_y0_vs_xi, *s_v_y_vs_xi, *s_L_y_vs_xi;
		};

		/// map: RP id --> optics data
		std::map<unsigned int, RPOpticsData> m_rp_optics;

		/// class for calculation of chi^2
		class ChiSquareCalculator
		{
			public:
				const std::vector<CTPPSLocalTrackLite> *tracks;
				const std::map<unsigned int, RPOpticsData> *m_rp_optics;

				ChiSquareCalculator() {}

				double operator() (const double *parameters) const;
		};

		/// instance of ChiSquareCalculator used by the fitter
		ChiSquareCalculator *chiSquareCalculator = NULL;

		/// fitter object
		ROOT::Fit::Fitter *fitter = NULL;
};

#endif
