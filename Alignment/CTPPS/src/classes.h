#include "Alignment/CTPPS/interface/AlignmentTask.h"
#include "Alignment/CTPPS/interface/AlignmentConstraint.h"
#include "Alignment/CTPPS/interface/AlignmentGeometry.h"
#include "Alignment/CTPPS/interface/IdealResult.h"
#include "Alignment/CTPPS/interface/JanAlignmentAlgorithm.h"

#include <TVectorD.h>
#include <TGraph.h>
#include <TH1D.h>

namespace {
  namespace {
    AlignmentTask at;

    AlignmentConstraint ac;
	std::map<unsigned int, TVectorD> muitvd;

	std::vector<AlignmentConstraint> vac;

	DetGeometry dg;
    AlignmentGeometry ag;

	std::map<unsigned int, DetGeometry> muidg;
	std::set<unsigned int> sui;

    IdealResult ir;
	
	SingularMode sm;
	std::vector<SingularMode> vsm;

    JanAlignmentAlgorithm jaa;
	JanAlignmentAlgorithm::ScatterPlot jaasp;
	JanAlignmentAlgorithm::DetStat jaads;
	std::map<unsigned int, JanAlignmentAlgorithm::DetStat> muids;
	std::vector<TH1D*> vth;
    std::vector<TGraph*> vtg;  
    std::map< std::set<unsigned int>, JanAlignmentAlgorithm::ScatterPlot> msuisp;
  }
}
