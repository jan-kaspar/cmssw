/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "Alignment/CTPPS/interface/AlignmentTask.h"
#include "Alignment/CTPPS/interface/AlignmentConstraint.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"

#include <algorithm>

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;
using namespace CLHEP;

//----------------------------------------------------------------------------------------------------


AlignmentTask::AlignmentTask(const ParameterSet& ps) :
  resolveShR(ps.getParameter<bool>("resolveShR")),
  resolveShZ(ps.getParameter<bool>("resolveShZ")),
  resolveRotZ(ps.getParameter<bool>("resolveRotZ")),

  useExtendedRotZConstraint(ps.getParameter<bool>("useExtendedRotZConstraint")),
  useZeroThetaRotZConstraint(ps.getParameter<bool>("useZeroThetaRotZConstraint")),
  useExtendedShZConstraints(ps.getParameter<bool>("useExtendedShZConstraints")),
  useEqualMeanUMeanVRotZConstraint(ps.getParameter<bool>("useEqualMeanUMeanVRotZConstraint")),
  oneRotZPerPot(ps.getParameter<bool>("oneRotZPerPot")),

  homogeneousConstraints(ps.getParameterSet("homogeneousConstraints")),
  fixedDetectorsConstraints(ps.getParameterSet("fixedDetectorsConstraints")),
  standardConstraints(ps.getParameterSet("standardConstraints"))
{
  if (resolveShR)
    quantityClasses.push_back(qcShR);

  if (resolveShZ)
    quantityClasses.push_back(qcShZ);

  if (resolveRotZ)
    quantityClasses.push_back(qcRotZ);
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildGeometry(const vector<unsigned int> &rpDecIds,
  const vector<unsigned int> &excludedSensors, const TotemRPGeometry *input, double z0, AlignmentGeometry &geometry)
{
  geometry.clear();
  geometry.z0 = z0; 

  // traverse full known geometry
  for (auto it = input->beginDet(); it != input->endDet(); ++it)
  {
    // skip excluded sensors
    if (find(excludedSensors.begin(), excludedSensors.end(), it->first) != excludedSensors.end())
      continue;

    // is RP selected?
    const CTPPSDetId detId(it->first);
    const unsigned int rpDecId = 100*detId.arm() + 10*detId.station() + detId.rp();
    if (find(rpDecIds.begin(), rpDecIds.end(), rpDecId) == rpDecIds.end())
      continue;

    // extract geometry data
    Hep3Vector c = input->LocalToGlobal(detId, Hep3Vector(0., 0., 0.));
    Hep3Vector d1 = input->LocalToGlobal(detId, Hep3Vector(1., 0., 0.)) - c;
    Hep3Vector d2 = input->LocalToGlobal(detId, Hep3Vector(0., 1., 0.)) - c;

    // for strips: is it U plane?
    bool isU = false;
    if (detId.subdetId() == CTPPSDetId::sdTrackingStrip)
    {
      TotemRPDetId stripDetId(it->first);
      unsigned int rpNum = stripDetId.rp();
      unsigned int plNum = stripDetId.plane();
      isU = (plNum % 2 != 0);
      if (rpNum == 2 || rpNum == 3)
        isU = !isU;
    }

    DetGeometry dg(c.z() - z0, c.x(), c.y(), isU);
    dg.SetDirection(1, d1.x(), d1.y(), d1.z());
    dg.SetDirection(2, d2.x(), d2.y(), d2.z());
    geometry.Insert(it->first, dg);
  }

  // set matrix and rpMatrix indeces
  unsigned int index = 0;
  for (AlignmentGeometry::iterator it = geometry.begin(); it != geometry.end(); ++it, ++index)
  {
    it->second.matrixIndex = index;
  }
}

//----------------------------------------------------------------------------------------------------

string AlignmentTask::QuantityClassTag(QuantityClass qc) const
{
  switch (qc)
  {
    case qcShR: return "ShR"; 
    case qcShZ: return "ShZ"; 
    case qcRotZ: return "RotZ"; 
  }

  throw cms::Exception("AlignmentTask::QuantityClassTag") << "Unknown quantity class " << qc << ".";
}

//----------------------------------------------------------------------------------------------------

unsigned int AlignmentTask::QuantitiesOfClass(QuantityClass qc) const
{
  // TODO
  return geometry.Detectors();
}

//----------------------------------------------------------------------------------------------------

unsigned int AlignmentTask::ConstraintsForClass(QuantityClass qc) const
{
  switch (qc)
  {
    case qcShR: return 4; 
    case qcShZ: return (useExtendedShZConstraints) ? 4 : 2;
    case qcRotZ:
      if (oneRotZPerPot)
      {
        return 9*geometry.RPs() + 1;
      } else {
        unsigned int count = (useZeroThetaRotZConstraint) ? 2 : 1;
        if (useExtendedRotZConstraint)
          count *= 2;
        return count;
      }
  }

  throw cms::Exception("AlignmentTask::ConstraintsForClass") << "Unknown quantity class " << qc << ".";
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildHomogeneousConstraints(vector<AlignmentConstraint> &constraints) const
{
  printf(">> AlignmentTask::BuildHomogeneousConstraints > not yet ported, sorry!\n");
  throw 1;

#if 0
  for (unsigned int cl = 0; cl < quantityClasses.size(); cl++)
  {
    unsigned int size = ConstraintsForClass(quantityClasses[cl]);
    const string &tag = QuantityClassTag(quantityClasses[cl]);
    
    // just one basic constraint
    if (oneRotZPerPot && quantityClasses[cl] == qcRotZ)
      size = 1;

    // get constraint values
    char buf[20];
    sprintf(buf, "%s_values", tag.c_str());
    vector<double> values(homogeneousConstraints.getParameter< vector<double> >(buf));
    if (values.size() < size)
      throw cms::Exception("AlignmentTask::BuildHomogeneousConstraints") <<
        "Invalid number of constraint values for " << tag << ". Given " << values.size() <<
        ", expected " << size << ".";

    for (unsigned int j = 0; j < size; j++)
    {
      // prepare a constraint with coefficient vectors
      AlignmentConstraint ac;
      ac.forClass = quantityClasses[cl];
      for (unsigned int i = 0; i < quantityClasses.size(); i++)
      {
        ac.coef[quantityClasses[i]].ResizeTo(QuantitiesOfClass(quantityClasses[i]));
        ac.coef[quantityClasses[i]].Zero();
      }
      ac.val = values[j];
      ac.extended = false;

      unsigned int indeces = QuantitiesOfClass(quantityClasses[cl]);

      for (unsigned int idx = 0; idx < indeces; ++idx)
      {
        double &coef = ac.coef[quantityClasses[cl]][idx];
        const DetGeometry &dt = geometry.FindByMatrixIndex(idx)->second;
        double sc = -dt.d2x*dt.sy + dt.d2y*dt.sx;

        switch (quantityClasses[cl])
        {
          case qcShR:
            switch (j)
            {
              case 0: ac.name = "ShR: z*dx"; coef = dt.z * dt.d2x; break;
              case 1: ac.name = "ShR: dx"; coef = dt.d2x; break;
              case 2: ac.name = "ShR: z*dy"; coef = dt.z * dt.d2y; break;
              case 3: ac.name = "ShR: dy"; coef = dt.d2y; break;
            }
            break;

          case qcShZ:
            switch (j)
            {
              case 0: ac.name = "ShZ: z"; coef = dt.z; break;
              case 1: ac.name = "ShZ: 1"; coef = 1.; break;
              case 2: ac.name = "ShZ: z for V-det"; coef = (dt.isU) ? 0. : dt.z; ac.extended = true; break;
              case 3: ac.name = "ShZ: 1 for V-det"; coef = (dt.isU) ? 0. : 1.; ac.extended = true; break;
            }
            break;
        
          case qcRPShZ:
            switch (j)
            {
              case 0: ac.name = "RPShZ: 1"; coef = 1.; break;
              case 1: ac.name = "RPShZ: z"; coef = dt.z; ac.extended = true; break;
            }
            break;
          
          case qcRotZ:
            unsigned int je = j;
            if (!useExtendedRotZConstraint)
              je *= 2;
            switch (je)
            {
              case 0: ac.name = "RotZ: 1 all det"; coef = 1.; break;
              case 1: ac.name = "RotZ: 1 V-det"; coef = (dt.isU) ? 0. : 1.; ac.extended = true; break;
              case 2: ac.name = "RotZ: z all det"; coef = dt.z; ac.extended = true; break;
              case 3: ac.name = "RotZ: z V-det"; coef = (dt.isU) ? 0. : dt.z; ac.extended = true; break;
            }
            ac.coef[qcShR][idx] = sc * ac.coef[qcRotZ][idx];
            break;
        }
      }

      constraints.push_back(ac);
    } 
    
    // only 1 rot_z per RP
    if (oneRotZPerPot && quantityClasses[cl] == qcRotZ)
      BuildOneRotZPerPotConstraints(constraints);
  }
#endif
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildFixedDetectorsConstraints(vector<AlignmentConstraint> &constraints) const
{
  // TODO
  //printf(">> AlignmentTask::BuildFixedDetectorsConstraints > not yet ported.\n");
  //throw 1;

  for (unsigned int cl = 0; cl < quantityClasses.size(); cl++)
  {
    const string &tag = QuantityClassTag(quantityClasses[cl]);

    unsigned int basicNumber = 0;
    switch (quantityClasses[cl])
    {
      case qcShR: basicNumber = 4; break;
      case qcRotZ: basicNumber = 1; break;
      case qcShZ: basicNumber = 2; break;
    }

    const ParameterSet &classSettings = fixedDetectorsConstraints.getParameterSet(tag.c_str());
    vector<unsigned int> ids(classSettings.getParameter< vector<unsigned int> >("ids"));
    vector<double> values(classSettings.getParameter< vector<double> >("values"));

    if (ids.size() != values.size())
      throw cms::Exception("AlignmentTask::BuildFixedDetectorsConstraints") << 
        "Different number of constraint ids and values for " << tag << ".";
    
    unsigned int size = ConstraintsForClass(quantityClasses[cl]);
    
    // just one basic constraint
    if (oneRotZPerPot && quantityClasses[cl] == qcRotZ)
      size = 1;
    
    if (ids.size() < size)
      throw cms::Exception("AlignmentTask::BuildFixedDetectorsConstraints") << 
        "Too few constrainted ids for " << tag << ". Given " << ids.size() <<
        ", while " << size << " expected.";
    
    for (unsigned int j = 0; j < size; j++)
    {
      AlignmentConstraint ac;
      ac.forClass = quantityClasses[cl];
      for (unsigned int i = 0; i < quantityClasses.size(); i++)
      {
        ac.coef[quantityClasses[i]].ResizeTo(QuantitiesOfClass(quantityClasses[i]));
        ac.coef[quantityClasses[i]].Zero();
      }
      ac.val = values[j] * 1E-3;
      ac.extended = (j >= basicNumber);

      char buf[40];
      sprintf(buf, "%s: fixed plane %4u", tag.c_str(), ids[j]);
      ac.name = buf;

      // is the detector in geometry?
      if (!geometry.ValidSensorId(ids[j]))
        throw cms::Exception("AlignmentTask::BuildFixedDetectorsConstraints") <<
          "Detector with id " << ids[j] << " is not in the geometry.";

      const auto git = geometry.find(ids[j]);
      unsigned int idx = git->second.matrixIndex;
      ac.coef[quantityClasses[cl]][idx] = 1.;

      constraints.push_back(ac);
    }
    
    if (oneRotZPerPot && quantityClasses[cl] == qcRotZ)
        BuildOneRotZPerPotConstraints(constraints);
  }
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildOneRotZPerPotConstraints(std::vector<AlignmentConstraint> &constraints) const
{
  printf(">> AlignmentTask::BuildOneRotZPerPotConstraints > not yet ported\n");
  throw 1;

#if 0
  // geometry is sorted by the detector number
  unsigned int prev_rp = 12345;
  for (AlignmentGeometry::iterator it = geometry.begin(); it != geometry.end(); ++it)
  {
    // do not mix different pots
    if (it->first / 10 != prev_rp)
    {
      prev_rp = it->first / 10;
      continue;
    }

    AlignmentConstraint ac;
    ac.forClass = qcRotZ;
    for (unsigned int i = 0; i < quantityClasses.size(); i++)
    {
      ac.coef[quantityClasses[i]].ResizeTo(QuantitiesOfClass(quantityClasses[i]));
      ac.coef[quantityClasses[i]].Zero();
    }
    ac.val = 0;
    ac.extended = true;
    ac.coef[qcRotZ][geometry[it->first].matrixIndex] = 1.;
    ac.coef[qcRotZ][geometry[it->first-1].matrixIndex] = -1.;

    char buf[40];
    sprintf(buf, "RotZ: %u = %u", it->first, it->first-1);
    ac.name = buf;

    constraints.push_back(ac);
  }
#endif
}

//----------------------------------------------------------------------------------------------------

void AlignmentTask::BuildStandardConstraints(vector<AlignmentConstraint> &constraints) const
{
  // TODO

#if 0
  const vector<unsigned int> &units = standardConstraints.getParameter<vector<unsigned int>>("units");

  // count planes in RPs
  struct PlaneCount { unsigned int all, U, V; };
  map<unsigned int, PlaneCount> planesPerPot;
  set<unsigned int> rpsToConstrainShR, rpsToConstrainRotZ;
  for (const auto it : geometry)
  {
    unsigned int rp = it.first / 10;
    planesPerPot[rp].all++;
    if (it.second.isU)
      planesPerPot[rp].U++;
    else
      planesPerPot[rp].V++;
  }

  // ShR constraints
  if (resolveShR)
  {
    for (const auto &coUnit : units)
    {
      // prepare empty constraints
      AlignmentConstraint ac_X;
      ac_X.forClass = qcRotZ;
      for (unsigned int i = 0; i < quantityClasses.size(); i++)
      {
        ac_X.coef[quantityClasses[i]].ResizeTo(QuantitiesOfClass(quantityClasses[i]));
        ac_X.coef[quantityClasses[i]].Zero();
      }
      ac_X.val = 0;
      ac_X.extended = false;
    
      AlignmentConstraint ac_Y(ac_X);
  
      char buf[50];
      sprintf(buf, "ShR: unit %u, SmX=0", coUnit);
      ac_X.name = buf;
      sprintf(buf, "ShR: unit %u, SmY=0", coUnit);
      ac_Y.name = buf;
  
      unsigned int indeces = QuantitiesOfClass(qcShR);
      for (unsigned int idx = 0; idx < indeces; ++idx)
      {
        const auto &entry = geometry.FindByMatrixIndex(idx);
        unsigned int plId = entry->first;
        const DetGeometry &dt = entry->second;
  
        unsigned int rpId = plId / 10;
        unsigned int rpNum = rpId % 10;
        unsigned int stNum = (rpId / 10) % 10;
        unsigned int unNum = (rpNum > 2) ? 1 : 0;     
        unsigned int unitIdx = 10*stNum + unNum;
        bool hor = (rpNum == 2 || rpNum == 3);
  
        if (unitIdx == coUnit && !hor)
        {
          const auto &pl = planesPerPot[rpId];
          ac_X.coef[qcShR][idx] = dt.d2x / pl.all;
          ac_Y.coef[qcShR][idx] = dt.d2y / pl.all;
        }
      }
  
      constraints.push_back(ac_X);
      constraints.push_back(ac_Y);
    }
  }

  // RotZ constraints
  if (resolveRotZ)
  {
    for (const auto &coUnit : units)
    {
      // prepare empty constraints
      AlignmentConstraint ac_U;
      ac_U.forClass = qcRotZ;
      for (unsigned int i = 0; i < quantityClasses.size(); i++)
      {
        ac_U.coef[quantityClasses[i]].ResizeTo(QuantitiesOfClass(quantityClasses[i]));
        ac_U.coef[quantityClasses[i]].Zero();
      }
      ac_U.val = 0;
      ac_U.extended = true;
    
      AlignmentConstraint ac_V(ac_U);
  
      char buf[50];
      sprintf(buf, "RotZ: unit %u, SmU=0", coUnit);
      ac_U.name = buf;
      sprintf(buf, "RotZ: unit %u, SmV=0", coUnit);
      ac_V.name = buf;
  
      unsigned int indeces = QuantitiesOfClass(qcRotZ);
      for (unsigned int idx = 0; idx < indeces; ++idx)
      {
        const auto &entry = geometry.FindByMatrixIndex(idx);
        unsigned int plId = entry->first;
        const DetGeometry &dt = entry->second;
  
        unsigned int rpId = plId / 10;
        unsigned int rpNum = rpId % 10;
        unsigned int stNum = (rpId / 10) % 10;
        unsigned int unNum = (rpNum > 2) ? 1 : 0;     
        unsigned int unitIdx = 10*stNum + unNum;
  
        if (unitIdx == coUnit)
        {
          const auto &pl = planesPerPot[rpId];
          ac_U.coef[qcRotZ][idx] = (dt.isU) ? 1. / pl.U : 0.;
          ac_V.coef[qcRotZ][idx] = (!dt.isU) ? 1. / pl.V : 0.;
        }
      }
  
      constraints.push_back(ac_U);
  
      if (!useEqualMeanUMeanVRotZConstraint)
        constraints.push_back(ac_V);
    }
  }

  // mean U = mean V RotZ constraints
  if (resolveRotZ && useEqualMeanUMeanVRotZConstraint)
  {
    map<unsigned int, pair<unsigned int, unsigned int>> planesPerPot;
    for (const auto it : geometry)
    {
      unsigned int rp = it.first / 10;
      if (it.second.isU)
        planesPerPot[rp].first++;
      else
        planesPerPot[rp].second++;
    }

    for (const auto it : planesPerPot)
    {
      const unsigned int &rp = it.first;
      const unsigned int &rp_u_planes = it.second.first;
      const unsigned int &rp_v_planes = it.second.second;

      // prepare empty constraint
      AlignmentConstraint ac;
      ac.forClass = qcRotZ;
      for (unsigned int i = 0; i < quantityClasses.size(); i++)
      {
        ac.coef[quantityClasses[i]].ResizeTo(QuantitiesOfClass(quantityClasses[i]));
        ac.coef[quantityClasses[i]].Zero();
      }
      ac.val = 0;
      ac.extended = true;
  
      char buf[40];
      sprintf(buf, "RotZ: RP%u, mU=mV", rp);
      ac.name = buf;

      // set the coefficients
      const unsigned int indeces = QuantitiesOfClass(qcRotZ);
      for (unsigned int idx = 0; idx < indeces; idx++)
      {
        double &coef = ac.coef[qcRotZ][idx];
        const auto &entry = geometry.FindByMatrixIndex(idx);
        unsigned int idx_det = entry->first;
        unsigned int idx_rp = idx_det / 10;
        bool idx_isU = entry->second.isU;

        if (rp == idx_rp)
          coef = (idx_isU) ? +1./rp_u_planes : -1./rp_v_planes;
      }
  
      // push the constraint
      constraints.push_back(ac);
    }
  }

  // print constraints
  for (const auto &c : constraints)
  {
    printf("* name=%s, value=%.3f, forClass=%i\n", c.name.c_str(), c.val, c.forClass);

    for (const auto &ccl : c.coef)
    {
      printf("    class %u (%s)\n", ccl.first, QuantityClassTag(QuantityClass(ccl.first)).c_str());

      for (int i = 0; i < ccl.second.GetNrows(); i++)
        printf("    %u -> %+.4f\n", i, ccl.second(i));
    }  
  }
#endif
 
  // TODO: when oneRotZPerPot=true, make use of BuildOneRotZPerPotConstraints(constraints) ??
}
