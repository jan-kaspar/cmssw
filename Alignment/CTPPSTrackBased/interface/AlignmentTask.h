/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#ifndef Alignment_CTPPSTrackBased_AlignmentTask_h
#define Alignment_CTPPSTrackBased_AlignmentTask_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Alignment/CTPPSTrackBased/interface/AlignmentGeometry.h"
class AlignmentConstraint;
class CTPPSGeometry;

#include <vector>


/**
 *\brief Represents an alignment task.
 **/
class AlignmentTask
{
  public:
    // -------------------- config file parameters --------------------

    /// whether to resolve detector shifts in readout direction(s)
    bool resolveShR;

    /// whether to resolve detector shifts in z
    bool resolveShZ;

    /// whether to resolve detector rotations around z
    bool resolveRotZ;


    /// whether the per-group constraint shall be applied
    bool useExtendedRotZConstraint;
    
    /// whether the per-group constraint shall be applied
    bool useZeroThetaRotZConstraint;

    /// whether the per-group constraints shall be applied
    bool useExtendedShZConstraints;

    /// whether to apply the constraint mean U = mean V RotZ ("standard" set of constraints only)
    bool useEqualMeanUMeanVRotZConstraint;

    /// whether to resolve only 1 rot_z per RP
    bool oneRotZPerPot;


    /// homogeneous constraints from config file
    edm::ParameterSet homogeneousConstraints;

    /// fixed detectors constraints from config file
    edm::ParameterSet fixedDetectorsConstraints;

    /// settings of "standard" constraints from config file
    edm::ParameterSet standardConstraints;


    // -------------------- geometry-related members --------------------

    /// the geometry for this task
    AlignmentGeometry geometry;

    /// builds the alignment geometry
    static void BuildGeometry(const std::vector<unsigned int> &rpDecIds,
        const std::vector<unsigned int> &excludedSensors, const CTPPSGeometry *,
        double z0, AlignmentGeometry &geometry);


    // -------------------- quantity-class-related members --------------------

    /// quantity classes
    enum QuantityClass
    {
      qcShR1,   ///< detector shifts in first readout direction
      qcShR2,   ///< detector shifts in second readout direction
      qcShZ,    ///< detector shifts in z
      qcRotZ,   ///< detector rotations around z
    };

    /// list of quantity classes to be optimized
    std::vector<QuantityClass> quantityClasses;

    /// returns a string tag for the given quantity class
    std::string QuantityClassTag(QuantityClass) const;

    struct DetIdDirIdxPair
    {
      unsigned int detId;
      unsigned int dirIdx;
    
      bool operator< (const DetIdDirIdxPair &other) const
      {
        if (detId < other.detId)
          return true;
        if (detId > other.detId)
          return false;
        if (dirIdx < other.dirIdx)
          return true;

        return false;
      }
    };

    /// for each quantity class contains mapping (detector id, direction) --> measurement index
    std::map<QuantityClass, std::map<DetIdDirIdxPair, unsigned int>> mapMeasurementIndeces;

    /// for each quantity class contains mapping detector id --> quantity index
    std::map<QuantityClass, std::map<unsigned int, unsigned int>> mapQuantityIndeces;

    /// builds "mapMatrixIndeces" from "geometry"
    void BuildIndexMaps();

    /// returns the number of quantities of the given class
    unsigned int MeasurementsOfClass(QuantityClass) const;

    /// returns the number of quantities of the given class
    unsigned int QuantitiesOfClass(QuantityClass) const;

    /// returns measurement index (if non-existent, returns -1)
    signed int GetMeasurementIndex(QuantityClass cl, unsigned int detId, unsigned int dirIdx) const;

    /// returns measurement index (if non-existent, returns -1)
    signed int GetQuantityIndex(QuantityClass cl, unsigned int detId) const;


    // -------------------- constraint-related members --------------------

    /// builds a set of homogeneous constraints
    void BuildHomogeneousConstraints(std::vector<AlignmentConstraint>&) const;
    
    /// builds a set of fixed-detector constraints
    void BuildFixedDetectorsConstraints(std::vector<AlignmentConstraint>&) const;
    
    /// builds the standard constraints
    void BuildStandardConstraints(std::vector<AlignmentConstraint>&) const;

    /// adds constraints such that only 1 rot_z per RP is left
    void BuildOneRotZPerPotConstraints(std::vector<AlignmentConstraint>&) const;


    // -------------------- constructors --------------------
 
    /// dummy constructor (not to be used)
    AlignmentTask()
    {
    }
    
    /// normal constructor
    AlignmentTask(const edm::ParameterSet& ps);
};

#endif
