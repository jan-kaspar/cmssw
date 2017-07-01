/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#ifndef _AlignmentTask_h_
#define _AlignmentTask_h_

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Alignment/CTPPS/interface/AlignmentGeometry.h"
class AlignmentConstraint;
class TotemRPGeometry;

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
        const std::vector<unsigned int> &excludedSensors, const TotemRPGeometry *,
        double z0, AlignmentGeometry &geometry);


    // -------------------- quantity-class-related members --------------------

    /// quantity classes
    enum QuantityClass
    {
      // TODO: change to ShR1, ShR2
      qcShR,    ///< detector shifts in readout direction
      qcShZ,    ///< detector shifts in z
      qcRotZ,   ///< detector rotations around z
    };

    /// list of quantity classes to be optimized
    std::vector<QuantityClass> quantityClasses;

    /// returns a string tag for the given quantity class
    std::string QuantityClassTag(QuantityClass) const;

    /// returns the number of quantities of the given class
    unsigned int QuantitiesOfClass(QuantityClass) const;

    // TODO: add matrixIndex mappings and methods


    // -------------------- constraint-related members --------------------

    /// returns the number of constraints of the given class
    unsigned int ConstraintsForClass(QuantityClass) const;
    
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
