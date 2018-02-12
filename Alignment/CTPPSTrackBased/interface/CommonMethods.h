/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *	Jan Kašpar (jan.kaspar@gmail.com)
 *
 ****************************************************************************/

#ifndef Alignment_CTPPSTrackBased_CommonMethods_h
#define Alignment_CTPPSTrackBased_CommonMethods_h

class RPAlignmentCorrectionsData;
class AlignmentGeometry;

class CommonMethods
{

  public:
    CommonMethods() {}

    /// factors out the common shifts and rotations for every RP and saves these values as RPalignment
    /// (factored variable), the expanded alignments are created as a by-product
    static void FactorRPFromSensorCorrections(const RPAlignmentCorrectionsData &input,
      RPAlignmentCorrectionsData &expanded,
      RPAlignmentCorrectionsData &factored,
      const AlignmentGeometry &, bool equalWeights=false, unsigned int verbosity = 0);
};

#endif
