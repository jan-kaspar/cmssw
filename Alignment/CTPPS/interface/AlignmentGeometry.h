/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#ifndef _AlignmentGeometry_h_
#define _AlignmentGeometry_h_

#include <map>
#include <set>
#include <string>

//----------------------------------------------------------------------------------------------------

/**
 *\brief A structure to hold relevant geometrical information about one detector/sensor.
 **/
struct DetGeometry
{
  double z;                     ///< z postion at detector centre; mm

  double sx, sy;                ///< detector nominal shift = detector center in global coordinates; in mm

  double d1x, d1y;              ///< local (1, 0, 0) in global coordinates
  double d2x, d2y;              ///< local (0, 1, 0) in global coordinates

  double s1;                    ///< detector nominal shift in direction d1
  double s2;                    ///< detector nominal shift in direction d2

  unsigned int matrixIndex;     ///< index (0 ... AlignmentGeometry::Detectors()) within a S matrix block (for detector-related quantities)

  bool isU;                     ///< only relevant for strips: true for U detectors, false for V detectors
                                ///< global U, V frame is used - that matches with u, v frame of the 1200 detector

  DetGeometry(double _z = 0., double _sx = 0., double _sy = 0.,
    double _d1x = 0., double _d1y = 0.,
    double _d2x = 0., double _d2y = 0.,
    bool _isU = false) : z(_z), sx(_sx), sy(_sy), d1x(_d1x), d1y(_d1y), d2x(_d2x), d2y(_d2y), matrixIndex(0), isU(_isU)
  {
    s1 = sx*d1x + sy*d1y;
    s2 = sx*d2x + sy*d2y;
  }
};

//----------------------------------------------------------------------------------------------------

/**
 *\brief A collection of geometrical information.
 * Map: (decimal) detector ID --> DetGeometry
 **/
class AlignmentGeometry : public std::map<unsigned int, DetGeometry>
{
  protected:
    std::set<unsigned int> rps;

  public:
    /// a characteristic z in mm
    double z0;

    /// puts an element to the map
    void Insert(unsigned int id, const DetGeometry &g);

    /// returns the number of detectors in the collection
    unsigned int Detectors()
    {
      return size();
    }

    /// returns the number of RPs in the collection
    unsigned int RPs()
    {
      return rps.size();
    }
    
    /// returns detector id corresponding to the given matrix index
    unsigned int MatrixIndexToDetId(unsigned int) const;

    /// returns reference the the geometry of the detector with the given matrix index
    const_iterator FindByMatrixIndex(unsigned int) const;
    
    /// check whether the sensor Id is valid (present in the map)
    bool ValidSensorId(unsigned int id) const
    {
      return (find(id) != end());
    }

    /// check whether the RP Id is valid (present in the set)
    bool ValidRPId(unsigned int id) const
    {
      return (rps.find(id) != rps.end());
    }

    /// Prints the geometry.
    void Print() const;

    /// loads geometry from a text file of 5 columns:
    /// id | center x, y, z (all in mm) | read-out direction x projection, y projection
    void LoadFromFile(const std::string &filename);
};

#endif
