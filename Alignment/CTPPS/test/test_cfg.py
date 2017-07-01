import FWCore.ParameterSet.Config as cms

process = cms.Process("FastAlignmentSimulation")

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
  statistics = cms.untracked.vstring(),
  destinations = cms.untracked.vstring('cout'),
  cout = cms.untracked.PSet(
    threshold = cms.untracked.string('WARNING')
  )
)

# set number of events
process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2)
)

# random seeds

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
  ctppsFastLocalSimulation = cms.PSet(
    initialSeed = cms.untracked.uint32(81)
  )
)

# base geometry
process.load("geometry_CTPPS_alaTotem_RECO_cfi")

# misalignments
process.load("Geometry.VeryForwardGeometryBuilder.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.MisalignedFiles = cms.vstring("./alignment.xml")

# geometry printer
process.geomInfo = cms.EDAnalyzer("GeometryInfoModule",
  geometryType = cms.untracked.string("misaligned"),
  printRPInfo = cms.untracked.bool(True),
  printSensorInfo = cms.untracked.bool(True),
  printMeanSensorInfo = cms.untracked.bool(False)
)

# station simulation
process.load("Alignment.CTPPS.ctppsFastLocalSimulation_cfi")
process.ctppsFastLocalSimulation.verbosity = 10
process.ctppsFastLocalSimulation.roundToPitch = False # TODO: eventually change to True
process.ctppsFastLocalSimulation.z0 = 210000
#process.ctppsFastLocalSimulation.angular_distribution.x_width = 0
#process.ctppsFastLocalSimulation.angular_distribution.y_width = 0

# alignment
process.load("Alignment.CTPPS.ctppsStraightTrackAligner_cfi")
process.ctppsStraightTrackAligner.verbosity = 10
process.ctppsStraightTrackAligner.rpIds = cms.vuint32(103) # TODO: eventually add 116 and 123
process.ctppsStraightTrackAligner.z0 = process.ctppsFastLocalSimulation.z0
process.ctppsStraightTrackAligner.maxResidualToSigma = 10
process.ctppsStraightTrackAligner.chiSqPerNdfCut = 100
process.ctppsStraightTrackAligner.algorithms = cms.vstring("Jan")

process.ctppsStraightTrackAligner.constraintsType = cms.string("fixedDetectors")
process.ctppsStraightTrackAligner.fixedDetectorsConstraints = cms.PSet(
      ShR1 = cms.PSet(
        ids = cms.vuint32(),
        values = cms.vdouble()
      ),
      ShR2 = cms.PSet(
        ids = cms.vuint32(1998061568, 1998094336, 1998323712, 1998356480),
        values = cms.vdouble(0, 0, 0, 0)
      )
    )

process.eca = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.geomInfo
    * process.ctppsFastLocalSimulation
    * process.ctppsStraightTrackAligner
)
