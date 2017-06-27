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
    input = cms.untracked.int32(10)
)

# random seeds

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
  ctppsFastLocalSimulation = cms.PSet(
    initialSeed = cms.untracked.uint32(81)
  )
)

# base geometry
process.load("Configuration.Geometry.geometry_CTPPS_alaTotem_cfi")

# misalignments
process.load("Geometry.VeryForwardGeometryBuilder.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.MisalignedFiles = cms.vstring()

# station simulation
process.load("Alignment.CTPPS.CTPPSFastLocalSimulation_cfi")
process.ctppsFastLocalSimulation.verbosity = 10

process.p = cms.Path(
    process.ctppsFastLocalSimulation
)
