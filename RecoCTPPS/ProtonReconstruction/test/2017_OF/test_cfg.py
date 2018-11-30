import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("CTPPSTestProtonReconstruction", eras.ctpps_2016)

process.load("test_base_cff")
from test_base_cff import *

# run 302654
process.source.fileNames = cms.untracked.vstring("/store/data/Run2017D/DoubleEG/MINIAOD/17Nov2017-v1/70000/90E543DB-42D5-E711-9F5B-0242AC1C0501.root")

#process.source.fileNames = cms.untracked.vstring("root://eostotem.cern.ch//eos/totem/data/ctpps/reconstruction/2017/preTS2_alignment_data/version4/fill5839_xangle140_DoubleEG.root")
#process.source.fileNames = cms.untracked.vstring("root://eostotem.cern.ch//eos/totem/data/ctpps/reconstruction/2017/preTS2_alignment_data/version4/fill6091_xangle140_DoubleEG.root")
#process.source.fileNames = cms.untracked.vstring("root://eostotem.cern.ch//eos/totem/data/ctpps/reconstruction/2017/preTS2_alignment_data/version4/fill6192_xangle140_DoubleEG.root")

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(100000)
)

UseCrossingAngle140(process.source)
