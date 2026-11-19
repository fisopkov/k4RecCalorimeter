#ifndef RECFCCEECALORIMETER_SHOWERAXISFCCEE_H
#define RECFCCEECALORIMETER_SHOWERAXISFCCEE_H

// Key4HEP
#include "k4FWCore/DataHandle.h"
class IGeoSvc;

// Gaudi
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// DD4HEP
// #include "DDSegmentation/Segmentation.h"

// EDM4HEP
namespace edm4hep {
class ClusterCollection;
}

// DD4HEP
namespace dd4hep {
namespace DDSegmentation {
  class BitFieldCoder;
  class Segmentation;
} // namespace DDSegmentation
} // namespace dd4hep

/** @class ShowerAxisFCCee
 *
 *  Add to the cluster shape parameters the sum of the cluster cells energy and barycenter theta/phi coordinates
 *  per layer. The theta position is calculated with a log(E) weighting
 *
 *  @author Alexis Maloizel
 *  @author Giovanni Marchiori
 *  @author Tong Li
 */

class ShowerAxisFCCee : public Gaudi::Algorithm {

public:
  ShowerAxisFCCee(const std::string& name, ISvcLocator* svcLoc);
  StatusCode initialize();
  StatusCode execute(const EventContext& evtCtx) const;
  StatusCode finalize();

private:
    
  std::vector<size_t> m_thetaShapeParamIndices;
  std::vector<size_t> m_phiShapeParamIndices;

  /// Handle for input clusters
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_inClusters{"inClusters", Gaudi::DataHandle::Reader, this};
  /// Handle for output clusters
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_outClusters{"outClusters", Gaudi::DataHandle::Writer,
                                                                         this};

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  /// ID of the detectors
  Gaudi::Property<std::vector<int>> m_systemIDs{this, "systemIDs", {4}, "IDs of systems"};
  /// Name of the detectors (for the metadata)
  Gaudi::Property<std::vector<std::string>> m_detectorNames{
      this, "systemNames", {"EMB"}, "Names of the detectors, corresponding to systemIDs"};
  /// systemID of EMB, should be retrieved from mydetector.constantAsDouble("DetID_ECAL_Barrel")
  int systemID_EMB;
  /// Numbers of layers of the detectors
  Gaudi::Property<std::vector<size_t>> m_numLayers{this, "numLayers", {11}, "Numbers of layers of the systems"};
  /// Weights for each detector layer for theta position log-weighting
  
  /// Name of the detector readouts, corresponding to system IDs in m_systemIDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames{
      this,
      "readoutNames",
      {"ECalBarrelModuleThetaMerged"},
      "Names of the detector readouts, corresponding to systemIDs"};
  /// Name of the layer/cell field
  Gaudi::Property<std::vector<std::string>> m_layerFieldNames{
      this, "layerFieldNames", {"layer"}, "Identifiers of layers, corresponding to systemIDs"};
  Gaudi::Property<std::vector<std::string>> m_thetaFieldNames{
      this, "thetaFieldNames", {"theta"}, "Identifiers of theta, corresponding to systemIDs"};
  Gaudi::Property<std::vector<std::string>> m_moduleFieldNames{
      this, "moduleFieldNames", {"module"}, "Identifiers of module, corresponding to systemIDs"};
 
  

  // the number of grouped theta and phi cells
  //std::vector<int> nMergedThetaCells;
  //std::vector<int> nMergedModules;
  // the number of modules / phi cells
  //std::vector<int> nModules;

  /// Limit of debug printing
  mutable uint debugIter = 0;
  Gaudi::Property<uint> m_maxDebugPrint{this, "maxDebugPrint", 10,
                                        "maximum number of debug/warning messages from execute()"};

  void PrintDebugMessage(MsgStream stream, const std::string& text) const;
};

#endif /* RECFCCEECALORIMETER_SHOWERAXISFCCEE_H */
