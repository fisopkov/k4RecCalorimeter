#include "ShowerAxisFCCee.h"

// k4FWCore
#include "k4FWCore/MetadataUtils.h"
#include "k4Interface/IGeoSvc.h"

// edm4hep
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/Constants.h"

// DD4hep
#include "DD4hep/Detector.h"

// k4geo
#include "detectorSegmentations/FCCSWGridModuleThetaMerged_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"

// ROOT
//#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
//#include "TVector3.h"
#include <Eigen/Dense>
#include <algorithm>
#include <unordered_map>

DECLARE_COMPONENT(ShowerAxisFCCee)

ShowerAxisFCCee::ShowerAxisFCCee(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc), m_geoSvc("GeoSvc", "ShowerAxisFCCee") {
  declareProperty("inClusters", m_inClusters, "Input clusters");
  declareProperty("outClusters", m_outClusters, "Output clusters");
}

void ShowerAxisFCCee::PrintDebugMessage(MsgStream stream, const std::string& text) const {
  if (debugIter < m_maxDebugPrint) {
    stream << text << endmsg;
    debugIter++;
  } else if (debugIter == m_maxDebugPrint) {
    stream << "Maximum number of messages reached, suppressing further output" << endmsg;
    debugIter++;
  } else
    debugIter++;
}

StatusCode ShowerAxisFCCee::initialize() {
  {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure())
      return sc;
  }

  // check if readouts exist
  for (size_t k = 0; k < m_readoutNames.size(); k++) {
    std::string readoutName = m_readoutNames[k];
    if (m_geoSvc->getDetector()->readouts().find(readoutName) == m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << readoutName << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Check if readout related variables have the same size
  if (m_systemIDs.size() != m_readoutNames.size()) {
    error() << "Sizes of the systemIDs vector and readoutNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_numLayers.size()) {
    error() << "Sizes of systemIDs vector and numLayers vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_layerFieldNames.size()) {
    error() << "Sizes of systemIDs vector and layerFieldNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_thetaFieldNames.size()) {
    error() << "Sizes of systemIDs vector and thetaFieldNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_moduleFieldNames.size()) {
    error() << "Sizes of systemIDs vector and moduleFieldNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_detectorNames.size()) {
    error() << "Sizes of systemIDs vector and detectorNames vector do not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }

  // retrieve systemID of EMB from "DetID_ECAL_Barrel"
  systemID_EMB = m_geoSvc->getDetector()->constantAsDouble("DetID_ECAL_Barrel");

  // retrieve some information from the segmentation for later use
  // - number of modules/phi bins => needed to account for module/phi periodicity
  // - number of merged modules and theta cells per layer => needed to check if two cells are neighbours
 

  // initialise the list of metadata for the clusters
  // append to the metadata of the input clusters (if any)
  //auto showerShapeDecorations = k4FWCore::getCollectionParameter<std::vector<std::string>>(
  //                                  m_inClusters.objKey(), edm4hep::labels::ShapeParameterNames, this)
  //                                  .value_or(std::vector<std::string>{});

  

  //k4FWCore::putCollectionParameter(m_outClusters.objKey(), edm4hep::labels::ShapeParameterNames, showerShapeDecorations,
  //                                 this);

  auto shapeParameterNames =
    k4FWCore::getCollectionParameter<std::vector<std::string>>(
        m_inClusters.objKey(),
        edm4hep::labels::ShapeParameterNames,
        this
    ).value_or(std::vector<std::string>{});


    m_thetaShapeParamIndices.clear();
m_phiShapeParamIndices.clear();

for (unsigned layer = 0; layer < m_numLayers[0]; ++layer) {

    std::string thetaName = Form("theta_EMB_layer_%d", layer);
    std::string phiName   = Form("phi_EMB_layer_%d", layer);

    auto thetaIt = std::find(shapeParameterNames.begin(),
                             shapeParameterNames.end(), thetaName);

    auto phiIt = std::find(shapeParameterNames.begin(),
                           shapeParameterNames.end(), phiName);

    if (thetaIt == shapeParameterNames.end() ||
        phiIt == shapeParameterNames.end()) {
        error() << "Cannot find theta/phi shape parameters for layer "
                << layer << endmsg;
        return StatusCode::FAILURE;
    }

    m_thetaShapeParamIndices.push_back(
        std::distance(shapeParameterNames.begin(), thetaIt));

    m_phiShapeParamIndices.push_back(
        std::distance(shapeParameterNames.begin(), phiIt));
}

  return StatusCode::SUCCESS;
}

// for all cells in a certain layer: vec A is theta_id, vec B is E_cell
// make a 1D projection to theta: sum up the E_cell over modules with the same theta_id
// then sort the theta vec (for finding theta neighbors), update the E vec simultaneously
// return a pair of vectors: theta_id (in ascending order) and E (summed over modules)


StatusCode ShowerAxisFCCee::finalize() { return Gaudi::Algorithm::finalize(); }

StatusCode ShowerAxisFCCee::execute([[maybe_unused]] const EventContext& evtCtx) const {
  // get the input collection with clusters
  const edm4hep::ClusterCollection* inClusters = m_inClusters.get();
  // create the new output collection
  edm4hep::ClusterCollection* outClusters = m_outClusters.createAndPut();

  // total number of layers
  size_t numLayersTotal = 0;
  for (size_t i = 0; i < m_numLayers.size(); ++i) {
    numLayersTotal += m_numLayers[i];
  }
  
  const double maxAxisAngle = 20.0 * dd4hep::deg;

  // loop over the clusters
  for (const auto& cluster : *inClusters) {
   
    // clone original cluster
    auto newCluster = cluster.clone();
    std::unordered_map<dd4hep::DDSegmentation::CellID, Eigen::Vector3d> correctedPositions;
    //outClusters->push_back(newCluster);
    //std::cout << "  Cluster with energy " << cluster.getEnergy() << " and " << cluster.hits_size() << " cells." << std::endl;

    
    // loop over all cells to:
    // - calculate the cluster invariant mass
    // - calculate the energy deposited in each layer
    // - find out the energy of the cells with largest energy in each layer
    // - find out if cluster is around -pi..pi transition and/or max module .. 0 transition

    
    unsigned int nCells(0);
    
    int startPositionToFill = 0;
    // loop over each system/readout
    
    for (size_t k = 0; k < m_readoutNames.size(); k++) {

      if (k > 0)
        startPositionToFill += m_numLayers[k - 1];
      double E(0.0);
      double deltaPhi = 0.;
      double deltaTheta = 0.;
      double angle = 0.;

      
      Eigen::Vector3d v3clW = Eigen::Vector3d::Zero();

      Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
      Eigen::Vector3d prAxis = Eigen::Vector3d::Zero();
      Eigen::Vector3d showerCenter = Eigen::Vector3d::Zero();
      Eigen::Vector3d showerAxis = Eigen::Vector3d::Zero();

      

      int systemID = m_systemIDs[k];
     
      std::string layerField = m_layerFieldNames[k];
      
      std::string readoutName = m_readoutNames[k];
      dd4hep::DDSegmentation::BitFieldCoder* decoder = m_geoSvc->getDetector()->readout(readoutName).idSpec().decoder();


      const auto& shapeParameters = cluster.getShapeParameters();

      std::vector<double> thetaLayer(m_numLayers[0], 0.);
      std::vector<double> phiLayer(m_numLayers[0], 0.);

      for (unsigned layer = 0; layer < m_numLayers[0]; ++layer) {

        thetaLayer[layer] = shapeParameters[m_thetaShapeParamIndices[layer]];
        phiLayer[layer] = shapeParameters[m_phiShapeParamIndices[layer]];
      }


      // loop over the cells to get E_layer, E_cluster, max and min phi of cluster cells
      for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++) {
        dd4hep::DDSegmentation::CellID cID = cell->getCellID();

        // skip cells with wrong system ID
        int sysId = decoder->get(cID, "system");
        if (sysId != systemID)
          continue;

        // retrieve layer
        uint layer = decoder->get(cID, layerField);

        // retrieve cell energy, update sum of cell energies per layer, total energy, and max cell energy
        double eCell = cell->getEnergy();
        E += eCell;


        Eigen::Vector3d vOrg(cell->getPosition().x,
                          cell->getPosition().y,
                          cell->getPosition().z);

        double r = vOrg.norm();

        double theta = std::atan2(std::sqrt(vOrg.x() * vOrg.x() + vOrg.y() * vOrg.y()), vOrg.z());
        double phi = std::atan2(vOrg.y(), vOrg.x());

        // Get average angles of the corresponding layer
        double thetaLayerValue = thetaLayer[layer];
        double phiLayerValue = phiLayer[layer];

        // Calculate differences
        double deltaTheta = theta - thetaLayerValue;
        double deltaPhi = std::remainder(phi - phiLayerValue, 2.0 * M_PI);

        // Scale according to cell width
        //double correctionScale = (layer == 2) ? 0.5 : 1.0;

        // Apply scaled differences
        theta = thetaLayerValue; //( correctionScale * deltaTheta);
        //phi = phiLayerValue; //( correctionScale * deltaPhi);

        // Normalize angles
        theta = std::clamp(theta, 0.0, M_PI);
        phi = std::remainder(phi, 2.0 * M_PI);

        Eigen::Vector3d v(r * std::sin(theta) * std::cos(phi),
                             r * std::sin(theta) * std::sin(phi),
                             r * std::cos(theta));


       

        // compute phi of cell to determine min/max phi of cluster

        //v3clW += cell->getEnergy() * Eigen::Vector3d(cell->getPosition().x,
        //                                            cell->getPosition().y,
        //                                            cell->getPosition().z);
       
        v3clW += cell->getEnergy() * v;                                              
      
        // Store corrected position and corresponding cell ID
        correctedPositions[cID] = v;
 
        nCells++;
        
      } // end of loop over cells
    
        

    if (E == 0.0) {
      warning() << "No cells found for readout " << readoutName
              << " in this cluster." << endmsg;
      continue;
    }


    v3clW /= E;


    showerCenter = v3clW;
    showerAxis = showerCenter.normalized();

    double w = 0;
  

    for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++) { //fs
      
      dd4hep::DDSegmentation::CellID cID = cell->getCellID();

      int sysId = decoder->get(cID, "system");
      if (sysId != systemID)
        continue;

      auto it = correctedPositions.find(cell->getCellID());

      if (it == correctedPositions.end())
          continue;

      const Eigen::Vector3d& v = it->second;  

      if (nCells > 2){
        
        

        double e2 = cell->getEnergy() * cell->getEnergy(); //fs
      
        /*C(0,0) += e2 * (cell->getPosition().x-v3clW[0]) * (cell->getPosition().x-v3clW[0]);
        C(1,0) += e2 * (cell->getPosition().x-v3clW[0]) * (cell->getPosition().y-v3clW[1]);
        C(2,0) += e2 * (cell->getPosition().x-v3clW[0]) * (cell->getPosition().z-v3clW[2]);

        C(1,1) += e2 * (cell->getPosition().y-v3clW[1]) * (cell->getPosition().y-v3clW[1]);
        C(2,1) += e2 * (cell->getPosition().y-v3clW[1]) * (cell->getPosition().z-v3clW[2]);

        C(2,2) += e2 * (cell->getPosition().z-v3clW[2]) * (cell->getPosition().z-v3clW[2]);*/

        C(0,0) += e2 * (v.x()-v3clW[0]) * (v.x()-v3clW[0]);
        C(1,0) += e2 * (v.x()-v3clW[0]) * (v.y()-v3clW[1]);
        C(2,0) += e2 * (v.x()-v3clW[0]) * (v.z()-v3clW[2]);

        C(1,1) += e2 * (v.y()-v3clW[1]) * (v.y()-v3clW[1]);
        C(2,1) += e2 * (v.y()-v3clW[1]) * (v.z()-v3clW[2]);

        C(2,2) += e2 * (v.z()-v3clW[2]) * (v.z()-v3clW[2]);
        w += e2;
      }
    }
    C/=(w != 0 ? w : 1.0);

    
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(C);
    if (eigensolver.info() != Eigen::Success) {
        debug() << "Failed to compute Eigenvalues -> Can't determine shower axis" << endmsg;
    }
    else {
         // don't use the principal axes if at least one of the 3 
         // diagonal elements is 0
 
         const Eigen::Vector3d& S=eigensolver.eigenvalues(); 
         const Eigen::Matrix3d& U=eigensolver.eigenvectors();
 
         const double epsilon = 1.E-6;
 
         if ( std::abs(S[0]) >= epsilon && std::abs(S[1]) >= epsilon && std::abs(S[2]) >= epsilon ) { 
        
           prAxis = showerAxis;
           int iEigen = -1;
           
           for (int i=0;i<3;i++) {
              Eigen::Vector3d tmpAxis=U.col(i);
 
            // calculate the angle          
              double cosAngle = tmpAxis.normalized().dot(showerAxis.normalized());
              cosAngle = std::clamp(cosAngle, -1.0, 1.0);

              double tmpAngle = std::acos(cosAngle);


              if ( tmpAngle > 90*dd4hep::deg ) { 
                tmpAngle = 180*dd4hep::deg - tmpAngle;
                tmpAxis = -tmpAxis;
              }
         
              if ( iEigen == -1 || tmpAngle < angle ) {
                iEigen = i;
                angle = tmpAngle;
                prAxis = tmpAxis;
              }
           }//end for loop     
         
           // calculate theta and phi angle differences
           deltaPhi = std::remainder(std::atan2(showerAxis.y(), showerAxis.x()) -std::atan2(prAxis.y(), prAxis.x()),2.0 * M_PI);
           deltaTheta = std::acos(showerAxis.normalized().z()) - std::acos(prAxis.normalized().z());
         
           
           // check the angle
           
           if ( angle < maxAxisAngle ) {
              showerAxis = prAxis;
           }
           else{ 
                   debug() << "     principal Direction (" << prAxis[0] << ", " 
                   << prAxis[1] << ", " << prAxis[2] << ") deviates MORE than " << maxAxisAngle*(1./dd4hep::deg)
                   << " deg from IP-to-ClusterCenter-axis (" << showerAxis[0] << ", "
                   << showerAxis[1] << ", " << showerAxis[2] << ")"
                  <<endmsg;
           }      
         }//end if std::abs(S)<epsilon
         else {
          debug() << "Eigenvalues close to 0, do not use principal axis" << endmsg;
         }
       }//end got eigenvalues else

 
    // Store shower center in the cluster position
    newCluster.setPosition({
      static_cast<float>(showerCenter.x()),
      static_cast<float>(showerCenter.y()),
      static_cast<float>(showerCenter.z())
    });

    // Normalize the principal axis
    if (showerAxis.norm() > 0.0) {
      showerAxis.normalize();

      // Convert the principal axis to spherical angles
      
      const double phi = std::atan2(showerAxis.y(), showerAxis.x());

      
      const double theta = std::atan2(
        std::sqrt(showerAxis.x() * showerAxis.x() + showerAxis.y() * showerAxis.y()), showerAxis.z());

      // Store the direction in the existing cluster fields
      newCluster.setIPhi(static_cast<float>(phi));
      newCluster.setITheta(static_cast<float>(theta));
  
    } else {
        warning() << "Principal axis has zero norm, cannot normalize." << endmsg;
    }

    // Add the modified cluster to the output collection
    outClusters->push_back(newCluster);

   } // end of loop over system / readout 
  } // end of loop over clusters

  return StatusCode::SUCCESS;
}
