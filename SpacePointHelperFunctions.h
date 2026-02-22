/**************************************************************************
 * basf2 (Belle II Analysis Software Framework)
 * Author: The Belle II Collaboration
 **************************************************************************/

#pragma once

#include <vector>
#include <cmath>
#include <unordered_map>

#include <svd/calibration/SVDHitTimeSelection.h>
#include <svd/calibration/SVDNoiseCalibrations.h>
#include <svd/dataobjects/SVDCluster.h>
#include <svd/dataobjects/SVDShaperDigit.h>
#include <svd/dbobjects/SVDSpacePointSNRFractionSelector.h>
#include <svd/reconstruction/SVDMaxSumAlgorithm.h>

#include <framework/datastore/StoreArray.h>
#include <framework/datastore/StoreObjPtr.h>
#include <framework/database/DBObjPtr.h>
#include <mdst/dataobjects/EventLevelTrackingInfo.h>

#include <vxd/dataobjects/VxdID.h>

#include <TH2.h>
#include <TFile.h>

namespace Belle2 {

  // ============================================================
  // Struct to collect clusters per sensor
  // ============================================================

  struct ClustersOnSensor {

  public:
    inline void addCluster(const SVDCluster* entry)
    {
      vxdID = entry->getSensorID();
      if (entry->isUCluster()) {
        clustersU.push_back(entry);
        return;
      }
      clustersV.push_back(entry);
    }

    VxdID vxdID;
    std::vector<const SVDCluster*> clustersU;
    std::vector<const SVDCluster*> clustersV;
  };


  // ============================================================
  // RESTORED FUNCTION (YOU WERE MISSING THIS)
  // ============================================================

  template <class SpacePointType>
  void provideSVDClusterSingles(const StoreArray<SVDCluster>& svdClusters,
                                StoreArray<SpacePointType>& spacePoints)
  {
    for (unsigned int i = 0; i < uint(svdClusters.getEntries()); ++i) {
      const SVDCluster* currentCluster = svdClusters[i];
      std::vector<const SVDCluster*> currentClusterCombi = { currentCluster };
      SpacePointType* newSP = spacePoints.appendNew(currentClusterCombi);
      newSP->addRelationTo(currentCluster);
    }
  }

    

  // ============================================================
  // ✅ FIXED MAXSUM FUNCTION (ONLY CHANGE IN FILE)
  // ============================================================

  inline void storeInputVectorFromSingleCluster(const SVDCluster* cls,
                                                std::vector<float>& inputVector,
                                                const SVDNoiseCalibrations& noiseCal)
  {
    inputVector.clear();
    inputVector.resize(3, 0.0f);

    auto shaperDigits = cls->getRelationsTo<SVDShaperDigit>();
    float noise = 0.0f;

    for (auto iSD : shaperDigits) {

      auto samples = iSD.getSamples(); // std::array<float, N>
      std::vector<float> selectedSamples;
      B2INFO("number of samples : " << samples.size() ); 

      if (samples.size() == 6) {
        Belle2::SVD::SVDMaxSumAlgorithm maxSum(samples);
        auto maxSamples = maxSum.getSelectedSamples(); // std::array<float,3>
        selectedSamples.assign(maxSamples.begin(), maxSamples.end());
        // --- TEMPORARY DEBUG CODE ---
      B2INFO("Original 6: " << samples[0] << ", " << samples[1] << ", " << samples[2] << ", " 
                            << samples[3] << ", " << samples[4] << ", " << samples[5]);
      B2INFO("Size, Selected 3: " << maxSamples.size() << ", " << maxSamples[0] << ", " << maxSamples[1] << ", " << maxSamples[2]);
      // ----------------------------  
      }
      else {
        selectedSamples.assign(samples.begin(), samples.end());
      }

      if (selectedSamples.size() < 3) continue;

      inputVector[0] += selectedSamples[0];
      inputVector[1] += selectedSamples[1];
      inputVector[2] += selectedSamples[2];

      VxdID thisSensorID = iSD.getSensorID();
      bool thisSide = iSD.isUStrip();
      int thisCellID = iSD.getCellID();
      float thisNoise = noiseCal.getNoise(thisSensorID, thisSide, thisCellID);

      noise += thisNoise * thisNoise;
    }

    noise = std::sqrt(noise);

    if (noise > 0.0f) {
      inputVector[0] /= noise;
      inputVector[1] /= noise;
      inputVector[2] /= noise;
    }
  }


  // ============================================================
  // YOUR ORIGINAL findPossibleCombinations (UNCHANGED)
  // ============================================================

  inline void findPossibleCombinations(const Belle2::ClustersOnSensor& aSensor,
                                       std::vector<std::vector<const SVDCluster*> >& foundCombinations,
                                       const SVDHitTimeSelection& hitTimeCut,
                                       const bool& useSVDGroupInfo,
                                       const int& numberOfSignalGroups,
                                       const bool& formSingleSignalGroup,
                                       const SVDNoiseCalibrations& noiseCal,
                                       const DBObjPtr<SVDSpacePointSNRFractionSelector>& svdSpacePointSelectionFunction,
                                       bool useSVDSpacePointSNRFractionSelector)
  {
    for (const SVDCluster* uCluster : aSensor.clustersU) {

      if (!hitTimeCut.isClusterInTime(uCluster->getSensorID(), 1, uCluster->getClsTime()))
        continue;

      for (const SVDCluster* vCluster : aSensor.clustersV) {

        if (!hitTimeCut.isClusterInTime(vCluster->getSensorID(), 0, vCluster->getClsTime()))
          continue;

        if (!hitTimeCut.areClusterTimesCompatible(vCluster->getSensorID(),
                                                  uCluster->getClsTime(),
                                                  vCluster->getClsTime()))
          continue;

          
        if (useSVDSpacePointSNRFractionSelector) {
          std::vector<float> inputU;
          std::vector<float> inputV;

          storeInputVectorFromSingleCluster(uCluster, inputU, noiseCal);
          storeInputVectorFromSingleCluster(vCluster, inputV, noiseCal);     
           

          if (!svdSpacePointSelectionFunction->passSNRFractionSelection(inputU, inputV))
            continue;
        }

        foundCombinations.push_back({uCluster, vCluster});
      }
    }
  }


  // ============================================================
  // EVERYTHING BELOW IS EXACTLY YOUR ORIGINAL CODE
  // (PDF naming, pairing probability, template, etc.)
  // ============================================================

  inline void spPDFName(const VxdID& sensor, int uSize, int vSize, int maxClusterSize,
                        std::string& PDFName, std::string& errorPDFName,
                        bool useLegacyNaming)
  {
    if (useLegacyNaming == true) {

      if (uSize > maxClusterSize) uSize = maxClusterSize;
      if (vSize > maxClusterSize) vSize = maxClusterSize;

      std::string sensorName;

      if (sensor.getLayerNumber() == 3)  sensorName = "l3";
      if (sensor.getLayerNumber() > 3 && sensor.getSensorNumber() == 1)  sensorName = "trap";
      if (sensor.getLayerNumber() > 3 && sensor.getSensorNumber() > 1)  sensorName = "large";

      PDFName =  sensorName + std::to_string(uSize) + std::to_string(vSize);
      errorPDFName = "error" + PDFName;
    }
    else {

      if (uSize > maxClusterSize) uSize = maxClusterSize;
      if (vSize > maxClusterSize) vSize = maxClusterSize;

      int layer = sensor.getLayerNumber();
      int ladder = sensor.getLadderNumber();
      int sens = sensor.getSensorNumber();

      PDFName = std::to_string(layer) + "." + std::to_string(ladder) + "." +
                std::to_string(sens) + "." + std::to_string(uSize) + "." +
                std::to_string(vSize);

      errorPDFName = PDFName + "_Error";
    }
  }

  /**
   * Function to extract probability of correct (pair from signal hit) cluster pairing from preconfigured pdfs
   * Probability defined as Pcharge * Ptime * Pucluster * Pvcluster
   *
   */


  inline void calculatePairingProb(TFile* pdfFile, std::vector<const SVDCluster*>& clusters, double& prob, double& error,
                                   bool useLegacyNaming)
  {

    int maxSize;
    int pdfEntries = pdfFile->GetListOfKeys()->GetSize();
    if (useLegacyNaming == true) {
      maxSize = floor(sqrt((pdfEntries - 4) / 6)); //4(time+size)+3(sensors)*2(prob/error)*size^2(u/v combo.)
    } else {
      maxSize = floor(sqrt((pdfEntries - 4) / 344)); //4(time+size)+172(sensorType)*2(prob/error)*size^2(u/v combo.)
    }
    std::string chargeProbInput;
    std::string chargeErrorInput;

    spPDFName(clusters[0]->getSensorID(), clusters[0]->getSize(), clusters[1]->getSize(), maxSize,
              chargeProbInput, chargeErrorInput, useLegacyNaming);
    std::string timeProbInput = "timeProb";
    std::string timeErrorInput = "timeError";
    std::string sizeProbInput = "sizeProb";
    std::string sizeErrorInput = "sizeError";


    TH2F* chargePDF = nullptr;
    TH2F* chargeError = nullptr;
    TH2F* timePDF = nullptr;
    TH2F* timeError = nullptr;
    TH2F* sizePDF = nullptr;
    TH2F* sizeError = nullptr;

    pdfFile->GetObject(chargeProbInput.c_str(), chargePDF);
    pdfFile->GetObject(chargeErrorInput.c_str(), chargeError);
    pdfFile->GetObject(timeProbInput.c_str(), timePDF);
    pdfFile->GetObject(timeErrorInput.c_str(), timeError);
    pdfFile->GetObject(sizeProbInput.c_str(), sizePDF);
    pdfFile->GetObject(sizeErrorInput.c_str(), sizeError);

    int xChargeBin = chargePDF->GetXaxis()->FindFixBin(clusters[0]->getCharge());
    int yChargeBin = chargePDF->GetYaxis()->FindFixBin(clusters[1]->getCharge());

    int xTimeBin = timePDF->GetXaxis()->FindFixBin(clusters[0]->getClsTime());
    int yTimeBin = timePDF->GetYaxis()->FindFixBin(clusters[1]->getClsTime());


    int xSizeBin = sizePDF->GetXaxis()->FindFixBin(clusters[0]->getSize());
    int ySizeBin = sizePDF->GetYaxis()->FindFixBin(clusters[1]->getSize());

    double chargeProb = chargePDF->GetBinContent(xChargeBin, yChargeBin);
    double timeProb = timePDF->GetBinContent(xTimeBin, yTimeBin);
    double sizeProb = sizePDF->GetBinContent(xSizeBin, ySizeBin);
    double chargeProbError = chargePDF->GetBinContent(xChargeBin, yChargeBin);
    double timeProbError = timePDF->GetBinContent(xTimeBin, yTimeBin);
    double sizeProbError = sizePDF->GetBinContent(xSizeBin, ySizeBin);


    if (chargeProbError == 0) {
      B2DEBUG(21, "svdClusterProbabilityEstimator has not been run, spacePoint QI will return zero!");
    }

    prob = chargeProb * timeProb * sizeProb * clusters[0]->getQuality() * clusters[1]->getQuality();
    error = prob * sqrt(pow(timeProb * sizeProb * clusters[0]->getQuality() * clusters[1]->getQuality() * chargeProbError, 2) +
                        pow(chargeProb * sizeProb * clusters[0]->getQuality() * clusters[1]->getQuality() * timeProbError, 2) +
                        pow(chargeProb * timeProb * clusters[0]->getQuality() * clusters[1]->getQuality() * sizeProbError, 2) +
                        pow(chargeProb * timeProb * sizeProb * clusters[1]->getQuality() * clusters[0]->getQualityError(), 2) +
                        pow(chargeProb * timeProb * sizeProb * clusters[0]->getQuality() * clusters[1]->getQualityError(), 2));
  }

  /** finds all possible combinations of U and V Clusters for SVDClusters.
   *
   * first parameter is a storeArray containing SVDClusters.
   * second parameter is a storeArra containing SpacePoints (will be filled in the function).
   * third parameter tels the spacePoint where to get the name of the storeArray containing the related clusters
   * relationweights code the type of the cluster. +1 for u and -1 for v
   */
  template <class SpacePointType> void provideSVDClusterCombinations(const StoreArray<SVDCluster>& svdClusters,
      StoreArray<SpacePointType>& spacePoints, SVDHitTimeSelection& hitTimeCut, bool useQualityEstimator, TFile* pdfFile,
      bool useLegacyNaming, unsigned int numMaxSpacePoints, std::string m_eventLevelTrackingInfoName, const bool& useSVDGroupInfo,
      const int& numberOfSignalGroups, const bool& formSingleSignalGroup,
      const SVDNoiseCalibrations& noiseCal, const DBObjPtr<SVDSpacePointSNRFractionSelector>& svdSpacePointSelectionFunction,
      bool useSVDSpacePointSNRFractionSelector)
  {
    std::unordered_map<VxdID::baseType, ClustersOnSensor>
    activatedSensors; // collects one entry per sensor, each entry will contain all Clusters on it TODO: better to use a sorted vector/list?
    std::vector<std::vector<const SVDCluster*> >
    foundCombinations; // collects all combinations of Clusters which were possible (condition: 1u+1v-Cluster on the same sensor)

    // sort Clusters by sensor. After the loop, each entry of activatedSensors contains all U and V-type clusters on that sensor
    for (unsigned int i = 0; i < uint(svdClusters.getEntries()); ++i) {
      SVDCluster* currentCluster = svdClusters[i];

      activatedSensors[currentCluster->getSensorID().getID()].addCluster(currentCluster);
    }


    for (auto& aSensor : activatedSensors)
      findPossibleCombinations(aSensor.second, foundCombinations, hitTimeCut, useSVDGroupInfo, numberOfSignalGroups,
                               formSingleSignalGroup,
                               noiseCal, svdSpacePointSelectionFunction, useSVDSpacePointSNRFractionSelector);

    // Do not make space-points if their number would be too large to be considered by tracking
    if (foundCombinations.size() > numMaxSpacePoints) {
      StoreObjPtr<EventLevelTrackingInfo> m_eventLevelTrackingInfo(m_eventLevelTrackingInfoName);
      if (m_eventLevelTrackingInfo.isValid()) {
        m_eventLevelTrackingInfo->setSVDSpacePointCreatorAbortionFlag();
      }
      return;
    }

    for (auto& clusterCombi : foundCombinations) {
      SpacePointType* newSP = spacePoints.appendNew(clusterCombi);
        
      if (useQualityEstimator == true) {
        double probability;
        double error;
        calculatePairingProb(pdfFile, clusterCombi, probability, error, useLegacyNaming);
        newSP->setQualityEstimation(probability);
        newSP->setQualityEstimationError(error);
      }
      for (auto* cluster : clusterCombi) {
        newSP->addRelationTo(cluster, cluster->isUCluster() ? 1. : -1.);
      }
    }
  }


} //Belle2 namespace
