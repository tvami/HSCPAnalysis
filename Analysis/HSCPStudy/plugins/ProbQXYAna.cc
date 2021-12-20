// -*- C++ -*-
//
// Package:    ProbQXYAna
// Class:      ProbQXYAna
//
/**\class ProbQXYAna  ProbQXYAna.cc RecoTracker/DeDx/plugins/ProbQXYAna.cc

 Description: SiPixel Charge and shape probabilities combined for tracks

*/
//
// Original Author:  Tamas Almos Vami
//         Created:  Mon Nov 17 14:09:02 CEST 2021
//

// system include files
#include <memory>

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/TrackReco/interface/SiPixelTrackProbQXY.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TTree.h"

//
// class declaration
//

class ProbQXYAna : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ProbQXYAna(const edm::ParameterSet&);
  ~ProbQXYAna() override = default;
  float combineProbs(float probOnTrackWMulti, int numRecHits) const;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  const bool debugFlag_;
  const std::string pixelCPE_;
  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<reco::DeDxHitInfoAss> gt2dedxHitInfo_;
  const edm::EDPutTokenT<edm::ValueMap<reco::SiPixelTrackProbQXY>> putProbQXYVMToken_;
  const edm::EDPutTokenT<edm::ValueMap<reco::SiPixelTrackProbQXY>> putProbQXYNoLayer1VMToken_;
  TTree* smalltree;
  int numTrack; 
  float tree_track_probQ[100000];
  float tree_track_probQNew[100000];
};

using namespace reco;
using namespace std;
using namespace edm;

ProbQXYAna::ProbQXYAna(const edm::ParameterSet& iConfig)
    : debugFlag_(iConfig.getUntrackedParameter<bool>("debug", false)),
      pixelCPE_(iConfig.getParameter<std::string>("pixelCPE")),
      trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      gt2dedxHitInfo_(consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("dedxInfos"))) {
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  smalltree = fs->make<TTree>("ttree", "ttree");
  smalltree->Branch("numTrack", &numTrack);
  smalltree->Branch("track_probQ", tree_track_probQ, "track_probQ[numTrack]/F");
  smalltree->Branch("track_probQNew", tree_track_probQNew, "track_probQNew[numTrack]/F");
}

void ProbQXYAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Retrieve track collection from the event
  auto trackCollectionHandle = iEvent.getHandle(trackToken_);
  const TrackCollection& trackCollection(*trackCollectionHandle.product());
  numTrack = 0;
  int numTrackWSmallProbQ = 0;

  auto const& gt2dedxHitInfo = iEvent.get(gt2dedxHitInfo_);

  // Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> TopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(TopoHandle);
  const TrackerTopology* tTopo = TopoHandle.product();

  edm::ESHandle<PixelClusterParameterEstimator> pixelCPE;
  iSetup.get<TkPixelCPERecord>().get(pixelCPE_, pixelCPE);

  edm::ESHandle<TrackerGeometry> tkGeometry;
  iSetup.get<TrackerDigiGeometryRecord>().get(tkGeometry);

  // Creates the output collection
  auto resultSiPixelTrackProbQXYColl = std::make_unique<reco::SiPixelTrackProbQXYCollection>();
  auto resultSiPixelTrackProbQXYNoLayer1Coll = std::make_unique<reco::SiPixelTrackProbQXYCollection>();

  // Loop through the tracks
  for (const auto& track : trackCollection) {
    if (!(track.pt() > 10)) continue;
    reco::TrackRef trackRef(trackCollectionHandle, numTrack);
    numTrack++;
    float probQonTrack = 0.0;
    float probXYonTrack = 0.0;
    float probQonTrackNoLayer1 = 0.0;
    float probXYonTrackNoLayer1 = 0.0;
    int numRecHits = 0;
    int numRecHitsNoLayer1 = 0;
    float probQonTrackWMulti = 1;
    float probXYonTrackWMulti = 1;
    float probQonTrackWMultiNoLayer1 = 1;
    float probXYonTrackWMultiNoLayer1 = 1;

    auto const& dedxRef = gt2dedxHitInfo[trackRef];
    const reco::DeDxHitInfo* dedxHits = nullptr;
    if (!dedxRef.isNull()) {
      dedxHits = &(*dedxRef);
    } else {
      LogDebug("ProbQXYAna") << " dedxRef.isNull(), skipping it";
      continue;
    }

    if (debugFlag_) {
      LogPrint("ProbQXYAna") << "  -----------------------------------------------";
      LogPrint("ProbQXYAna") << "  For track " << numTrack << " dedxinfoRef is "
                             << (!dedxRef.isNull());
      if (!dedxRef.isNull()) {
        LogPrint("ProbQXYAna") << " Track pt : " << track.pt() << " eta: " << track.eta() << " phi: " << track.phi() << " dedx: size "
                               << dedxRef->size();
      }
    }

    // Loop through the rechits on the given track
    for (unsigned int iHit = 0; iHit < dedxHits->size(); iHit++) {
      if (!dedxRef.isNull()) {
        DetId dedxId = dedxRef->detId(iHit);
        auto const* dedxClu = dedxRef->pixelCluster(iHit);
        if (dedxClu == nullptr) {
          LogDebug("ProbQXYAna") << "  >>>>>> FAILED to get dedxClu";
        } else {
          LogDebug("ProbQXYAna") << "  >>>>>> dedx " << dedxClu->x() << " " << dedxClu->y() << " " << dedxClu->charge()
                                 << " " << dedxClu->size();
          const GeomDetUnit& geomDet = *tkGeometry->idToDetUnit(dedxId);
          //lazy option, no propagation
          LocalVector lv = geomDet.toLocal(GlobalVector(track.px(), track.py(), track.pz()));
          auto qual_2 = std::get<2>(pixelCPE->getParameters(
              *dedxClu, geomDet, LocalTrajectoryParameters(dedxRef->pos(iHit), lv, track.charge())));
          float probQ = SiPixelRecHitQuality::thePacking.probabilityQ(qual_2);
          float probXY = SiPixelRecHitQuality::thePacking.probabilityXY(qual_2);
          if (probQ == 0) {
            continue;  // if any of the rechits have zero probQ, skip them
          }
          numRecHits++;

          if (debugFlag_) {
            LogPrint("ProbQXYAna") << "    >>>> For rechit: " << numRecHits << " ProbQ = " << probQ;
          }

          // Calculate alpha term needed for the combination
          probQonTrackWMulti *= probQ;
          probXYonTrackWMulti *= probXY;

        }  // end of if for dedxClu's existence
      }    // end of if for DeDx being non null

      // Have a separate variable that excludes Layer 1
      // Layer 1 was very noisy in 2017/2018
      /*  if ((dedxId.subdetId() == PixelSubdetector::PixelEndcap) ||
         (dedxId.subdetId() == PixelSubdetector::PixelBarrel &&
          tTopo->pxbLayer(dedxId) != 1)) {
        float probQNoLayer1 = 0;
        float probXYNoLayer1 = 0;  //pixhit->probabilityXY();
        if (probQNoLayer1 > 0.f) {  // only save the non-zero rechits
          numRecHitsNoLayer1++;
          // Calculate alpha term needed for the combination
          probQonTrackWMultiNoLayer1 *= probQNoLayer1;
          probXYonTrackWMultiNoLayer1 *= probXYNoLayer1;
          if(debugFlag_) {
            LogPrint("ProbQXYAna") << "    >>>> For rechit excluding Layer 1: " << numRecHitsNoLayer1 << " ProbQ = " << probQNoLayer1;
      }

        }
      }*/
    }  // end looping on the rechits
    probQonTrack = combineProbs(probQonTrackWMulti, numRecHits);
    probXYonTrack = combineProbs(probXYonTrackWMulti, numRecHits);
    LogPrint("ProbQXYAna") << "probQonTrack: " << probQonTrack;
    LogDebug("ProbQXYAna") << probXYonTrack;
   tree_track_probQ[numTrack] = probQonTrack;
   tree_track_probQNew[numTrack] = probQonTrack; // probQonTrackNew;
    // probQonTrackNoLayer1 = combineProbs(probQonTrackWMultiNoLayer1, numRecHitsNoLayer1);
    // probXYonTrackNoLayer1 = combineProbs(probXYonTrackWMultiNoLayer1, numRecHitsNoLayer1);
  }  // end loop on track collection
  smalltree->Fill();
}

float ProbQXYAna::combineProbs(float probOnTrackWMulti, int numRecHits) const {
  float logprobOnTrackWMulti = probOnTrackWMulti > 0 ? log(probOnTrackWMulti) : 0;
  float factQ = -logprobOnTrackWMulti;
  float probOnTrackTerm = 0.f;

  if (numRecHits == 1) {
    probOnTrackTerm = 1.f;
  } else if (numRecHits > 1) {
    probOnTrackTerm = 1.f + factQ;
    for (int iTkRh = 2; iTkRh < numRecHits; ++iTkRh) {
      factQ *= -logprobOnTrackWMulti / float(iTkRh);
      probOnTrackTerm += factQ;
    }
  }
  float probOnTrack = probOnTrackWMulti * probOnTrackTerm;

  //  LogPrint("SiPixelTrackProbQXYProducer")
  //      << "  >> probOnTrack = " << probOnTrack << " = " << probOnTrackWMulti << " * " << probOnTrackTerm;
  return probOnTrack;
}

void ProbQXYAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Producer that creates SiPixel Charge and shape probabilities combined for tracks");
  desc.addUntracked<bool>("debug", false);
  desc.add<std::string>("pixelCPE", {"PixelCPEClusterRepair"})->setComment("PixelCPE name");
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"))
      ->setComment("Input track collection for the producer");
  desc.add<edm::InputTag>("dedxInfos", edm::InputTag("dedxHitInfo"))->setComment("Input track dedxHitInfo");
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ProbQXYAna);
