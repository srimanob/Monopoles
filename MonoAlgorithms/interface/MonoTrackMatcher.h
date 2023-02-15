#ifndef Monopoles_MonoAlgorithms_MonoTrackMatcher_h
#define Monopoles_MonoAlgorithms_MonoTrackMatcher_h

////////////////////////////////////
// Use MonoTrackExtrapolator functions
// to match monopole tracks to calorimeter
// clusters.
////////////////////////////////////

#include "Monopoles/MonoAlgorithms/interface/MonoTrackExtrapolator.h"

#include "Monopoles/MonoAlgorithms/interface/MonoEcalCluster.h"
#include "Monopoles/MonoAlgorithms/interface/MonoTrack.h"
#include "Monopoles/MonoAlgorithms/interface/MonoEcalObs0.h"

namespace reco {
class CaloCluster;
}

namespace Mono {

class MonoTrackMatcher {

public:
  inline MonoTrackMatcher(double dRcut)
    :m_dRcut(dRcut)
  { }

  inline virtual ~MonoTrackMatcher() { }

  // match
  void match(unsigned nClusters, const MonoEcalCluster *clusters
    ,const EBmap &map
    ,unsigned nTracks, const MonoTrack *tracks
    ,std::vector<int> &matchMap, std::vector<double> &distances);

  void match(unsigned nClusters, const reco::CaloCluster **clusters
    ,unsigned nTracks, const MonoTrack *tracks
    ,std::vector<int> &matchMap, std::vector<double> &distances
    ,const bool isBarrel=true);


  class MatchInfo {
    public:
    inline MatchInfo():dist(999.),ic(0),it(0){}
    inline MatchInfo(double d,double c, double t):dist(d),ic(c),it(t){ }

    inline bool operator >(const MatchInfo &r) const {
      return dist > r.dist;
    }
   
    inline bool operator <(const MatchInfo &r) const {
      return dist < r.dist;
    }

    inline double getDist() const { return dist;}
    inline unsigned getic() const { return ic; }
    inline unsigned getit() const { return it; }

    private:
    double dist;
    unsigned ic;
    unsigned it;


  };

private:
  inline MonoTrackMatcher() { }

  static constexpr double s_ecalRad = 129.;
  static constexpr double s_EEz = 3.144;

  double m_dRcut;
  
};

} // end mono namespace

#endif
