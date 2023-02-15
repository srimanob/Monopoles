
#include "Monopoles/MonoAlgorithms/interface/MonoTrackMatcher.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <algorithm>

namespace Mono {

void Useless_asdf(){
  if(CLHEP::electron_charge==0) std::cout << "useless" << std::endl;
}

void MonoTrackMatcher::match(const unsigned nClusters, const MonoEcalCluster *clusters
  ,const EBmap &map
  ,const unsigned nTracks, const MonoTrack *tracks
  ,std::vector<int> & matchMap, std::vector<double> & distances)
{

  matchMap.clear();
  distances.clear();

  matchMap.resize(nTracks);
  distances.resize(nTracks);
  for ( unsigned c=0; c != nTracks; c++){
    matchMap[c] = -1;
    distances[c] = 999;
  }

  const unsigned nMatch = nClusters*nTracks;
  std::vector<MatchInfo> matchInfoMap(nMatch);

  MonoTrackExtrapolator extrap;

//  unsigned M=0;

  const double myInf = log(0);

  // cycle over all pairs of clusters and tracks
  // build map of distances
  for ( unsigned c=0; c != nClusters; c++ ) {
    const MonoEcalCluster & clust = clusters[c];
    const double cphi = map.phi( clust.iphi() );
    const double ceta = map.eta( clust.ieta() ); 

    for ( unsigned t=0; t != nTracks; t++ ) {
      const MonoTrack & track = tracks[t];
      const double tz = extrap.zVr(track.rzp0(),track.rzp1(),track.rzp2(),s_ecalRad);
      const double teta = extrap.eta(tz,s_ecalRad);
      const double tphi = extrap.phiVr(track.xyp0(),track.xyp1(),track.xyp2(),s_ecalRad);

      // check for NAN in tphi
      // if it is continue default distance is 999
      if ( tphi != tphi ) continue;
      if ( teta != teta ) continue;
      if ( teta == myInf || teta == -myInf ) continue;

      const double dR = reco::deltaR(ceta,cphi,teta,tphi);
      matchInfoMap[c*nTracks+t] = MatchInfo(dR,c,t);  

      assert(matchInfoMap[c*nTracks+t].getic()<nClusters);
      assert(matchInfoMap[c*nTracks+t].getit()<nTracks);
      assert(matchInfoMap[c*nTracks+t].getDist()<999.);

    }
  }


  std::vector<bool> taken(nClusters);
  for ( unsigned c=0; c != nClusters; c++ )
    taken[c] = false;

  unsigned nUsed = 0;
  std::sort(matchInfoMap.begin(),matchInfoMap.end());
  for ( unsigned i=0; i != nMatch && nUsed < nClusters && nUsed < nTracks; i++ ) {
    unsigned ic = matchInfoMap[i].getic();
    unsigned it = matchInfoMap[i].getit();
    double dist = matchInfoMap[i].getDist();

    if ( dist > m_dRcut ) break;

    if ( matchMap[it] < 0 && !taken[ic] ) {
      matchMap[it] = ic;
      taken[ic] = 1;
      distances[it] = dist;
      nUsed++;
    } 
    
  }
  
 

}


void MonoTrackMatcher::match(const unsigned nClusters, const reco::CaloCluster **clusters
  ,const unsigned nTracks, const MonoTrack *tracks
  ,std::vector<int> & matchMap, std::vector<double> & distances, const bool isBarrel)
{

  matchMap.clear();
  distances.clear();

  matchMap.resize(nTracks);
  distances.resize(nTracks);
  for ( unsigned c=0; c != nTracks; c++){
    matchMap[c] = -1;
    distances[c] = 999;
  }

  const unsigned nMatch = nClusters*nTracks;
  std::vector<MatchInfo> matchInfoMap(nMatch);

  MonoTrackExtrapolator extrap;

  //unsigned M=0;

  const double myInf = log(0);

  // cycle over all pairs of clusters and tracks
  // build map of distances
  for ( unsigned c=0; c != nClusters; c++ ) {
    const reco::CaloCluster * clust = clusters[c];
    const double cphi = clust->phi();
    const double ceta = clust->eta();

    for ( unsigned t=0; t != nTracks; t++ ) {
      const MonoTrack & track = tracks[t];
      const double tz = isBarrel ? extrap.zVr(track.rzp0(),track.rzp1(),track.rzp2(),s_ecalRad) : extrap.rVz(track.rzp0(),track.rzp1(),track.rzp2(),s_EEz);
      const double teta = isBarrel ? extrap.eta(tz,s_ecalRad) : extrap.eta(s_EEz,tz);
      const double tphi = isBarrel ? extrap.phiVr(track.xyp0(),track.xyp1(),track.xyp2(),s_ecalRad) : extrap.phiVr(track.xyp0(),track.xyp1(),track.xyp2(),tz);

      // check for NAN in tphi
      // if it is continue default distance is 999
      if ( tphi != tphi ) continue;
      if ( teta != teta ) continue;
      if ( teta == myInf || teta == -myInf ) continue;

      const double dR = reco::deltaR(cphi,ceta,teta,tphi);
      matchInfoMap[c*nTracks+t] = MatchInfo(dR,c,t);  

      assert(matchInfoMap[c*nTracks+t].getic()<nClusters);
      assert(matchInfoMap[c*nTracks+t].getit()<nTracks);
      assert(matchInfoMap[c*nTracks+t].getDist()<999.);

    }
  }


  std::vector<bool> taken(nClusters);
  for ( unsigned c=0; c != nClusters; c++ )
    taken[c] = false;

  unsigned nUsed = 0;
  std::sort(matchInfoMap.begin(),matchInfoMap.end());
  for ( unsigned i=0; i != nMatch && nUsed < nClusters && nUsed < nTracks; i++ ) {
    unsigned ic = matchInfoMap[i].getic();
    unsigned it = matchInfoMap[i].getit();
    double dist = matchInfoMap[i].getDist();

    if ( dist > m_dRcut ) break;

    if ( matchMap[it] < 0 && !taken[ic] ) {
      matchMap[it] = ic;
      taken[ic] = 1;
      distances[it] = dist;
      nUsed++;
    } 
    
  }
  
 

}

} // end mono namespace
