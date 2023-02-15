#ifndef MONOANALYS_ECALMAPPER_H
#define MONOANALYS_ECALMAPPER_H


///////////////////////////////////////////////////////////////////
// C S Cowden 					15 February 2012  
// Helper class to fill array/histogram of Ecal energy/timing 
// (or any other information) mapped to eta phi locations.
///////////////////////////////////////////////////////////////////


#include <vector>
#include <cstring>



// forward class declarations
class TFileDirectory;
namespace edm { 
  class Event; 
  class EventSetup; 
}
class CaloGeometry;



namespace Mono {


struct marker {
  double eta;
  double phi;
  char text[100];

  marker(const double e, const double p, const char * t):
    eta(e),phi(p) { strncpy(text,t,100); }

};


class EcalMapper
{

  public:
  explicit EcalMapper(const edm::EventSetup &es );
  ~EcalMapper();

 
  // fill the maps.  If dir is passed, place maps in it
  // otherwise use default TFileService directory. 
  void fillMap(const edm::Event &ev, TFileDirectory *dir = 0);  

  // set TLatex marker (eta, phi, text);
  void setMarker(const double eta, const double phi, const char *);


  private:

  // calorimetry geometry
  const CaloGeometry* m_caloGeo;

  std::vector<marker> m_markers;


}; // end Ecal Mapper class



} // end Mono namespace 

#endif
