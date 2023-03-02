#ifndef _CANDIDATE_DATA_H_
#define _CANDIDATE_DATA_H_
class MonoCandidate
{
public:
  MonoCandidate(){}
  //This will be used in main function and to absord the data in root file
    MonoCandidate(double sh, double satsh,double dist,double hiso, double xyp0, double xyp1, double xyp2,
    double rzp0, double rzp1, double rzp2,
    double e55, double f51, double dedxsig,double eta, double event):
  subHits_(sh),subSatHits_(satsh),dist_(dist),hIso_(hiso),xyp0_(xyp0),xyp1_(xyp1),xyp2_(xyp2),
  rzp0_(rzp0),rzp1_(rzp1),rzp2_(rzp2),e55_(e55),f51_(f51),dEdXSig_(dedxsig),eta_(eta),
  event_(event) { }
  //This will be used in comparing with cut
  MonoCandidate(const MonoCandidate &mc) : 
    subHits_(mc.subHits_),subSatHits_(mc.subSatHits_),dist_(mc.dist_),hIso_(mc.hIso_),
    xyp0_(mc.xyp0_),xyp1_(mc.xyp1_),xyp2_(mc.xyp2_),
    rzp0_(mc.rzp0_),rzp1_(mc.rzp1_),rzp2_(mc.rzp2_),
    e55_(mc.e55_),f51_(mc.f51_),dEdXSig_(mc.dEdXSig_),eta_(mc.eta_),
    event_(mc.event_) { } 
        
  ~MonoCandidate() {}
  bool operator<(const MonoCandidate &mc)const{
   if(dEdXSig_>mc.dEdXSig_) return true;
   else if(dEdXSig_==mc.dEdXSig_){
        if(f51_>mc.f51_) return true;
        else return false;
        }
    else return false;
  }
  //All candidates variable
  double subHits_;
  double subSatHits_;
  double dEdXSig_;
  double xyp0_;
  double xyp1_;
  double xyp2_;
  double rzp0_;
  double rzp1_;
  double rzp2_;

  double dist_;
  double f51_;
  double e55_;
  double hIso_;
  double eta_;
  double event_;

};
#endif	
