#ifndef _PLOTSET_H_
#define _PLOTSET_H_
static const unsigned nPlot = 20U;
enum PlotName{
  FracSatVNstrips = 0, // fraction of saturated strips vs. number of strips
  DedXSig,             // dE/dX significance
  XYPar0,//d0
  XYPar1,//phi0
  XYPar2,//radius
  RZPar0,
  RZPar1,
  RZcurv,              // RZ curvature/uncertainty
  E55,
  F51,                 // frac 51
  HcalIso,             // Hcal Iso
  Dist,
  ABCD,
  Spike,
  EcalBarrel,
  EcalEndCup,
  EcalAll,
  PileUp_DedXSig,
  PileUp_f51
};
class PlotSet
{
public:
  PlotSet(){
	plots_.resize(nPlot);
	profile_.resize(nPlot);
  }
  ~PlotSet(){}
  void CreatPlot(const PlotName pn, TH1* h){ 
	plots_[pn] = h;
  }
  void CreatProfile(const PlotName pn, TProfile* p){ 
	profile_[pn] = p;
  }
  TH1 * GetPlot(const PlotName pn){ return plots_[pn]; }
  TProfile * GetProfile(const PlotName pn){ return profile_[pn]; }
  void WritePlot(){
	for(int pn=0;pn<nPlot;pn++){
	TH1 *h = plots_[pn];
	if(h){ 
	  h->Write();
	  }
	}
  }
  void WriteProfile(){
	for(int pn=0;pn<nPlot;pn++){
	TProfile *p = profile_[pn];
	if(p){ 
	  p->Write();
	  }
	}
  }
private:
  vector<TH1*> plots_;
  vector<TProfile*> profile_;
};
#endif

