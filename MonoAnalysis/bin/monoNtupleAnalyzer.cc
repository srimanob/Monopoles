//////////////////////////////////////////////
// monoNtupleAnalyzer.cc  C Cowden
// 11/5/2015
//---------------------------------------------
// This program is intended to analyze the ntuple 
// produced by MonoNtupleDumper.
// Make the following plots for several cases:
// + Ecal
// |`+ energy in 5x5 cluster
// |`+ frac51
// + Track
// |`+ dE/dX significance
//


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unistd.h>

#include "TH1.h"
#include "TH1D.h"
#include "TTree.h"
#include "TH2D.h"
#include "TFile.h"

using namespace std;

//================================================
//  helper classes and functions

// enumeration to help number plots
enum PlotNum {
  FracSatVNstrips = 0, // fraction of saturated strips vs. number of strips
  DedXSig,             // dE/dX significance
  RZcurv,              // RZ curvature/uncertainty
  MCtag,               // MC tag to Ecal cluster
  F51,                 // frac 51
  F15,                 // frac 15 (orthogonal to analysis)
  Energy,              // energy in 5x5 cluster
  HcalIso,             // Hcal Iso
  ABCD
};

// number of plots in PlotNum Enum
static const unsigned gNplots = 9U;


// PlotSet - a class to help hold the several plots to keep
// at various stages of selections
class PlotSet {

public:
  PlotSet() { plots_.resize(gNplots); }
  PlotSet( const PlotSet &ps ): plots_(ps.plots_) { }
 
  ~PlotSet() { }
  
  void AddPlot(const PlotNum pn, TH1 *h) { plots_[pn] = h; }
  TH1 * GetPlot(const PlotNum pn) { return plots_[pn]; }

  void WritePlots() {
    for ( unsigned p=0; p != gNplots; p++ ) {
      TH1 * h = plots_[p];
      if ( h ) h->Write();
    }
  }
  
private:
  vector<TH1*> plots_;
};

// Candidate
// The candidate class is used to create a sortable list of monopole candidates.
class MonoCandidate {
public:
  MonoCandidate() { }
  
  MonoCandidate(double sh, double satsh, double dedxsig, double tiso, 
		double xyp0, double xyp1, double xyp2,
		double rzp0, double rzp1, double rzp2, 
		double dist, double f51, double f15, double e55, double met, double hiso, 
		double eta, double phi ):
    subHits_(sh),subSatHits_(satsh),dEdXSig_(dedxsig),tIso_(tiso),
    xyp0_(xyp0),xyp1_(xyp1),xyp2_(xyp2),
    rzp0_(rzp0),rzp1_(rzp1),rzp2_(rzp2),
    dist_(dist),f51_(f51),f15_(f15),e55_(e55),met_(met),hIso_(hiso),
    eta_(eta),phi_(phi) { }

  MonoCandidate(const MonoCandidate & mc) :
    subHits_(mc.subHits_),subSatHits_(mc.subSatHits_),dEdXSig_(mc.dEdXSig_),tIso_(mc.tIso_),
    xyp0_(mc.xyp0_),xyp1_(mc.xyp1_),xyp2_(mc.xyp2_),
    rzp0_(mc.rzp0_),rzp1_(mc.rzp1_),rzp2_(mc.rzp2_),
    dist_(mc.dist_),f51_(mc.f51_),f15_(mc.f15_),e55_(mc.e55_),met_(mc.met_),hIso_(mc.hIso_),
    eta_(mc.eta_),phi_(mc.phi_) { }
  
  ~MonoCandidate() { }
  
  // comparitor operator
  bool operator <(const MonoCandidate &mc) const {
    if ( dEdXSig_ > mc.dEdXSig_ ) return true;
    else if ( dEdXSig_ == mc.dEdXSig_ ) {
      if ( f51_ > mc.f51_ ) return true;
      else return false;
    } 
    return false;
  }
  
  double subHits_;
  double subSatHits_;
  double dEdXSig_;
  double tIso_;
  double xyp0_;
  double xyp1_;
  double xyp2_;
  double rzp0_;
  double rzp1_;
  double rzp2_;  
  double dist_;
  double f51_;
  double f15_;
  double e55_;
  double met_;
  double hIso_;
  double eta_;
  double phi_;
};

// MonoCutFlow - a class to help hold several PlotSets and tables of 
// events passing cuts and triggers.
// This is expected to exist for each trigger of interest.  
// The write function will add a directory to output root file and
// dump the latex tables.
class MonoCutFlow {

public:
  // This constructor will create the necessary plots (it assumes an output root file
  // exists and there is no problem with directory structures.
  MonoCutFlow(const string & trName, TFile *oFile) : trigName_(trName) { 
    
    // create a new directory in the output tree "tree"
    oFile->mkdir(trigName_.c_str());
    
    // create table for nm1Count (this has the same number scheme as nm1Plots_)
    nm1Count_.resize(nCuts_+1U,0U);
    
    //  create n-1 plots
    char name[50];
    sprintf(name,"%d",series_);
    string nm1base = "nm1_" + string(name) + "_";
    nm1Plots_.resize(nCuts_+1U); 
    
    string trgnm1name = nm1base + "HLT_";
    PlotSet & z = nm1Plots_[0];
    
    z.AddPlot(FracSatVNstrips,new TH2D((trgnm1name+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
    z.AddPlot(DedXSig,new TH1D((trgnm1name+"dEdXSig").c_str(),"",100,0,30));
    z.AddPlot(RZcurv,new TH1D((trgnm1name+"RZcurv").c_str(),"",100,-10,10));
    z.AddPlot(MCtag,new TH1D((trgnm1name+"MCtag").c_str(),"",2,-0.5,1.5));
    z.AddPlot(F51,new TH1D((trgnm1name+"F51").c_str(),"",100,0,1.1));
    z.AddPlot(F15,new TH1D((trgnm1name+"F15").c_str(),"",100,0,1.1));
    z.AddPlot(Energy,new TH1D((trgnm1name+"E55").c_str(),"",100,0,1500));
    z.AddPlot(HcalIso,new TH1D((trgnm1name+"HcalIso").c_str(),"",100,0,30));
    z.AddPlot(ABCD,new TH2D((trgnm1name+"ABCD").c_str(),"",100,0,1.1,100,0,30));
    
    for ( unsigned c=0; c != nCuts_; c++ ) {
      // create plots for a particular cut
      string cutName( cutNames_[c] );
      PlotSet & z = nm1Plots_[c+1U];
      
      string cutnm1name = nm1base + cutName + "_";
      
      z.AddPlot(FracSatVNstrips,new TH2D((cutnm1name+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
      z.AddPlot(DedXSig,new TH1D((cutnm1name+"dEdXSig").c_str(),"",100,0,30));
      z.AddPlot(RZcurv,new TH1D((cutnm1name+"RZcurv").c_str(),"",100,-10,10));
      z.AddPlot(MCtag,new TH1D((cutnm1name+"MCtag").c_str(),"",2,-0.5,1.5));
      z.AddPlot(F51,new TH1D((cutnm1name+"F51").c_str(),"",100,0,1.1));
      z.AddPlot(F15,new TH1D((cutnm1name+"F15").c_str(),"",100,0,1.1));
      z.AddPlot(Energy,new TH1D((cutnm1name+"E55").c_str(),"",100,0,1500));
      z.AddPlot(HcalIso,new TH1D((cutnm1name+"HcalIso").c_str(),"",100,0,30));
      z.AddPlot(ABCD,new TH2D((cutnm1name+"ABCD").c_str(),"",100,0,1.1,100,0,30));
    }

    // create a table to cut flow counts (this has the same numbering as cutFlow_)
    cutFlowCount_.resize(nCuts_+1U,0U);  // the four cuts + trigger
    
    // create cut flow plots
    cutFlow_.resize(nCuts_+1U);  // the four cuts + trigger 
    string cutFlowBase = "flow_" + string(name) + "_";
    string noCut = "HLT";
    string cutFlowName = cutFlowBase + noCut + "_";
    
    PlotSet & y = cutFlow_[0];
    y.AddPlot(FracSatVNstrips,new TH2D((cutFlowName+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
    y.AddPlot(DedXSig,new TH1D((cutFlowName+"dEdXSig").c_str(),"",100,0,30));
    y.AddPlot(RZcurv,new TH1D((cutFlowName+"RZcurv").c_str(),"",100,-10,10));
    y.AddPlot(MCtag,new TH1D((cutFlowName+"MCtag").c_str(),"",2,-0.5,1.5));
    y.AddPlot(F51,new TH1D((cutFlowName+"F51").c_str(),"",100,0,1.1));
    y.AddPlot(F15,new TH1D((cutFlowName+"F15").c_str(),"",100,0,1.1));
    y.AddPlot(Energy,new TH1D((cutFlowName+"E55").c_str(),"",100,0,1500));
    y.AddPlot(HcalIso,new TH1D((cutFlowName+"HcalIso").c_str(),"",100,0,30));
    y.AddPlot(ABCD,new TH2D((cutFlowName+"ABCD").c_str(),"",100,0,1.1,100,0,30));
    
    for ( unsigned c=0; c != nCuts_; c++ ) {
      // create plots for a particular cut
      string cutName( cutNames_[c] );
      PlotSet &z = cutFlow_[c+1U];
      
      cutFlowName = cutFlowBase + cutName + "_";
      
      z.AddPlot(FracSatVNstrips,new TH2D((cutFlowName+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
      z.AddPlot(DedXSig,new TH1D((cutFlowName+"dEdXSig").c_str(),"",100,0,30));
      z.AddPlot(RZcurv,new TH1D((cutFlowName+"RZcurv").c_str(),"",100,-10,10));
      z.AddPlot(MCtag,new TH1D((cutFlowName+"MCtag").c_str(),"",2,-0.5,1.5));
      z.AddPlot(F51,new TH1D((cutFlowName+"F51").c_str(),"",100,0,1.1));
      z.AddPlot(F15,new TH1D((cutFlowName+"F15").c_str(),"",100,0,1.1));
      z.AddPlot(Energy,new TH1D((cutFlowName+"E55").c_str(),"",100,0,1500));
      z.AddPlot(HcalIso,new TH1D((cutFlowName+"HcalIso").c_str(),"",100,0,30));
      z.AddPlot(ABCD,new TH2D((cutFlowName+"ABCD").c_str(),"",100,0,1.1,100,0,30));

    }

    // increment the series
    series_++;

  }

  MonoCutFlow( const MonoCutFlow & mf ): trigName_(mf.trigName_), nm1Plots_(mf.nm1Plots_),
  cutFlow_(mf.cutFlow_), nm1Count_(mf.nm1Count_), cutFlowCount_(mf.cutFlowCount_) { }

  ~MonoCutFlow() { }

  // perform the analysis
  void doAnalysis(bool trg, const vector<MonoCandidate> & cands) {

    // apply all cuts
    if ( cands.size() > 0 ) {
      const MonoCandidate & cand = cands[0];
  
      // get the results of the cuts (or sets of cuts)
      const bool qualRes = evalQuality(cand);
      const bool eRes = evalE(cand);
      const bool f51Res = evalF51(cand);
      const bool dEdXRes = evaldEdX(cand);

      // apply all but trigger
      if ( qualRes && eRes && f51Res && dEdXRes ) {

	// count the event
	nm1Count_[0]++;
	
	// get the plot HLT is 0
	PlotSet & z = nm1Plots_[0];
	z.GetPlot(FracSatVNstrips)->Fill(cand.subHits_,cand.subSatHits_);
	z.GetPlot(DedXSig)->Fill(cand.dEdXSig_);
	z.GetPlot(RZcurv)->Fill(cand.rzp2_);
	// z.GetPlot(MCtag) needs to be filled
	z.GetPlot(F51)->Fill(cand.f51_);
	z.GetPlot(F15)->Fill(cand.f15_);
	z.GetPlot(Energy)->Fill(cand.e55_);
	z.GetPlot(HcalIso)->Fill(cand.hIso_);
	z.GetPlot(ABCD)->Fill(cand.f51_,cand.dEdXSig_);

      } 

      // apply all but quality cuts
      if ( eRes && f51Res && dEdXRes && trg ) {

	// count the event
	nm1Count_[1]++;

	// get the plot quality is 1
	PlotSet & z = nm1Plots_[1];
	z.GetPlot(FracSatVNstrips)->Fill(cand.subHits_,cand.subSatHits_);
	z.GetPlot(DedXSig)->Fill(cand.dEdXSig_);
	z.GetPlot(RZcurv)->Fill(cand.rzp2_);
	// z.GetPlot(MCtag) needs to be filled
	z.GetPlot(F51)->Fill(cand.f51_);
	z.GetPlot(F15)->Fill(cand.f15_);
	z.GetPlot(Energy)->Fill(cand.e55_);
	z.GetPlot(HcalIso)->Fill(cand.hIso_);
	z.GetPlot(ABCD)->Fill(cand.f51_,cand.dEdXSig_);

      }

      // apply all but energy cut
      if ( qualRes && f51Res && dEdXRes && trg ) {

	// count the event
	nm1Count_[2]++;

	// get the plot energy is 2
	PlotSet & z = nm1Plots_[2];
	z.GetPlot(FracSatVNstrips)->Fill(cand.subHits_,cand.subSatHits_);
	z.GetPlot(DedXSig)->Fill(cand.dEdXSig_);
	z.GetPlot(RZcurv)->Fill(cand.rzp2_);
	// z.GetPlot(MCtag) needs to be filled
	z.GetPlot(F51)->Fill(cand.f51_);
	z.GetPlot(F15)->Fill(cand.f15_);
	z.GetPlot(Energy)->Fill(cand.e55_);
	z.GetPlot(HcalIso)->Fill(cand.hIso_);
	z.GetPlot(ABCD)->Fill(cand.f51_,cand.dEdXSig_);

      }

      // apply all but f51 cut
      if ( qualRes && eRes && dEdXRes && trg ) {

	// count the event
	nm1Count_[3]++;

	// get the plot f51 is 3
	PlotSet & z = nm1Plots_[3];
	z.GetPlot(FracSatVNstrips)->Fill(cand.subHits_,cand.subSatHits_);
	z.GetPlot(DedXSig)->Fill(cand.dEdXSig_);
	z.GetPlot(RZcurv)->Fill(cand.rzp2_);
	// z.GetPlot(MCtag) needs to be filled
	z.GetPlot(F51)->Fill(cand.f51_);
	z.GetPlot(F15)->Fill(cand.f15_);
	z.GetPlot(Energy)->Fill(cand.e55_);
	z.GetPlot(HcalIso)->Fill(cand.hIso_);
	z.GetPlot(ABCD)->Fill(cand.f51_,cand.dEdXSig_);

      }

      // apply all but dedx cut
      if ( qualRes && eRes && f51Res && trg ) {

	// count the event
	nm1Count_[4]++;

	// get the plot dedx is 4
	PlotSet & z = nm1Plots_[4];
	z.GetPlot(FracSatVNstrips)->Fill(cand.subHits_,cand.subSatHits_);
	z.GetPlot(DedXSig)->Fill(cand.dEdXSig_);
	z.GetPlot(RZcurv)->Fill(cand.rzp2_);
	// z.GetPlot(MCtag) needs to be filled
	z.GetPlot(F51)->Fill(cand.f51_);
	z.GetPlot(F15)->Fill(cand.f15_);
	z.GetPlot(Energy)->Fill(cand.e55_);
	z.GetPlot(HcalIso)->Fill(cand.hIso_);
	z.GetPlot(ABCD)->Fill(cand.f51_,cand.dEdXSig_);

      }
     
      // --- make plots along the cut flow ( in the order listed in cutNames_ )
      // **** MAKE SURE cutNames_ IS THE SAME ORDERING AS accepts!!!! ****
      // first is the trigger case
      bool accept = trg;

      if ( accept ) {

      	// count event in cut flow
      	cutFlowCount_[0]++;

  	// fill the HLT plots
  	PlotSet & y = cutFlow_[0];
	y.GetPlot(FracSatVNstrips)->Fill(cand.subHits_,cand.subSatHits_);
	y.GetPlot(DedXSig)->Fill(cand.dEdXSig_);
	y.GetPlot(RZcurv)->Fill(cand.rzp2_);
	// y.GetPlot(MCtag) needs to be filled
	y.GetPlot(F51)->Fill(cand.f51_);
	y.GetPlot(F15)->Fill(cand.f15_);
	y.GetPlot(Energy)->Fill(cand.e55_);
	y.GetPlot(HcalIso)->Fill(cand.hIso_);
	y.GetPlot(ABCD)->Fill(cand.f51_,cand.dEdXSig_);

      }

      bool accepts[nCuts_] = { qualRes, eRes, f51Res, dEdXRes };
      
      for ( unsigned c=0; c != nCuts_; c++ ) {
  	accept = accept && accepts[c];
	if ( accept )  {

      	  // count event in cut flow
      	  cutFlowCount_[c+1U]++;

  	  // fill the HLT plots
  	  PlotSet & y = cutFlow_[c+1U];
	  y.GetPlot(FracSatVNstrips)->Fill(cand.subHits_,cand.subSatHits_);
	  y.GetPlot(DedXSig)->Fill(cand.dEdXSig_);
	  y.GetPlot(RZcurv)->Fill(cand.rzp2_);
	  // y.GetPlot(MCtag) needs to be filled
	  y.GetPlot(F51)->Fill(cand.f51_);
	  y.GetPlot(F15)->Fill(cand.f15_);
	  y.GetPlot(Energy)->Fill(cand.e55_);
	  y.GetPlot(HcalIso)->Fill(cand.hIso_);
	  y.GetPlot(ABCD)->Fill(cand.f51_,cand.dEdXSig_);

	}
      }
    }

  }


  // write histograms to the root file
  void WritePlots(TFile *oFile ) {
    oFile->cd(trigName_.c_str());
  
    // write N-1 plots  
    nm1Plots_[0].WritePlots();
    for ( unsigned c=0; c != nCuts_; c++ ) nm1Plots_[c+1U].WritePlots();

    // write cut flow plots
    cutFlow_[0].WritePlots();
    for ( unsigned c=0; c != nCuts_; c++ ) cutFlow_[c+1U].WritePlots();

  }

  // Dump tabulated event counts for cut flows
  void DumpTables() {
    // Dump tables with a simple list format
    // TrgName N-1 CutMissing NumEvents
    // ...
    // TrgName CutFlow Cut NumEvents
    cout << trigName_ << " N-1 TRG " << nm1Count_[0] << endl;
    for ( unsigned c=0; c != nCuts_; c++ ) {
      cout << trigName_ << " N-1 " << cutNames_[c] << " " << nm1Count_[c+1U] << endl;
    }
    cout << trigName_ << " CutFlow TRG " << cutFlowCount_[0] << endl;
    for ( unsigned c=0; c != nCuts_; c++ ) {
      cout << trigName_ << " CutFlow " << cutNames_[c] << " " << cutFlowCount_[c+1U] << endl;
    }

  }

private:
  MonoCutFlow() { }

  // keep the trigger name
  string trigName_;

  // keep a set of plots for each n-1 and step in the cut flow
  vector<PlotSet>  nm1Plots_; // n-1 plots
  vector<PlotSet>  cutFlow_; // plot set for progressive cuts

  // keep information to dump tables of number of events passing cuts (N-1) and cut flow
  vector<unsigned> nm1Count_;
  vector<unsigned> cutFlowCount_;

  // private inline members to evaluate sets of cuts
  inline bool evalQuality(const MonoCandidate & c) { return (c.dist_ < dCut_ && c.hIso_ < hIsoCut_ 
							     && fabs(c.xyp0_) < xyp0Cut_ && fabs(c.xyp2_) > xyp2Cut_ 
							     && fabs(c.rzp0_) < rzp0Cut_ && fabs(c.rzp1_) < rzp1Cut_ && fabs(c.rzp2_) < rzp2Cut_); }
  inline bool evalE(const MonoCandidate & c ) { return c.e55_ > eCut_; }
  inline bool evalF51(const MonoCandidate & c) { return c.f51_ > f51Cut_ ; }
  inline bool evaldEdX(const MonoCandidate & c) { return c.dEdXSig_ > dEdXSigCut_; }

  // add a number to make plot names unique
  static unsigned series_;

  //  number of cuts
  static const unsigned nCuts_ = 6U; // quality cuts + 5 analysis cuts
  
  static const char *cutNames_[nCuts_];
  
  // cuts
  // track/quality cuts will be counted as one
  static const double dCut_; 
  static const double hIsoCut_;
  static const double xyp0Cut_;
  static const double xyp2Cut_;
  static const double rzp0Cut_;
  static const double rzp1Cut_;
  static const double rzp2Cut_;

  // analysis cuts (5)
  static const double eCut_;
  static const double metCut_;
  static const double f51Cut_;
  static const double f15Cut_;
  static const double dEdXSigCut_;

  // The loose cuts will not be counted as special cuts
  // these values are only used in the ABCD plot boundary/region
  // definitions.
  static const double f51LooseCut_;
  static const double dEdXSigLooseCut_;
};

unsigned MonoCutFlow::series_ = 0U;

const double MonoCutFlow::dCut_ = 0.5;
const double MonoCutFlow::hIsoCut_ = 10.;
const double MonoCutFlow::xyp0Cut_ = 0.6;
const double MonoCutFlow::xyp2Cut_ = 1000;
const double MonoCutFlow::rzp0Cut_ =  10;
const double MonoCutFlow::rzp1Cut_ = 999;
const double MonoCutFlow::rzp2Cut_ = 0.005;
const double MonoCutFlow::eCut_ = 200.;
const double MonoCutFlow::metCut_ = 0.;
const double MonoCutFlow::f51Cut_ = 0.95;
const double MonoCutFlow::f15Cut_ = 0.;
const double MonoCutFlow::dEdXSigCut_ = 9;
const double MonoCutFlow::f51LooseCut_ = 0.8;
const double MonoCutFlow::dEdXSigLooseCut_ = 7;

const char *MonoCutFlow::cutNames_[MonoCutFlow::nCuts_] = {
 "QualityCuts",
 "ECut",
 "METCut",
 "F51Cut",
 "F15Cut",
 "dEdXSigCut"};

//
//  print usage
void printUsage() {
  cout << "monoNtupleAnalyzer -i ntuple.root -o outFile.root" << endl;
}


int main ( int argc, char **argv ) {

  // trigger path list
  vector<string> trigNameList;

  trigNameList.push_back("HLT_Photon175_v");
  trigNameList.push_back("HLT_Photon200_v");
  trigNameList.push_back("HLT_PFMET300_v");
  trigNameList.push_back("HLT_MET200_v");
  trigNameList.push_back("HLT_PFMET250_HBHECleaned_v");
  trigNameList.push_back("HLT_CaloMET350_HBHECleaned_v");

  string inFileName;
  string outFileName;
  
  // parse command line options
  //
  int c;
  opterr = 0;
  while (( c=getopt(argc,argv,"hi:o:")) != -1 )
    switch (c) {
      case 'h':
	printUsage();
	return 0;
      case 'i':
	inFileName = optarg;
      	break;
      case 'o':
	outFileName = optarg;
	break;
    }

  if ( inFileName == "" ) {
    cerr << "You must provide an ntuple file to analyze!" << endl;
    return 1;
  }

  const unsigned nTriggerNames = trigNameList.size();


  // trigger branches
  vector<bool> * trigResults = 0;
  vector<string> * trigNames = 0; 

  unsigned nCandidates;
  vector<double> * subHits = 0;
  vector<double> * subSatHits = 0;
  vector<double> * dEdXSig = 0;
  vector<double> * tIso = 0;
  vector<double> * xyp0 = 0;
  vector<double> * xyp1 = 0;
  vector<double> * xyp2 = 0;
  vector<double> * rzp0 = 0;
  vector<double> * rzp1 = 0;
  vector<double> * rzp2 = 0;

  vector<double> * dist = 0;
  vector<double> * f51 = 0;
  vector<double> * f15 = 0;  
  vector<double> * e55 = 0;
  double met = 0;
  vector<double> * hIso = 0;
  vector<double> * eta = 0;
  vector<double> * phi = 0;
  
  // open the ntuple
  TFile inFile(inFileName.c_str());
  TTree *tree = (TTree*)inFile.Get("monopoles");

  tree->SetBranchAddress("trigResult",&trigResults); 
  tree->SetBranchAddress("trigNames",&trigNames);

  tree->SetBranchAddress("cand_N",&nCandidates);
  tree->SetBranchAddress("cand_dist",&dist);
  tree->SetBranchAddress("cand_SubHits",&subHits);
  tree->SetBranchAddress("cand_SatSubHits",&subSatHits);
  tree->SetBranchAddress("cand_dEdXSig",&dEdXSig);
  tree->SetBranchAddress("cand_TIso",&tIso);
  tree->SetBranchAddress("cand_f51",&f51);
  tree->SetBranchAddress("cand_f15",&f15);
  tree->SetBranchAddress("cand_e55",&e55);
  //tree->SetBranchAddress("mpt_pt",met);
  tree->SetBranchAddress("cand_HIso",&hIso);
  tree->SetBranchAddress("cand_XYPar0",&xyp0);
  tree->SetBranchAddress("cand_XYPar1",&xyp1);
  tree->SetBranchAddress("cand_XYPar2",&xyp2);
  tree->SetBranchAddress("cand_RZPar0",&rzp0);
  tree->SetBranchAddress("cand_RZPar1",&rzp1);
  tree->SetBranchAddress("cand_RZPar2",&rzp2);
  tree->SetBranchAddress("cand_eta",&eta);
  tree->SetBranchAddress("cand_phi",&phi);

  tree->SetBranchStatus("*",1); 

  // prepare the output file
  TFile *oFile = new TFile(outFileName.c_str(),"RECREATE");

  // create a MonoCutFlow for the case of no trigger
  MonoCutFlow noTriggerAnalysis("NoTRG",oFile);  

  // create a cut flow analysis object for each
  // interesting trigger path.
  vector<MonoCutFlow> triggerAnalyses;
  for ( unsigned i=0; i != nTriggerNames; i++ ) triggerAnalyses.push_back( MonoCutFlow(trigNameList[i],oFile) );
  
  vector<MonoCandidate> candVec(10); 

  // begin the event loop
  const unsigned nEvents = tree->GetEntries();
  for ( unsigned ev=0; ev != nEvents; ev++ ) {
    tree->GetEntry(ev);
    
    //
    met=100000.; //temporary met

    // build the candidate
    if ( nCandidates > candVec.size() ) candVec.resize(nCandidates);
    for ( unsigned i=0; i != nCandidates; i++ ) {
      candVec[i] = MonoCandidate(
				 (*subHits)[i],
				 (*subSatHits)[i],
				 (*dEdXSig)[i],
				 (*tIso)[i],
				 (*xyp0)[i],
				 (*xyp1)[i],
				 (*xyp2)[i],
				 (*rzp0)[i],
				 (*rzp1)[i],
				 (*rzp2)[i],
				 (*dist)[i],
				 (*f51)[i],
				 (*f15)[i], 
				 (*e55)[i],
				 met,
				 (*hIso)[i],
				 (*eta)[i],
				 (*phi)[i]
				 );
    }

    // sort the monopole candidates
    sort(candVec.begin(),candVec.begin()+nCandidates);
    
    // do the analysis for no trigger selection
    noTriggerAnalysis.doAnalysis(true,candVec);
    
    // loop over triggers check if name is in the list
    const unsigned nTrigs = trigNames->size();
    for ( unsigned tr=0; tr != nTrigs; tr++ ) {
      const string & trName = (*trigNames)[tr];
      for ( unsigned tn=0; tn != nTriggerNames; tn++ ) {
	if ( trName.find(trigNameList[tn].c_str()) != string::npos ) {
	  // a very verbose debug statement to print out (just in case)
	  //cout << "trName " << trName << " " << trName.find(trigNameList[tn].c_str()) << " with result " << (*trigResults)[tr] << endl;	
	  triggerAnalyses[tn].doAnalysis((*trigResults)[tr],candVec);

	}
      }
    } 
  }

  // write the plots to the output file
  noTriggerAnalysis.WritePlots(oFile);
  for ( unsigned i=0; i != nTriggerNames; i++ ) {
    triggerAnalyses[i].WritePlots(oFile);
  }
  
  // dump summary table of cut results for each trigger
  // TrgName N-1 CutMissing nEvents
  // ...
  // TrgName cut nEvents
  cout << "=====================================================" << endl;
  cout << "===== Dumping Cut Tables" << endl;
  noTriggerAnalysis.DumpTables();
  cout << "-----------------------------------------------------" << endl;
  for ( unsigned i=0; i != nTriggerNames; i++ ) {
    triggerAnalyses[i].DumpTables();
    cout << "-----------------------------------------------------" << endl;
  }
  cout << "==== End Cut Tables" << endl;
  cout << "=====================================================" << endl;
  
  oFile->Close();
  
  return 0;
}

