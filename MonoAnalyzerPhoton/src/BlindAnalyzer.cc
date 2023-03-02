//
//	BlindAnalyzer.cc
//	Created by  Shih Lin
//	
//	Just Draw some plots, not for do ABCD
//
#include <algorithm>
#include <string>
#include "../interface/Candidate_Data.h"
#include "../interface/PlotSet.h"
#include "../interface/MonoCuts_data.h"
using namespace std;

  void MonoCuts::doAnalysis_data(vector<MonoCandidate> &cand, unsigned nCandidates, bool TRG, unsigned ev)
  {
	Clear();

     for(unsigned c=0;c<nCandidates;c++){

	MonoCandidate &cands = cand[c];
	bool PreselectionCut = evalPreselection(cands);
	if(TRG && PreselectionCut ) Preselection.push_back(cands); 
	if(TRG ) TriggerSelection.push_back(cands); 
	
      }//for cand loop
    	
	//cut flow events calculating
        sort(Preselection.begin(),Preselection.begin()+Preselection.size());
        sort(TriggerSelection.begin(),TriggerSelection.begin()+TriggerSelection.size());
	if(Preselection.size()>0) 
	{
		PlotSet &z = CutFlow[0];
		PlotSet &x = Profile[0];
		     /*-----------------------
 			slice1 = 0 ~ 0.17
 			slice2 = 0.17 ~ 0.34
 			slice3 = 0.34 ~ 0.51
 			slice4 = 0.51 ~ 0.68
			slice5 = 0.68 ~ 0.85
			slice5 = 0.85 ~ 1.02
			---------------------*/
		Qual_count++;
		for(int i=0; i<Preselection.size();i++){
		//--- dEdxSig with slice f51 ---//
//		  if(Preselection[i].f51_ >= 0.00 && Preselection[i].f51_ < 0.17)	     z.GetPlot(F51_slice1)->Fill(Preselection[i].dEdXSig_);
//		  if(Preselection[i].f51_ >= 0.17 && Preselection[i].f51_ < 0.34)	     z.GetPlot(F51_slice2)->Fill(Preselection[i].dEdXSig_);
// 		  if(Preselection[i].f51_ >= 0.34 && Preselection[i].f51_ < 0.51)            z.GetPlot(F51_slice3)->Fill(Preselection[i].dEdXSig_);
//		  if(Preselection[i].f51_ >= 0.51 && Preselection[i].f51_ < 0.68)	     z.GetPlot(F51_slice4)->Fill(Preselection[i].dEdXSig_);
//		  if(Preselection[i].f51_ >= 0.68 && Preselection[i].f51_ < 0.85)	     z.GetPlot(F51_slice5)->Fill(Preselection[i].dEdXSig_);
//		  if(Preselection[i].f51_ >= 0.85 && Preselection[i].f51_ < 1.02)	     z.GetPlot(F51_slice6)->Fill(Preselection[i].dEdXSig_);
//		//--- f51 with slice dEdxSig ---//
//		  if(Preselection[i].dEdXSig_ >= 0.0 && Preselection[i].dEdXSig_ < 1.5)	     z.GetPlot(dEdXSig_slice1)->Fill(Preselection[i].f51_);
//		  if(Preselection[i].dEdXSig_ >= 1.5 && Preselection[i].dEdXSig_ < 3.0)	     z.GetPlot(dEdXSig_slice2)->Fill(Preselection[i].f51_);
// 		  if(Preselection[i].dEdXSig_ >= 3.0 && Preselection[i].dEdXSig_ < 4.5)      z.GetPlot(dEdXSig_slice3)->Fill(Preselection[i].f51_);
//		  if(Preselection[i].dEdXSig_ >= 4.5 && Preselection[i].dEdXSig_ < 6.0)	     z.GetPlot(dEdXSig_slice4)->Fill(Preselection[i].f51_);
//		  if(Preselection[i].dEdXSig_ >= 6.0 && Preselection[i].dEdXSig_ < 7.5)	     z.GetPlot(dEdXSig_slice5)->Fill(Preselection[i].f51_);
//		  if(Preselection[i].dEdXSig_ >= 7.5 && Preselection[i].dEdXSig_ < 9.0)	     z.GetPlot(dEdXSig_slice6)->Fill(Preselection[i].f51_);
//		  if(Preselection[i].dEdXSig_ >= 9.0 && Preselection[i].dEdXSig_ < 10.5)     z.GetPlot(dEdXSig_slice7)->Fill(Preselection[i].f51_);

		//--- Profile with f51(x-axis) dEdxSig(y-axis) slice in eta or energy ----//
		  if(TMath::Abs(Preselection[i].eta_) < 1.479)	x.GetProfile(EcalBarrel)->Fill(Preselection[i].f51_,Preselection[i].dEdXSig_);
		  if(TMath::Abs(Preselection[i].eta_) > 1.479 && TMath::Abs(Preselection[i].eta_) < 3.0)	x.GetProfile(EcalEndCup)->Fill(Preselection[i].f51_,Preselection[i].dEdXSig_);
		  if(TMath::Abs(Preselection[i].eta_) < 3.0)	x.GetProfile(EcalAll)->Fill(Preselection[i].f51_,Preselection[i].dEdXSig_);
		} 
	//	FillFlowHistogram(0,Preselection);
		FillFlowHistogram(0,TriggerSelection);

	}
	
  }
  void MonoCuts::FillFlowHistogram(int n, vector<MonoCandidate> CutFlowCand){
	PlotSet &z = CutFlow[n];
	  z.GetPlot(FracSatVNstrips)->Fill(CutFlowCand[0].subHits_,CutFlowCand[0].subSatHits_/CutFlowCand[0].subHits_);
	  z.GetPlot(DedXSig)->Fill(CutFlowCand[0].dEdXSig_);
	  z.GetPlot(XYPar0)->Fill(CutFlowCand[0].xyp0_);
	  z.GetPlot(XYPar1)->Fill(CutFlowCand[0].xyp1_);
	  z.GetPlot(XYPar2)->Fill(CutFlowCand[0].xyp2_);
	  z.GetPlot(RZPar0)->Fill(CutFlowCand[0].rzp0_);
	  z.GetPlot(RZPar1)->Fill(CutFlowCand[0].rzp1_);
	  z.GetPlot(RZcurv)->Fill(CutFlowCand[0].rzp2_);
	  z.GetPlot(F51)->Fill(CutFlowCand[0].f51_);
          z.GetPlot(ABCD)->Fill(CutFlowCand[0].f51_,CutFlowCand[0].dEdXSig_);
  }
  void MonoCuts::Clear(){
	TriggerSelection.clear();
	Preselection.clear();
  }

  void MonoCuts::WritePlots(TFile *oFile){
	oFile->cd(trigName_.c_str());
	CutFlow[0].WritePlot();
	Profile[0].WriteProfile();
	
  }

  //Cuts parameters
  const double MonoCuts::xyp0Cut_=0.6;
  const double MonoCuts::xyp2Cut_=1000;
  const double MonoCuts::rzp0Cut_=10;
  const double MonoCuts::rzp1Cut_=999;
  const double MonoCuts::rzp2Cut_=0.005;
  const double MonoCuts::distCut_ = 0.5;
  const double MonoCuts::hIsoCut_= 10;
  const double MonoCuts::dEdXSigCut_ = 0;//preselection cut
  const double MonoCuts::e55Cut_ = 200;
  const double MonoCuts::f51Cut_ = 0.85;

void BlindAnalyzer()
{
	TFile *oFile = new TFile("output/MonoData2018_Plot.root","recreate");
	TChain *Tree = new TChain("monopoles");
	Tree->Add("../../Data/Blind/BlindedData_2018AB.root");
	Tree->Add("../../Data/Blind/BlindedData_2018C.root");
	Tree->Add("../../Data/Blind/BlindedData_2018D.root");
//	TFile *fin = new TFile("../../Data/Blind/BlindedData_loose_2018all.root");
//        TTree *Tree = (TTree*)fin->Get("monopoles");

	int LastEvent = -1;
	vector<float> *SubHits = 0;   
	vector<float> *SatSubHits = 0;
	vector<float> *XYPar0 = 0;
	vector<float> *XYPar1 = 0;
	vector<float> *XYPar2 = 0;
	vector<float> *RZPar0 = 0;
	vector<float> *RZPar1 = 0;
	vector<float> *RZPar2 = 0;
	vector<float> *Eta = 0;
	vector<float> *seedFrac = 0;
	vector<float> *e55 = 0;
	vector<float> *Dist = 0;
	vector<float> *HIso = 0;
	vector<float> *dEdxSig = 0;   
	long long event;
	int run;
	int nPoint=0;
	int passHLT_Photon200;


	Tree->SetBranchAddress("SubHits",&SubHits);
	Tree->SetBranchAddress("SatSubHits",&SatSubHits);
	Tree->SetBranchAddress("XYPar0",&XYPar0);
	Tree->SetBranchAddress("XYPar1",&XYPar1);
	Tree->SetBranchAddress("XYPar2",&XYPar2);
	Tree->SetBranchAddress("RZPar0",&RZPar0);
	Tree->SetBranchAddress("RZPar1",&RZPar1);
	Tree->SetBranchAddress("RZPar2",&RZPar2);
	Tree->SetBranchAddress("Eta",&Eta);
	Tree->SetBranchAddress("seedFrac",&seedFrac);
	Tree->SetBranchAddress("e55",&e55);
	Tree->SetBranchAddress("Dist",&Dist);
	Tree->SetBranchAddress("HIso",&HIso);
	Tree->SetBranchAddress("dEdxSig",&dEdxSig);
	Tree->SetBranchAddress("Point_N",&nPoint);
	Tree->SetBranchAddress("event",&event);
//	Tree->SetBranchAddress("run",&run);
	Tree->SetBranchAddress("passHLT_Photon200", &passHLT_Photon200);
	
	vector<MonoCandidate> cand(10);
	MonoCuts TrgAnalysis("HLT_Photon",oFile);
	int nParticle = 0;
        for(int ev=0; ev<Tree->GetEntries();ev++){
		if(ev%1000000==0) cout<<ev<<"/"<<Tree->GetEntries()<<endl;
		Tree->GetEntry(ev);
		nParticle = nPoint;	
		if(nParticle > cand.size()) cand.resize(nParticle);
		for( int i=0;i < nParticle ;i++){
			if((*dEdxSig)[i]==-0) (*dEdxSig)[i]=0;
			cand[i] = MonoCandidate(
			(*SubHits)[i],
			(*SatSubHits)[i],
			(*Dist)[i],
			(*HIso)[i],
			(*XYPar0)[i],
			(*XYPar1)[i],
			(*XYPar2)[i],
			(*RZPar0)[i],
			(*RZPar1)[i],
			(*RZPar2)[i],
			(*e55)[i],
			(*seedFrac)[i],
			(*dEdxSig)[i], 
			(*Eta)[i],
			event
			);
			
		}
		TrgAnalysis.doAnalysis_data(cand,nParticle,passHLT_Photon200,ev);
	}
//	TrgAnalysis.SignalEff("HLT_Photon200");
	TrgAnalysis.WritePlots(oFile);

	oFile->Close();

}
