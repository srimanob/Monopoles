//
//	MonoAnalyzerPhoton.cc
//	Created by  Shih Lin
//	
//	Analysis code for Photon trigger(HLT_Photon175/200)
//
//      matching_option(bool):
//      0: NO truth matching
//      1: truth matching
//      analysis type (int):
//      0: MC
//      1: Delta-ray OFF
//      2: ECAL sys
//      3: DedxCrosstalk
//      
//      run:
//	root -l -q 'MonoAnalyzerPhoton.cc("year","mass",matching,analysis_type)'
//	ex. root -l -q 'MonoAnalyzerPhoton.cc("2018","1000",0,0)'
//
#include "iostream"
#include "TAttMarker.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "math.h"
#include <algorithm>
#include <string>
#include "../interface/Candidate.h"
#include "../interface/PlotSet.h"
#include "../interface/MonoCuts.h"
using namespace std;

void MonoCuts::doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, unsigned nCandidates,unsigned nPhoton, bool TRG, unsigned ev,bool matching_option,string year)
{
	Clear();

	for(unsigned c=0;c<nCandidates;c++){

		MonoCandidate &cands = cand[c];
		bool QualCut = evalQuality(cands);
		bool ECut = evalE(cands);
		bool ECut16 = evalE_16(cands);
		bool F51Cut = evalF51(cands);
		bool dEdXCut = evaldEdX(cands);


		//N-1 cut and relative efficiency
		if( ECut && F51Cut && dEdXCut &&TRG) N1CutCand_Qual.push_back(cands); 
		if( QualCut && F51Cut && dEdXCut &TRG) N1CutCand_Energy.push_back(cands);
		if( QualCut && ECut && dEdXCut &&TRG ) N1CutCand_F51.push_back(cands);
		if( QualCut && ECut && F51Cut &&TRG)  N1CutCand_Dedx.push_back(cands);
		if( QualCut && ECut & F51Cut && dEdXCut) N1CutCand_TRG.push_back(cands);

		//-----------------------------------------------------------------	
		//---Cutflow histograms--------------------------------------------
		//-----------------------------------------------------------------

		//signal efficiency

		//count for total events without TRG	
		if(TRG) CutFlowCand_TRG.push_back(cands);
		if(TRG && QualCut ) CutFlowCand_Qual.push_back(cands); 
		if(year == "2016" || year == "2016APV"){
			if(TRG && QualCut && ECut16 ) CutFlowCand_Energy.push_back(cands);
		}
		else{ 
			if(TRG && QualCut && ECut ) CutFlowCand_Energy.push_back(cands);
		}

	}//for cand loop

	sort(CutFlowCand_TRG.begin(),CutFlowCand_TRG.begin()+CutFlowCand_TRG.size());
	if(CutFlowCand_TRG.size()>0) 
	{
		count++;
		FillNoCutHistogram(0,CutFlowCand_TRG,matching_option);
	}
	sort(CutFlowCand_Qual.begin(),CutFlowCand_Qual.begin()+CutFlowCand_Qual.size());
	if(CutFlowCand_Qual.size()>0) 
	{
		Qual_count++;	
		FillFlowHistogram(0,CutFlowCand_Qual,matching_option);
	}

	sort(CutFlowCand_Energy.begin(),CutFlowCand_Energy.begin()+CutFlowCand_Energy.size());
	if(CutFlowCand_Energy.size()>0)
	{
		E_count++;	
		FillFlowHistogram(1,CutFlowCand_Energy,matching_option);

		MonoCandidate &Cand = CutFlowCand_Energy[0];
		bool F51Cut = evalF51(Cand);
		if(F51Cut) CutFlowCand_F51.push_back(Cand);
	}
	sort(CutFlowCand_F51.begin(),CutFlowCand_F51.begin()+CutFlowCand_F51.size());
	if(CutFlowCand_F51.size()>0)
	{
		f51_count++;	
		FillFlowHistogram(2,CutFlowCand_F51,matching_option);

		MonoCandidate &SelectedCand = CutFlowCand_F51[0];
		bool dEdXCut = evaldEdX(SelectedCand);
		if(dEdXCut) CutFlowCand_Dedx.push_back(SelectedCand);

	}
	sort(CutFlowCand_Dedx.begin(),CutFlowCand_Dedx.begin()+CutFlowCand_Dedx.size());
	if(CutFlowCand_Dedx.size()>0)
	{
		dEdX_count++;	
		FillFlowHistogram(3,CutFlowCand_Dedx,matching_option);
	}




	///////////////////////////////////////////////////
	/////////  N1 Cut Plots and  Count   //////////////
	///////////////////////////////////////////////////


	//count n_1Plot TRG+Qual and get plots
	sort(N1CutCand_Qual.begin(),N1CutCand_Qual.begin()+N1CutCand_Qual.size());
	if(N1CutCand_Qual.size()>0) 
	{
		FillN1Histogram(0,N1CutCand_Qual);
		NoQual++;
	}
	sort(N1CutCand_Energy.begin(),N1CutCand_Energy.begin()+N1CutCand_Energy.size());
	if(N1CutCand_Energy.size()>0)
	{
		FillN1Histogram(1,N1CutCand_Energy);
		NoE++;	
	}
	sort(N1CutCand_F51.begin(),N1CutCand_F51.begin()+N1CutCand_F51.size());
	if(N1CutCand_F51.size()>0)
	{
		FillN1Histogram(2,N1CutCand_F51);
		NoF51++;	
	}
	sort(N1CutCand_Dedx.begin(),N1CutCand_Dedx.begin()+N1CutCand_Dedx.size());
	if(N1CutCand_Dedx.size()>0)
	{

		FillN1Histogram(3,N1CutCand_Dedx);
		NodEdXCut++;	
	}
	sort(N1CutCand_TRG.begin(),N1CutCand_TRG.begin()+N1CutCand_TRG.size());
	if(N1CutCand_TRG.size()>0) 
	{
		FillN1Histogram(4,N1CutCand_TRG);
		NoTRG++;
	}
}


void MonoCuts::FillNoCutHistogram(int n,vector<MonoCandidate> Cand, bool matching){

	PlotSet &z = NoCutPlot[n];
	PlotSet &x = NoCutProfile[n];
	vector<MonoCandidate> Matched;
	if (matching == 1){
		Matched = Matching(Cand);
		if(Matched.size() != 0){	
			z.GetPlot(FracSatVNstrips)->Fill(Matched[0].subHits_,Matched[0].subSatHits_/Matched[0].subHits_);
			z.GetPlot(DedXSig)->Fill(Matched[0].dEdXSig_);
			z.GetPlot(XYPar0)->Fill(Matched[0].xyp0_);
			z.GetPlot(XYPar1)->Fill(Matched[0].xyp1_);
			z.GetPlot(XYPar2)->Fill(Matched[0].xyp2_);
			z.GetPlot(RZPar0)->Fill(Matched[0].rzp0_);
			z.GetPlot(RZPar1)->Fill(Matched[0].rzp1_);
			z.GetPlot(RZcurv)->Fill(Matched[0].rzp2_);
			z.GetPlot(E55)->Fill(Matched[0].e55_);
			z.GetPlot(F51)->Fill(Matched[0].f51_);
			z.GetPlot(HcalIso)->Fill(Matched[0].hIso_);
			z.GetPlot(ABCD)->Fill(Matched[0].f51_,Matched[0].dEdXSig_);

		//	for(int i=0; i < Matched.size() ;i++){
		//		x.GetProfile(PileUp_f51)->Fill(Matched[i].NPV_,Matched[i].f51_);
		//		x.GetProfile(PileUp_DedXSig)->Fill(Matched[i].NPV_,Matched[i].f51_);
		//		if(TMath::Abs(Matched[i].eta_) < 1.479)	  x.GetProfile(EcalBarrel)->Fill(Matched[i].f51_,Matched[i].dEdXSig_);
		//		if(TMath::Abs(Matched[i].eta_) > 1.479 && TMath::Abs(Matched[i].eta_) < 3.0) 	  x.GetProfile(EcalEndCup)->Fill(Matched[i].f51_,Matched[i].dEdXSig_);
		//		if(TMath::Abs(Matched[i].eta_) < 3.0 ) x.GetProfile(EcalAll)->Fill(Matched[i].f51_,Matched[i].dEdXSig_);
		//	}
		}
	}
	else{
		for(int i=0; i < Cand.size() ;i++){
			z.GetPlot(FracSatVNstrips)->Fill(Cand[i].subHits_,Cand[i].subSatHits_/Cand[i].subHits_);
			z.GetPlot(DedXSig)->Fill(Cand[i].dEdXSig_);
			z.GetPlot(XYPar0)->Fill(Cand[i].xyp0_);
			z.GetPlot(XYPar1)->Fill(Cand[i].xyp1_);
			z.GetPlot(XYPar2)->Fill(Cand[i].xyp2_);
			z.GetPlot(RZPar0)->Fill(Cand[i].rzp0_);
			z.GetPlot(RZPar1)->Fill(Cand[i].rzp1_);
			z.GetPlot(RZcurv)->Fill(Cand[i].rzp2_);
			z.GetPlot(E55)->Fill(Cand[i].e55_);
			z.GetPlot(F51)->Fill(Cand[i].f51_);
			z.GetPlot(HcalIso)->Fill(Cand[i].hIso_);
			z.GetPlot(ABCD)->Fill(Cand[i].f51_,Cand[i].dEdXSig_);
		//	x.GetProfile(PileUp_f51)->Fill(Cand[i].NPV_,Cand[i].f51_);
		//	x.GetProfile(PileUp_DedXSig)->Fill(Cand[i].NPV_,Cand[i].f51_);
		//	if(TMath::Abs(Cand[i].eta_) < 1.479)	  x.GetProfile(EcalBarrel)->Fill(Cand[i].f51_,Cand[i].dEdXSig_);
		//	if(TMath::Abs(Cand[i].eta_) > 1.479 && TMath::Abs(Cand[i].eta_) < 3.0) 	  x.GetProfile(EcalEndCup)->Fill(Cand[i].f51_,Cand[i].dEdXSig_);
		//	if(TMath::Abs(Cand[i].eta_) < 3.0 ) x.GetProfile(EcalAll)->Fill(Cand[i].f51_,Cand[i].dEdXSig_);
		}
	}
	Matched.clear();
}


void MonoCuts::FillFlowHistogram(int n, vector<MonoCandidate> CutFlowCand, bool matching){

	PlotSet &z = CutFlow[n];
	PlotSet &x = Profile[n];
	vector<MonoCandidate> Matched;
	if (matching == 1){
		Matched = Matching(CutFlowCand);	
		if(Matched.size() != 0){	
			z.GetPlot(FracSatVNstrips)->Fill(Matched[0].subHits_,Matched[0].subSatHits_/Matched[0].subHits_);
			z.GetPlot(DedXSig)->Fill(Matched[0].dEdXSig_);
			z.GetPlot(XYPar0)->Fill(Matched[0].xyp1_);
			z.GetPlot(XYPar1)->Fill(Matched[0].xyp0_);
			z.GetPlot(XYPar2)->Fill(Matched[0].xyp2_);
			z.GetPlot(RZPar0)->Fill(Matched[0].rzp0_);
			z.GetPlot(RZPar1)->Fill(Matched[0].rzp1_);
			z.GetPlot(RZcurv)->Fill(Matched[0].rzp2_);
			z.GetPlot(E55)->Fill(Matched[0].e55_);
			z.GetPlot(F51)->Fill(Matched[0].f51_);
			z.GetPlot(HcalIso)->Fill(Matched[0].hIso_);
			z.GetPlot(ABCD)->Fill(Matched[0].f51_,Matched[0].dEdXSig_);

		//	for(int i=0; i < Matched.size() ;i++){
		//		x.GetProfile(PileUp_f51)->Fill(Matched[i].NPV_,Matched[i].f51_);
		//		x.GetProfile(PileUp_DedXSig)->Fill(Matched[i].NPV_,Matched[i].dEdXSig_);
		//		if(TMath::Abs(Matched[i].eta_) < 1.479)	  x.GetProfile(EcalBarrel)->Fill(Matched[i].f51_,Matched[i].dEdXSig_);
		//		if(TMath::Abs(Matched[i].eta_) > 1.479 && TMath::Abs(Matched[i].eta_) < 3.0) 	  x.GetProfile(EcalEndCup)->Fill(Matched[i].f51_,Matched[i].dEdXSig_);
		//		if(TMath::Abs(Matched[i].eta_) < 3.0 ) x.GetProfile(EcalAll)->Fill(Matched[i].f51_,Matched[i].dEdXSig_);
		//	}
		}
	}
	else{
		for(int i=0; i < CutFlowCand.size() ;i++){
			z.GetPlot(FracSatVNstrips)->Fill(CutFlowCand[i].subHits_,CutFlowCand[i].subSatHits_/CutFlowCand[i].subHits_);
			z.GetPlot(DedXSig)->Fill(CutFlowCand[i].dEdXSig_);
			z.GetPlot(XYPar0)->Fill(CutFlowCand[i].xyp0_);
			z.GetPlot(XYPar1)->Fill(CutFlowCand[i].xyp1_);
			z.GetPlot(XYPar2)->Fill(CutFlowCand[i].xyp2_);
			z.GetPlot(RZPar0)->Fill(CutFlowCand[i].rzp0_);
			z.GetPlot(RZPar1)->Fill(CutFlowCand[i].rzp1_);
			z.GetPlot(RZcurv)->Fill(CutFlowCand[i].rzp2_);
			z.GetPlot(E55)->Fill(CutFlowCand[i].e55_);
			z.GetPlot(F51)->Fill(CutFlowCand[i].f51_);
			z.GetPlot(HcalIso)->Fill(CutFlowCand[i].hIso_);
			z.GetPlot(ABCD)->Fill(CutFlowCand[0].f51_,CutFlowCand[0].dEdXSig_);

		//	x.GetProfile(PileUp_f51)->Fill(CutFlowCand[i].NPV_,CutFlowCand[i].f51_);
		//	x.GetProfile(PileUp_DedXSig)->Fill(CutFlowCand[i].NPV_,CutFlowCand[i].dEdXSig_);
		//	if(TMath::Abs(CutFlowCand[i].eta_) < 1.479)	  x.GetProfile(EcalBarrel)->Fill(CutFlowCand[i].f51_,CutFlowCand[i].dEdXSig_);
		//	if(TMath::Abs(CutFlowCand[i].eta_) > 1.479 && TMath::Abs(CutFlowCand[i].eta_) < 3.0) 	  x.GetProfile(EcalEndCup)->Fill(CutFlowCand[i].f51_,CutFlowCand[i].dEdXSig_);
		//	if(TMath::Abs(CutFlowCand[i].eta_) < 3.0 ) x.GetProfile(EcalAll)->Fill(CutFlowCand[i].f51_,CutFlowCand[i].dEdXSig_);
		}
	}
	Matched.clear();
}
void MonoCuts::FillN1Histogram(int n, vector<MonoCandidate> N1CutCand){
	PlotSet &z = n_1Plot[n];
	for(int i=0; i < N1CutCand.size() ;i++){
		z.GetPlot(FracSatVNstrips)->Fill(N1CutCand[i].subHits_,N1CutCand[i].subSatHits_/N1CutCand[i].subHits_);
		z.GetPlot(DedXSig)->Fill(N1CutCand[i].dEdXSig_);
		z.GetPlot(RZcurv)->Fill(N1CutCand[i].rzp2_);
		z.GetPlot(E55)->Fill(N1CutCand[i].e55_);
		z.GetPlot(F51)->Fill(N1CutCand[i].f51_);
		z.GetPlot(HcalIso)->Fill(N1CutCand[i].hIso_);
		z.GetPlot(ABCD)->Fill(N1CutCand[i].f51_,N1CutCand[i].dEdXSig_);
	}
}
vector<MonoCandidate> MonoCuts::Matching(vector<MonoCandidate> Cand){


	for(int i=0; i<Cand.size();i++){
		double m_deltaR=0;
		double am_deltaR=0;
		m_deltaR = sqrt(pow(Cand[i].eta_-Cand[0].mono_eta_,2)+
				pow(Cand[i].phi_-Cand[0].mono_phi_,2));
		am_deltaR= sqrt(pow(Cand[i].eta_-Cand[0].amon_eta_,2)+
				pow(Cand[i].phi_-Cand[0].amon_phi_,2));

		if(m_deltaR<0.15||am_deltaR<0.15){
			Matched.push_back(Cand[i]);		
		}

	}

	return Matched;
}
void MonoCuts::Clear(){

	CutFlowCand_TRG.clear();	
	CutFlowCand_Qual.clear();
	CutFlowCand_Energy.clear();
	CutFlowCand_F51.clear();
	CutFlowCand_Dedx.clear();
	CutFlowCand_looseF51.clear();
	CutFlowCand_looseDedx.clear();
	N1CutCand_TRG.clear();	
	N1CutCand_Qual.clear();
	N1CutCand_Energy.clear();
	N1CutCand_F51.clear();
	N1CutCand_Dedx.clear();
	Matched.clear();
	HighPtPhoton.clear();
}

void MonoCuts::WritePlots(TFile *oFile){
	oFile->cd(trigName_.c_str());
	NoCutPlot[0].WritePlot();
	//NoCutProfile[0].WriteProfile();
	for(unsigned c=0; c<nCut; c++) n_1Plot[c].WritePlot();

	for(unsigned c=0; c<nCut; c++){
		CutFlow[c].WritePlot();
		//Profile[c].WriteProfile();
	}

}
void MonoCuts::SignalEff(const string trName, double TotalEvents)
{
	//signal efficiency = no. of events after all selection cuts/all events
	cout<<trName<<" ================================="<<endl;
	cout<<"        TRG "<<count<<endl;
	cout<<"QualityCuts "<<Qual_count<<endl;
	cout<<"       ECut "<<E_count<<endl;
	cout<<"     F51Cut "<<f51_count<<endl;
	cout<<" dEdXSigCut "<<dEdX_count<<endl;
	cout<<"------------------------------------------"<<endl;
	cout<<"Signal efficiency = "<<(double)dEdX_count/(double)TotalEvents<<endl;
	cout<<endl;
	cout<<"Relative effciency count"<<endl;
	cout<<"------------------------------------------"<<endl;
	cout<<"     No TRG "<<NoTRG<<endl;
	cout<<" No Quality "<<NoQual<<endl;
	cout<<"    No ECut "<< NoE <<endl;
	cout<<"  No F51Cut "<<NoF51<<endl;
	cout<<" No dEdXSig "<<NodEdXCut<<endl;
	cout<<"------------------------------------------"<<endl;
	cout<<"Relative efficiency"<<endl;
	cout<<"------------------------------------------"<<endl;
	cout<<"        TRG "<<(double)dEdX_count/(double)NoTRG<<endl;
	cout<<"QualityCuts "<<(double)dEdX_count/(double)NoQual<<endl;
	cout<<"       ECut "<<(double)dEdX_count/(double)NoE<<endl;
	cout<<"     F51Cut "<<(double)dEdX_count/(double)NoF51<<endl;
	cout<<" dEdXSigCut "<<(double)dEdX_count/(double)NodEdXCut<<endl;
	cout<<endl;

	cout<<endl;
}

//Cuts parameters
const double MonoCuts::xyp0Cut_=0.6;
const double MonoCuts::xyp2Cut_=1000;
const double MonoCuts::rzp0Cut_=10;
const double MonoCuts::rzp1Cut_=999;
const double MonoCuts::rzp2Cut_=0.005;
const double MonoCuts::distCut_ = 0.5;
const double MonoCuts::hIsoCut_= 10;
const double MonoCuts::dEdXSigCut_ = 9;
const double MonoCuts::e55Cut_ = 200;
const double MonoCuts::e55Cut2016_ = 175;
const double MonoCuts::f51Cut_ = 0.85;
const double MonoCuts::photonCut_ = 200;
const double MonoCuts::dEdXSig_looseCut_ = 7;
const double MonoCuts::f51_looseCut_= 0.75 ;

void MonoAnalyzerPhoton(string year, string mass,bool matching_option, int sys_option=0)
{
	string matching;
	if(matching_option == 0){
		matching = "0";
	}
	else{
		matching = "1";
	}

	TChain *tree = new TChain("monopoles");
	string sys;
	//tree->Add(("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/MCNtuple"+year+"/spin0/*").c_str());
	//tree->Add(("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/MCNtuple"+year+"/2018-1000-SpinZero.root").c_str());
	//
	// remember to change the path befor you run
	if(sys_option == 0){
		sys = "";
		tree->Add(("/eos/cms/store/user/srimanob/monopole/13TeV/Legacy-NTUPLE-v2/merges/"+year+"-"+mass+".root").c_str());
	}
	else if(sys_option == 1){
		sys = "DeltaRayOff";
		cout<<"Processing for the DeltaRayOff...."<<endl;
		tree->Add(("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/Systematic/DeltaRayOff/"+year+"/"+mass+"/*.root").c_str());
	}
	else if(sys_option == 2){
		sys = "SpikeAlgo";
		cout<<"Processing for the Ecal...."<<endl;
		tree->Add(("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/Systematic/ECal/"+year+"/"+mass+"/*.root").c_str());
	}
	else if(sys_option == 3){
		sys = "CrossTalk";
		cout<<"Processing for the CrossTalk...."<<endl;
		tree->Add(("/wk_cms2/shihlin0314/CMSSW_8_0_29/src/Systematic/DedxCrossTalk/"+year+"/"+mass+"/*.root").c_str());
	}

	TFile *oFile = new TFile(("output/MonoPhotonAnalysis_"+year+"_"+mass+"_"+sys+"_"+matching+".root").c_str(),"recreate");

	Bool_t passHLT_Photon200;
	Bool_t passHLT_Photon175;
	Bool_t passHLT_DoublePhoton70;
	Bool_t passHLT_DoublePhoton60;
	unsigned nCandidates;
	unsigned event;
	unsigned NPV;
	double mono_eta;
	double mono_phi;
	double amon_eta;
	double amon_phi;
	vector<double> *subHits=0;
	vector<double> *subSatHits=0;
	vector<double> *dEdXSig=0;
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
	vector<double> * hIso = 0;
	vector<double> * eta = 0;
	vector<double> * phi = 0;
	vector<double> * pho_eta = 0;
	vector<double> * pho_phi = 0;
	vector<double> * pho_pt = 0;
	unsigned nPhoton;

	vector<double> * Cross = 0;

	tree->SetBranchAddress("event",&event);
	tree->SetBranchAddress("NPV",&NPV);
	tree->SetBranchAddress("passHLT_Photon200" , &passHLT_Photon200);
	tree->SetBranchAddress("passHLT_Photon175" , &passHLT_Photon175);
	tree->SetBranchAddress("passHLT_DoublePhoton70",&passHLT_DoublePhoton70);
	tree->SetBranchAddress("passHLT_DoublePhoton60",&passHLT_DoublePhoton60);
	tree->SetBranchAddress("cand_N",&nCandidates);
	tree->SetBranchAddress("cand_SubHits",&subHits);
	tree->SetBranchAddress("cand_SatSubHits",&subSatHits);
	tree->SetBranchAddress("cand_dEdXSig",&dEdXSig);
	tree->SetBranchAddress("cand_TIso",&tIso);
	tree->SetBranchAddress("cand_f51",&f51);
	tree->SetBranchAddress("cand_f15",&f15);
	tree->SetBranchAddress("cand_e55",&e55);
	tree->SetBranchAddress("cand_HIso",&hIso);
	tree->SetBranchAddress("cand_XYPar0",&xyp0);
	tree->SetBranchAddress("cand_XYPar1",&xyp1);
	tree->SetBranchAddress("cand_XYPar2",&xyp2);
	tree->SetBranchAddress("cand_RZPar0",&rzp0);
	tree->SetBranchAddress("cand_RZPar1",&rzp1);
	tree->SetBranchAddress("cand_RZPar2",&rzp2);
	tree->SetBranchAddress("cand_eta",&eta);
	tree->SetBranchAddress("cand_phi",&phi);
	tree->SetBranchAddress("cand_dist",&dist);
	tree->SetBranchAddress("mono_eta",&mono_eta);
	tree->SetBranchAddress("mono_phi",&mono_phi);
	tree->SetBranchAddress("amon_eta",&amon_eta);
	tree->SetBranchAddress("amon_phi",&amon_phi);
	tree->SetBranchAddress("pho_N",&nPhoton);
	tree->SetBranchAddress("pho_eta",&pho_eta);
	tree->SetBranchAddress("pho_phi",&pho_phi);
	tree->SetBranchAddress("pho_pt",&pho_pt);
	//tree->SetBranchAddress("cand_SwissCross",&Cross);


	const unsigned NEvents = tree->GetEntries();


	MonoCuts noTrgAnalysis("NOTRG",oFile);
	MonoCuts TrgAnalysis("HLT_Photon200",oFile);

	vector<MonoCandidate> cand(10);	
	vector<Photon> photon(0);

	for(unsigned ev=0; ev<NEvents;ev++){
		if(ev%1000==0)	cout<<ev<<"/"<<NEvents<<endl;
		tree->GetEntry(ev);

		if(nCandidates>cand.size()) cand.resize(nCandidates);

		if(nPhoton>photon.size()) photon.resize(nPhoton);

		for(unsigned i=0;i<nCandidates;i++){

			cand[i] = MonoCandidate(
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
					//(*Cross)[i],
					(*e55)[i],
					(*hIso)[i],
					(*eta)[i],
					(*phi)[i],
					mono_eta,
					mono_phi,
					amon_eta,
					amon_phi,
					event,
					NPV	
						);
		}
		if(nPhoton!=0){
			for(unsigned j=0;j<nPhoton;j++){
				photon[j] = Photon(
						(*pho_eta)[j],
						(*pho_phi)[j],
						(*pho_pt)[j]
						);
			}
		}
		noTrgAnalysis.doAnalysis(cand,photon,nCandidates,nPhoton,true,ev,matching_option,year);		
		if( year == "2016" || year == "2016APV")      TrgAnalysis.doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_Photon175,ev,matching_option,year);
		else      TrgAnalysis.doAnalysis(cand,photon,nCandidates,nPhoton,passHLT_Photon200,ev,matching_option,year);
	}//for every event loop

	noTrgAnalysis.WritePlots(oFile);
	noTrgAnalysis.SignalEff("NOTRG",NEvents);
	noTrgAnalysis.SaveAs_csv(("output/csv_file/Signaleff_"+year+"_"+mass+"_"+sys+"_"+matching+".csv").c_str(),NEvents,mass,"NoTrg");
	TrgAnalysis.WritePlots(oFile);
	if(year == "2016" || year == "2016APV"){
		TrgAnalysis.SignalEff("HLT_Photon175",NEvents);
		TrgAnalysis.SaveAs_csv(("output/csv_file/Signaleff_"+year+"_"+mass+"_HLT_"+sys+"_"+matching+".csv").c_str(),NEvents,mass,"Photon175");
	}
	else{
		TrgAnalysis.SignalEff("HLT_Photon200",NEvents);
		TrgAnalysis.SaveAs_csv(("output/csv_file/Signaleff_"+year+"_"+mass+"_HLT_"+sys+"_"+matching+".csv").c_str(),NEvents,mass,"Photon200");
	}
	oFile->Close();	
}
