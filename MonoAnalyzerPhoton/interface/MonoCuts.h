#ifndef _MONOCUTS_H_
#define _MONOCUTS_H_
#include "PlotSet.h"
class MonoCuts:public MonoCandidate, public Photon
{
public:
  MonoCuts(){}

  MonoCuts(const string &trName,TFile *openFile):trigName_(trName){
	//create a new directory in the output file
	openFile->mkdir(trName.c_str());

	//create histrogram N-1 cutf for No trigger
	NoCutPlot.resize(1U);
	PlotSet &x = NoCutPlot[0];
        x.CreatPlot(FracSatVNstrips,new TH2D("FracSatVNstrips","",100,0,1000,100,0,1));
        x.CreatPlot(DedXSig,new TH1D("DedXSig","",100,0,30));
        x.CreatPlot(XYPar0,new TH1D("XYPar0","",50,-1,1));
        x.CreatPlot(XYPar1,new TH1D("XYPar1","",100,-10,10));
        x.CreatPlot(XYPar2,new TH1D("XYPar2","",100,-20000,20000));
        x.CreatPlot(RZPar0,new TH1D("RZPar0","",100,-20,20));
        x.CreatPlot(RZPar1,new TH1D("RZPar1","",100,-15,15));
        x.CreatPlot(RZcurv,new TH1D("RZcurv","",100,-0.01,0.01));
	x.CreatPlot(E55,new TH1D("E55","",100,-1,1200));
        x.CreatPlot(F51,new TH1D("F51","",100,0.2,1.1));
        x.CreatPlot(HcalIso,new TH1D("HcalIso","",100,-1,10));
        x.CreatPlot(ABCD,new TH2D("ABCD","",100,0,1.1,100,0,30));

//        NoCutProfile.resize(1U);
//        PlotSet &p = NoCutProfile[0];
//        p.CreatProfile(EcalBarrel,new TProfile("EcalBarrel","",30,0,1.1,0,30));
//        p.CreatProfile(EcalEndCup,new TProfile("EcalEndCup","",30,0,1.1,0,30));
//        p.CreatProfile(EcalAll ,new TProfile("EcalAll","",30,0,1.1,0,30));
//        p.CreatProfile(PileUp_DedXSig ,new TProfile("PileUp_DedXSig","",6,0,60,-1,30));
//        p.CreatProfile(PileUp_f51 ,new TProfile("PileUp_f51","",6,0,60,-0.1,1));
//
 	cutName_[0] = "Quality_";
 	cutName_[1] = "Energy_";
	cutName_[2] = "F51_";
	cutName_[3] = "dEdxSig_";
	cutName_[4] = "HLT_";
	
        n_1Plot.resize(nCut);

	for( int c = 0; c <nCut ;c++){

	   string cutn1name = "N1_"+ cutName_[c];
	   PlotSet &z = n_1Plot[c];
	   z.CreatPlot(FracSatVNstrips,new TH2D((cutn1name+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
	   z.CreatPlot(DedXSig,new TH1D((cutn1name+"DedXSig").c_str(),"",100,0,30));
//           z.CreatPlot(XYPar0,new TH1D((cutn1name+"XYPar0").c_str(),"",100,-1,1));
//           z.CreatPlot(XYPar1,new TH1D((cutn1name+"XYPar1").c_str(),"",100,-10,10));
//           z.CreatPlot(XYPar2,new TH1D((cutn1name+"XYPar2").c_str(),"",100,-20000,20000));
//           z.CreatPlot(RZPar0,new TH1D((cutn1name+"RZPar0").c_str(),"",100,-20,20));
//           z.CreatPlot(RZPar1,new TH1D((cutn1name+"RZPar1").c_str(),"",100,-15,15));
	   z.CreatPlot(RZcurv,new TH1D((cutn1name+"RZcurv").c_str(),"",100,-0.01,0.01));
           z.CreatPlot(E55,new TH1D((cutn1name+"E55").c_str(),"",100,-1,1200));
           z.CreatPlot(F51,new TH1D((cutn1name+"F51").c_str(),"",100,0.2,1.1));
           z.CreatPlot(HcalIso,new TH1D((cutn1name+"HcalIso").c_str(),"",100,-1,10));
           z.CreatPlot(ABCD,new TH2D((cutn1name+"ABCD").c_str(),"",100,0,1.1,100,0,30));
	}

	CutFlow.resize(10); 
        Profile.resize(100);
	for( int c = 0;c<nCut-1;c++){
	   PlotSet &y = CutFlow[c];
	   string cutflowName = "Flow_"+cutName_[c];
           y.CreatPlot(FracSatVNstrips,new TH2D((cutflowName+"FracSatVNstrips").c_str(),"",100,0,1000,100,0,1));
           y.CreatPlot(DedXSig,new TH1D((cutflowName+"DedXSig").c_str(),"",100,0,30));
           y.CreatPlot(XYPar0,new TH1D((cutflowName+"XYPar0").c_str(),"",100,-1,1));
           y.CreatPlot(XYPar1,new TH1D((cutflowName+"XYPar1").c_str(),"",100,-10,10));
           y.CreatPlot(XYPar2,new TH1D((cutflowName+"XYPar2").c_str(),"",100,-20000,20000));
           y.CreatPlot(RZPar0,new TH1D((cutflowName+"RZPar0").c_str(),"",100,-20,20));
           y.CreatPlot(RZPar1,new TH1D((cutflowName+"RZPar1").c_str(),"",100,-15,15));
           y.CreatPlot(RZcurv,new TH1D((cutflowName+"RZcurv").c_str(),"",100,-0.01,0.01));
           y.CreatPlot(E55,new TH1D((cutflowName+"E55").c_str(),"",100,-1,1200));
           y.CreatPlot(F51,new TH1D((cutflowName+"F51").c_str(),"",100,0.2,1.1));
           y.CreatPlot(HcalIso,new TH1D((cutflowName+"HcalIso").c_str(),"",100,-1,10));
           y.CreatPlot(ABCD,new TH2D((cutflowName+"ABCD").c_str(),"",100,0,1.1,100,0,30));

       	   PlotSet &w = Profile[c];
       	   w.CreatProfile(EcalBarrel ,new TProfile((cutflowName+"EcalBarrel").c_str(),"",30,0,1.1,0,30));
       	   w.CreatProfile(EcalEndCup ,new TProfile((cutflowName+"EcalEndCup").c_str(),"",30,0,1.1,0,30));
       	   w.CreatProfile(EcalAll ,new TProfile((cutflowName+"EcalAll").c_str(),"",30,0,1.1,0,30));
           w.CreatProfile(PileUp_DedXSig ,new TProfile((cutflowName+"PileUp_DedXSig").c_str(),"",6,0,60,-1,30));
           w.CreatProfile(PileUp_f51 ,new TProfile((cutflowName+"PileUp_f51").c_str(),"",6,0,60,-0.1,1));

	}
  }
  ~MonoCuts(){}

  void doAnalysis(vector<MonoCandidate> &cand, vector<Photon> & pho, unsigned nCandidates,unsigned nPhoton, bool TRG, unsigned ev,bool matching_option,string year);
  void doAnalysis_data(vector<MonoCandidate> &cand,unsigned nParticle,bool passHLT_Photon200,unsigned ev);
  void FillNoCutHistogram(int n,vector<MonoCandidate> Cand,bool matching);
  void FillFlowHistogram(int n, vector<MonoCandidate> CutFlowCand,bool matching);
  void FillN1Histogram(int n, vector<MonoCandidate> N1CutCand);
  vector<MonoCandidate>  Matching(vector<MonoCandidate> CutFlowCand);
  void Clear();
  void WritePlots(TFile *oFile);
  void PrintInfo(vector<MonoCandidate> Cand){	
	cout<<"event tag"<<Cand[0].event_<<endl;
	for(int i=0;i<Cand.size();i++){
		double m_deltaR =0;
		double am_deltaR=0;
                m_deltaR = sqrt(pow(Cand[i].eta_-Cand[0].mono_eta_,2)+
                                pow(Cand[i].phi_-Cand[0].mono_phi_,2));
                am_deltaR= sqrt(pow(Cand[i].eta_-Cand[0].amon_eta_,2)+
                                pow(Cand[i].phi_-Cand[0].amon_phi_,2));

    	        cout<<"          candidate           monoGen         antimonoGen"<<endl;
    	   	cout<<"eta      "<<setprecision(5)<<Cand[i].eta_<<setw(20)<<Cand[i].mono_eta_
    	   			 <<setw(20)<<Cand[i].amon_eta_<<endl;
    	   	cout<<"phi      "<<setprecision(5)<<Cand[i].phi_<<setw(20)<<Cand[i].mono_phi_
    	   			 <<setw(20)<<Cand[i].amon_phi_<<endl;
    	   	cout<<"m deltaR "<<m_deltaR<<endl;
    	   	cout<<"a deltaR "<<am_deltaR<<endl;
		cout<<i<<endl;
		cout<<"     eta "<<setprecision(5)<<Cand[i].eta_<<endl;
		cout<<"     phi "<<setprecision(5)<<Cand[i].phi_<<endl;
		cout<<"     E55 "<<setprecision(5)<<Cand[i].e55_<<endl;
		cout<<" dEdxSig "<<setprecision(5)<<Cand[i].dEdXSig_<<endl;
		cout<<"     f51 "<<setprecision(5)<<Cand[i].f51_<<endl;
		cout<<"   Swiss "<<setprecision(5)<<Cand[i].Cross_<<endl;
		cout<<"----------------------------"<<endl;
	}
  }
  void PrintPhoton(){	
	//Just Print out all photon info in  one event
	for(int j=0;j<HighPtPhoton.size();j++){
		cout<<"pho  eta "<<HighPtPhoton[j].pho_eta_<<endl;
        	cout<<"pho  phi "<<HighPtPhoton[j].pho_phi_<<endl;
        	cout<<"pho   pt "<<HighPtPhoton[j].pho_pt_<<endl;
	}

	//define: photon-monopole angle<0.15 is photon-like
	for(int i=0;i<HighPtPhoton.size();i++){
		double photonMono_deltaR = 0;
	//	photonMono_deltaR = sqrt(pow(CutFlowCand_Dedx[i].eta_-HighPtPhoton[j].pho_eta_,2)+
        //	               	         pow(CutFlowCand_Dedx[i].phi_-HighPtPhoton[j].pho_phi_,2));
		if(photonMono_deltaR<0.15){
			cout<<"deltaR "<<photonMono_deltaR<<endl;
			cout<<"photon-like monopole"<<endl;
			photonLike++;// otherwise not photon-like monopole
		}
	}
	cout<<endl;
	cout<<"=================================="<<endl;
  }
  void SignalEff(const string trName, double TotalEvents);

  void SaveAs_csv_crosstalk(const string fileName, double TotalEvents, string mass,string X_option, string UD_option);
  void SaveAs_csv(const string fileName, double TotalEvents, string mass,string trName){
        ofstream fout(fileName);

        if(!fout) cout<<"Cannot open file"<<endl;
        else    cout<<"Open file successfully"<<endl;
	fout<<mass+" "+trName<<",CutFlow,"<<endl;
        fout<<"Generated ev,"<<TotalEvents<<endl;
        fout<<"        TRG, "<<count<<endl;
        fout<<"QualityCuts, "<<Qual_count<<endl;
        fout<<"       ECut, "<<E_count<<endl;
        fout<<"     F51Cut, "<<f51_count<<endl;
        fout<<" dEdXSigCut, "<<dEdX_count<<endl;
        fout<<"Signal efficiency, "<<(double)dEdX_count/(double)TotalEvents<<endl;
        fout<<endl;
	fout<<",N1Cuts,Relative eff"<<endl;
        fout<<"     No TRG, "<<NoTRG <<","<<(double)dEdX_count/(double)NoTRG<<endl;
        fout<<" No Quality, "<<NoQual<<","<<(double)dEdX_count/(double)NoQual<<endl;
        fout<<"    No ECut, "<< NoE  <<","<<(double)dEdX_count/(double)NoE<<endl;
        fout<<"  No F51Cut, "<<NoF51 <<","<<(double)dEdX_count/(double)NoF51<<endl;
        fout<<" No dEdXSig, "<<NodEdXCut<<","<<(double)dEdX_count/(double)NodEdXCut<<endl;
        fout<<endl;

  }
  bool operator<(const MonoCandidate &mc)const{
   if(dEdXSig_>mc.dEdXSig_) return true;
   else if(dEdXSig_==mc.dEdXSig_){
        if(f51_>mc.f51_) return true;
        else return false;
        }
    else return false;
  }
  //Cuts parameters
  static const double xyp0Cut_ ;
  static const double xyp2Cut_ ;
  static const double rzp0Cut_ ;
  static const double rzp1Cut_ ;
  static const double rzp2Cut_ ;
  static const double distCut_ ;
  static const double hIsoCut_ ;
  static const double dEdXSigCut_ ;
  static const double e55Cut_ ;
  static const double e55Cut2016_ ;
  static const double f51Cut_ ;
  static const double photonCut_ ;
  static const double dEdXSig_looseCut_ ;
  static const double f51_looseCut_ ;
   // if you want to set parameter in the class, you should add constexpr
  //  ex. constexpr static const double x=1;

	
private:
  string trigName_;

  // No Cut plot
  vector<PlotSet> NoCutPlot;
  vector<PlotSet> NoCutProfile;

  //N-1 cut plot box
  vector<PlotSet> n_1Plot;

  //cutflow plot box
  vector<PlotSet> CutFlow;
  vector<PlotSet> Profile;
  static const unsigned nCut = 5U;
  string cutName_[nCut];
  


  //cuts analysis
  bool evalQuality(MonoCandidate &mc) { return TMath::Abs(mc.xyp0_) < xyp0Cut_ && TMath::Abs(mc.xyp2_) > xyp2Cut_ 
			&& mc.dist_ < distCut_  && mc.hIso_ <hIsoCut_  
			&& TMath::Abs( mc.rzp2_) < rzp2Cut_ && TMath::Abs(mc.rzp1_) < rzp1Cut_ && TMath::Abs(mc.rzp0_) < rzp0Cut_;  }
  bool evalE(MonoCandidate &mc) { return mc.e55_ > e55Cut_; }
  bool evalE_16(MonoCandidate &mc) { return mc.e55_ > e55Cut2016_; }
  bool evalF51(MonoCandidate &mc) { return mc.f51_ > f51Cut_ ; }
  bool evaldEdX(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSigCut_ ;}
  bool evalPhoton(Photon &mc){ return mc.pho_pt_ > photonCut_; }
 
  bool evalPreselection(MonoCandidate &mc) { return TMath::Abs(mc.xyp0_) < xyp0Cut_ && TMath::Abs(mc.xyp2_) > xyp2Cut_ 
			&& mc.dist_ < distCut_  && mc.hIso_ <hIsoCut_  
			&& TMath::Abs( mc.rzp2_) < rzp2Cut_ && TMath::Abs(mc.rzp1_) < rzp1Cut_ && TMath::Abs(mc.rzp0_) < rzp0Cut_
			&& mc.e55_ > e55Cut_ && mc.dEdXSig_ > dEdXSigCut_;  }
  bool evalF51loose(MonoCandidate &mc) { return mc.f51_ > f51_looseCut_ ; }
  bool evaldEdXloose(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSig_looseCut_ ;}

  vector<MonoCandidate> CutFlowCand_TRG;
  vector<MonoCandidate> CutFlowCand_Qual;
  vector<MonoCandidate> CutFlowCand_Energy;
  vector<MonoCandidate> CutFlowCand_F51;
  vector<MonoCandidate> CutFlowCand_Dedx;
  vector<MonoCandidate> CutFlowCand_looseF51;
  vector<MonoCandidate> CutFlowCand_looseDedx;
  vector<MonoCandidate> Matched;

  //for cross talk
  vector<MonoCandidate> CutFlowCand_onlyDedx;
  vector<MonoCandidate> CutFlowCand_QualDedx;


  vector<MonoCandidate> N1CutCand_TRG;
  vector<MonoCandidate> N1CutCand_Qual;
  vector<MonoCandidate> N1CutCand_Energy;
  vector<MonoCandidate> N1CutCand_F51;
  vector<MonoCandidate> N1CutCand_Dedx;
	
  vector<Photon> HighPtPhoton;

  //no. of every selections 
  int Qual=0;
  int E=0;
  int f51=0;
  int dEdX=0;

  //count for cutflow 
  int count=0;
  int Qual_count=0;
  int E_count=0;
  int f51_count=0;
  int dEdX_count=0; 
  int dEdX_loose_count=0; 
  int f51_loose_count=0;
  
  //relative eff without HLT
  int NoTRG=0;
  int NoQual=0;
  int NoE=0;
  int NoF51=0;
  int NodEdXCut=0;
  //for cross talk
  int dEdXOnly_count=0;
  int QualdEdX_count=0;

  // for trigger study
  int photonLike=0;

};

#endif
