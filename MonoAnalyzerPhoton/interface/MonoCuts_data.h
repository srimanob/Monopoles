#ifndef _MONOCUTS_DATA_H_
#define _MONOCUTS_DATA_H_
class MonoCuts:public MonoCandidate
{
public:
  MonoCuts(){}

  MonoCuts(const string &trName,TFile *openFile):trigName_(trName){
	
	//create a new directory in the output file
	openFile->mkdir(trName.c_str());
	 
	cutName_[0] = "Quality_";
	 cutName_[1] = "Energy_";
	 cutName_[2] = "F51_";
	 cutName_[3] = "dEdxSig_";
	 cutName_[4] = "HLT_";
	
	// nCut-1 is for ignoring HLT_ ,omit extra plots	
	CutFlow.resize(100);
	PlotSet &y = CutFlow[0];
	string cutflowName = "";

        y.CreatPlot(FracSatVNstrips,new TH2D((cutflowName+"FracSatVNstrips").c_str(),"",30,0,1000,100,0,1));
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
        y.CreatPlot(ABCD,new TH2D((cutflowName+"ABCD").c_str(),"",105,0,1.05,60,0,30));

        //y.CreatPlot(DedXSig,new TH1D((cutflowName+"DedXSig").c_str(),"",30,0,30));

        //y.CreatPlot(F51_slice1,new TH1D((cutflowName+"F51_slice1").c_str(),"",30,0.,30));
        //y.CreatPlot(F51_slice2,new TH1D((cutflowName+"F51_slice2").c_str(),"",30,0.,30));
        //y.CreatPlot(F51_slice3,new TH1D((cutflowName+"F51_slice3").c_str(),"",30,0.,30));
        //y.CreatPlot(F51_slice4,new TH1D((cutflowName+"F51_slice4").c_str(),"",30,0.,30));
        //y.CreatPlot(F51_slice5,new TH1D((cutflowName+"F51_slice5").c_str(),"",30,0.,30));
        //y.CreatPlot(F51_slice6,new TH1D((cutflowName+"F51_slice6").c_str(),"",30,0.,30));

        //y.CreatPlot(F51,new TH1D((cutflowName+"F51").c_str(),"",30,0.,1.1));

        //y.CreatPlot(dEdXSig_slice1,new TH1D((cutflowName+"DedXSig_slice1").c_str(),"",30,0,1.1));
        //y.CreatPlot(dEdXSig_slice2,new TH1D((cutflowName+"DedXSig_slice2").c_str(),"",30,0,1.1));
        //y.CreatPlot(dEdXSig_slice3,new TH1D((cutflowName+"DedXSig_slice3").c_str(),"",30,0,1.1));
        //y.CreatPlot(dEdXSig_slice4,new TH1D((cutflowName+"DedXSig_slice4").c_str(),"",30,0,1.1));
        //y.CreatPlot(dEdXSig_slice5,new TH1D((cutflowName+"DedXSig_slice5").c_str(),"",30,0,1.1));
        //y.CreatPlot(dEdXSig_slice6,new TH1D((cutflowName+"DedXSig_slice6").c_str(),"",30,0,1.1));
        //y.CreatPlot(dEdXSig_slice7,new TH1D((cutflowName+"DedXSig_slice7").c_str(),"",30,0,1.1));

	Profile.resize(5);	
	PlotSet &x = Profile[0];
	x.CreatProfile(EcalBarrel ,new TProfile((cutflowName+"EcalBarrel").c_str(),"",30,0,1.1));
	x.CreatProfile(EcalEndCup ,new TProfile((cutflowName+"EcalEndCup").c_str(),"",30,0,1.1));
	x.CreatProfile(EcalAll ,new TProfile((cutflowName+"EcalAll").c_str(),"",30,0,1.1));
	//cout<<"end of create plots"<<endl;
  }

  ~MonoCuts(){}
  void doAnalysis_data(vector<MonoCandidate> &cand,unsigned nParticle,bool passHLT_Photon200,unsigned ev);
  void FillNoCutHistogram(int n,vector<MonoCandidate> Cand);
  void FillFlowHistogram(int n, vector<MonoCandidate> CutFlowCand);
  void Clear();
  void WritePlots(TFile *oFile);
  void SignalEff(const string trName, double TotalEvents);

  void SaveAs_csv(const string fileName, double TotalEvents){
        ofstream fout(fileName);

        if(!fout) cout<<"Cannot open file"<<endl;
        else    cout<<"Open file successfully"<<endl;
	fout<<",CutFlow,"<<endl;
        fout<<"Generated ev,"<<TotalEvents<<endl;
        fout<<"        TRG, "<<count<<endl;
        fout<<"QualityCuts, "<<Qual_count<<endl;
        fout<<"       ECut, "<<E_count<<endl;
        fout<<"     F51Cut, "<<f51_count<<endl;
        fout<<" dEdXSigCut, "<<dEdX_count<<endl;
        fout<<"Signal efficiency, "<<(double)dEdX_count/(double)TotalEvents<<endl;
        fout<<endl;
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
  static const double dEdXSig_looseCut_ ;
  static const double e55Cut_ ;
  static const double f51Cut_ ;
  static const double f51_looseCut_ ;
  static const double photonCut_ ;
   // if you want to set parameter in the class, you should add constexpr
  //  ex. constexpr static const double x=1;

	
private:
  string trigName_;
  // No Cut plot
  vector<PlotSet> NoCutPlot;
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
  bool evalF51(MonoCandidate &mc) { return mc.f51_ > f51Cut_ ; }
  bool evalF51loose(MonoCandidate &mc) { return mc.f51_ > f51_looseCut_ ; }
  bool evaldEdX(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSigCut_ ;}
  bool evaldEdXloose(MonoCandidate &mc) { return mc.dEdXSig_ > dEdXSig_looseCut_ ;}
  //---preselection for data
  bool evalPreselection(MonoCandidate &mc) { return TMath::Abs(mc.xyp0_) < xyp0Cut_ && TMath::Abs(mc.xyp2_) > xyp2Cut_ 
			&& mc.dist_ < distCut_  && mc.hIso_ <hIsoCut_  
			&& TMath::Abs( mc.rzp2_) < rzp2Cut_ && TMath::Abs(mc.rzp1_) < rzp1Cut_ && TMath::Abs(mc.rzp0_) < rzp0Cut_
			&& mc.e55_ > e55Cut_ ;  }

  //for ABCD region
  //
  //	7 8 9
  //	4 5 6
  //	1 2 3
  //
  bool evalRegion1(MonoCandidate &mc) { return  0 < mc.dEdXSig_  && mc.dEdXSig_ < dEdXSig_looseCut_ &&
						0 < mc.f51_ 	 && mc.f51_	< f51_looseCut_;}
  bool evalRegion2(MonoCandidate &mc) { return  0 < mc.dEdXSig_  && mc.dEdXSig_ < dEdXSig_looseCut_ &&
						f51_looseCut_ < mc.f51_ && mc.f51_ < f51Cut_;}
  bool evalRegion3(MonoCandidate &mc) { return  0 < mc.dEdXSig_  && mc.dEdXSig_ < dEdXSig_looseCut_ &&
								 mc.f51_	> f51Cut_;}
  bool evalRegion4(MonoCandidate &mc) { return  dEdXSig_looseCut_< mc.dEdXSig_  && mc.dEdXSig_ < dEdXSigCut_ &&
						0 < mc.f51_ 	 && mc.f51_	< f51_looseCut_;}
  bool evalRegion5(MonoCandidate &mc) { return  dEdXSig_looseCut_< mc.dEdXSig_  && mc.dEdXSig_ < dEdXSigCut_ && 
						f51_looseCut_ < mc.f51_ && mc.f51_ < f51Cut_;}
  bool evalRegion6(MonoCandidate &mc) { return  dEdXSig_looseCut_< mc.dEdXSig_  && mc.dEdXSig_ < dEdXSigCut_ &&
								  mc.f51_	> f51Cut_;}
  bool evalRegion7(MonoCandidate &mc) { return  mc.dEdXSig_  > dEdXSigCut_ &&
						0 < mc.f51_ 	 && mc.f51_	< f51_looseCut_;}
  bool evalRegion8(MonoCandidate &mc) { return   mc.dEdXSig_  > dEdXSigCut_ &&
						f51_looseCut_ < mc.f51_ && mc.f51_ < f51Cut_;}
  bool evalRegion9(MonoCandidate &mc) { return   mc.dEdXSig_  > dEdXSigCut_ &&
								mc.f51_ > f51Cut_;}

 
  vector<MonoCandidate> CutFlowCand_TRG;
  vector<MonoCandidate> CutFlowCand_Qual;
  vector<MonoCandidate> CutFlowCand_Energy;
  vector<MonoCandidate> CutFlowCand_F51;
  vector<MonoCandidate> CutFlowCand_Dedx;
  //for data
  vector<MonoCandidate> Preselection;
  vector<MonoCandidate> TriggerSelection;
  vector<MonoCandidate> selectedCand;
  vector<MonoCandidate> otherCand;
  vector<MonoCandidate> Region1;
  vector<MonoCandidate> Region2;
  vector<MonoCandidate> Region3;
  vector<MonoCandidate> Region4;
  vector<MonoCandidate> Region5;
  vector<MonoCandidate> Region6;
  vector<MonoCandidate> Region7;
  vector<MonoCandidate> Region8;
  vector<MonoCandidate> Region9;
	
  //ABCD region counting
  int region1=0;
  int region2=0;
  int region3=0;
  int region4=0;
  int region5=0;
  int region6=0;
  int region7=0;
  int region8=0;
  int region9=0;

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
  int SpikeLike=0;
  int EBarrel=0; 
  int Reco=0;
  int photonLike=0;
  //relative eff without HLT
  int NoTRG=0;
  int NoQual=0;
  int NoE=0;
  int NoF51=0;
  int NodEdXCut=0;

  int MatchedEvent=0;
  bool Flag=0;
};

#endif
