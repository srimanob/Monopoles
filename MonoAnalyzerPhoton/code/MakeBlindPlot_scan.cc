//////////////////////////////////
//
//	MakeBlindPlot_scan.cc
//	Shih Lin March 4 2022
//	
//	root -l 'MakeBlindPlot_scan.cc(Unblind*,"year**")'
// 	*Unblind(int) 
// 	 0 to show 5 regions
// 	 1 to show 8, 
// 	 2 to show full signal, 
// 	 3 for MC
// 	**year (string)
// 	 16 : 2016 + 2016 APV
// 	 1718 : 2017 + 2018
//	
/////////////////////////////////
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

float Run, Event, SatSubHits, SubHits, Dist, HIso, XYPar0, XYPar1,XYPar2, RZPar0, RZPar1, RZPar2,Eta, seedFrac, cand_e55, passHLT_Photon200;
void ChiSquare_scan(double *x_,double *y_, int nloop, vector<double> ChiSquare ,string year);
TChain *SetMonoAddress(TChain *Mono,string year){

    //Mono->Add(("../../Data/data_"+year+"/*").c_str());
    if(year == "1718"){
	    Mono->Add(("../../Data/data_"+year+"/data_2017/*").c_str());
	    Mono->Add(("../../Data/data_"+year+"/data_2018/*").c_str());
    }
    else{
	    Mono->Add(("../../Data/data_"+year+"/data_2016/*").c_str());
	    Mono->Add(("../../Data/data_"+year+"/data_2016APV/*").c_str());
    }	

    Mono->SetBranchAddress("run", &Run);
    Mono->SetBranchAddress("event", &Event);// every Event
    Mono->SetBranchAddress("Dist", &Dist);
    Mono->SetBranchAddress("HIso", &HIso);
    Mono->SetBranchAddress("XYPar0", &XYPar0);
    Mono->SetBranchAddress("XYPar1", &XYPar1);
    Mono->SetBranchAddress("XYPar2", &XYPar2);
    Mono->SetBranchAddress("RZPar0", &RZPar0);
    Mono->SetBranchAddress("RZPar1", &RZPar1);
    Mono->SetBranchAddress("RZPar2", &RZPar2);
    Mono->SetBranchAddress("Eta", &Eta);
    Mono->SetBranchAddress("SatSubHits", &SatSubHits);
    Mono->SetBranchAddress("SubHits", &SubHits);
    Mono->SetBranchAddress("seedFrac", &seedFrac);
    Mono->SetBranchAddress("cand_e55", &cand_e55);
    Mono->SetBranchAddress("passHLT_Photon200", &passHLT_Photon200);
    return Mono;
}

void doAnalysis(TChain *Mono,TH2F *Plot,string year){

    int LastEvent = -1;
    vector<pair<float,float>> Point;
    for(int ev=0; ev<Mono->GetEntries(); ev++){
	Mono->GetEntry(ev);
	if (ev%10000000==0) cout<<(float)ev<<"/"<<Mono->GetEntries()<<endl;	
	
	if(year == "16"){
	      if(!(
	              Dist < 0.5  
	            && HIso < 10 
	            && abs(XYPar0) < 0.6 
	            && abs(XYPar1) < 10 
	            && abs(XYPar2) > 1000 
	            && abs(RZPar0) < 10 
	            && abs(RZPar1) < 999 
	            && abs(RZPar2) < 0.005 
	            && cand_e55 > 175
//	            &&  sqrt(-TMath::Log(TMath::BinomialI(0.07, SubHits, SatSubHits))) >2 
	             ) ) continue;
	}
	else if(year == "1718"){
	      if(!(
	  	    passHLT_Photon200 == 1  &&
	              Dist < 0.5  
	            && HIso < 10 
	            && abs(XYPar0) < 0.6 
	            && abs(XYPar1) < 10 
	            && abs(XYPar2) > 1000 
	            && abs(RZPar0) < 10 
	            && abs(RZPar1) < 999 
	            && abs(RZPar2) < 0.005 
	            && cand_e55 > 200
//	            &&  sqrt(-TMath::Log(TMath::BinomialI(0.07, SubHits, SatSubHits))) >2 
	             ) ) continue;

	}
	if((int)Event!=LastEvent && Point.size()>0){
	//1. Only see different events to do the ABCD method
	   sort(Point.begin(),Point.end());

	   for(int i=Point.size()-1;i>=0;i--){
		Plot->Fill(Point[i].second,Point[i].first);
		//2nd->f51 for x-axis, 1st->dedx for y-axis)
	  	if(Point[i].second>0.0 && Point[i].first > 0.0 ) break;
		// only fill one ev at a point, not superposition	
	    }	
		Point.clear();
		LastEvent = Event;
		
	}

	float Significance = sqrt(-TMath::Log(TMath::BinomialI(0.07, SubHits, SatSubHits))) ;
        if (Significance > 29) Significance = 29;

        Point.push_back(pair<float,float>(Significance, seedFrac));
   }

}
void MakeBlindPlot_scan(int Unblind, string year){

    TCanvas *c = new TCanvas("c","",800,600);
    TChain *Mono  = new TChain("monopoles");
    
    SetMonoAddress(Mono,year);
	
	
    double x_[100] = {0};
    double y_[100] = {0};
    int i = 0;
    vector<double> eChiSquare;
    vector<double> ChiSquare;

    double x_tight = 85;
    double y_tight = 18;

    TH2F *Plot = new TH2F("Plot", "Data", 105, 0, 1.05, 60, 0, 30);
    Plot->SetXTitle("Frac51");
    Plot->SetYTitle("HighDeDx Significance");
    doAnalysis(Mono,Plot,year);

 for(int x_loose = 40 ; x_loose < x_tight; ){	//f51 40, 45, 50, 55, 60, 65, 70, 75, 80
    for(int y_loose = 8; y_loose < y_tight;){  //dedx 8,10,12,14,16 (4, 5, 6, 7,8)
	
    TH2F *newPlot = (TH2F*)Plot->Clone("newPlot");
    string LooseX = to_string(x_loose);
    string LooseY = to_string(y_loose);
    newPlot->SetName(("Plot_"+LooseX+"_"+LooseY).c_str());

    double x[4]={0,(double)x_loose, x_tight, 105};//0.6 0.85
    double y[4]={0,(double)y_loose, y_tight, 60};//7 9

    double xbins[4] = {x[0]/100, x[1]/100, x[2]/100, x[3]/100};//0, 0.60, 0.85,1.05 x-axis
    double ybins[4] = {y[0]/2, y[1]/2, y[2]/2, y[3]/2};//0, 7, 9, 30 y-axis
    
 
    TH2F *Actual = new TH2F("Actual","",3,xbins,3,ybins);
    TH2F *Expected = new TH2F("Expected","",3,xbins,3,ybins);


    int LowEdgeX, LowEdgeY;
    if(Unblind>=2){
	LowEdgeX = 999;
	LowEdgeY = 999;
    }
    else if(Unblind==1){
	LowEdgeX = x[2];//0.85
	LowEdgeY = y[2];//9
    }
    else{
	LowEdgeX = x[1];//0.60
	LowEdgeY = y[1];//7
    }

    for(int i=LowEdgeX+1; i<=x[3]+1; i++){
        for(int j=LowEdgeY+1; j<=y[3]+1; j++){
            newPlot->SetBinContent(i, j, 0);
        }
    }


    float a[10];//a[0] nothing

    a[7] = newPlot->Integral(x[0]+1,x[1],y[2]+1,y[3]);
    a[4] = newPlot->Integral(x[0]+1,x[1],y[1]+1,y[2]);
    a[1] = newPlot->Integral(x[0]+1,x[1],y[0]+1,y[1]);
    a[2] = newPlot->Integral(x[1]+1,x[2],y[0]+1,y[1]);
    a[3] = newPlot->Integral(x[2]+1,x[3],y[0]+1,y[1]);
    a[5] = newPlot->Integral(x[1]+1,x[2],y[1]+1,y[2]);
    a[8] = newPlot->Integral(x[1]+1,x[2],y[2]+1,y[3]);
    a[6] = newPlot->Integral(x[2]+1,x[3],y[1]+1,y[2]);
    a[9] = newPlot->Integral(x[2]+1,x[3],y[2]+1,y[3]);

    // fill actual data inti plots 
    for(int i=1;i<=9;i++){
	Actual->SetBinContent( ((i-1)%3)+1, (i-1)/3+1, a[i]);
    }

    //calculate expected value in signal region( region5,6,8,9)
    //
    //	7  8  9
    //  4  5  6
    //  1  2  3
    //
    //  SetBinContent(x,y,Content);
    //
    //region 8
    Expected->SetBinContent(2,3,((a[7]+1)*(a[2]+1)/(a[1]+1)));
    Expected->SetBinError(2,3, Expected ->GetBinContent(2,3)*sqrt(1/(a[7]+1)+1/(a[2]+1)+1/(a[1]+1)));
    //region 5
    Expected->SetBinContent(2, 2, (a[2]+1)*(a[4]+1)/(a[1]+1));
    Expected->SetBinError(2, 2, Expected->GetBinContent(2,2)*sqrt(1/(a[2]+1) + 1/(a[4]+1) + 1/(a[1]+1)) );
    //region 6
    Expected->SetBinContent(3, 2, (a[3]+1)*(a[4]+1)/(a[1]+1));
    Expected->SetBinError(3, 2, Expected->GetBinContent(3,2)*sqrt(1/(a[3]+1) + 1/(a[4]+1) + 1/(a[1]+1)) );
    //region 9 
    Expected->SetBinContent(3, 3, (a[3]+1)*(a[7]+1)/(a[1]+1));
    Expected->SetBinError(3, 3, Expected->GetBinContent(3,3)*sqrt(1/(a[3]+1) + 1/(a[7]+1) + 1/(a[1]+1)) );
    if(Unblind > 0){
        Expected->SetBinContent(3, 3, (a[3]+a[6]+1)*(a[7]+a[8]+1)/(a[1]+a[2]+a[4]+a[5]+1));
        Expected->SetBinError(3, 3, Expected->GetBinContent(3,3)*sqrt(1/(a[3]+a[6]+1) + 1/(a[7]+a[8]+1) + 1/(a[1]+a[2]+a[4]+a[5]+1)) );
    }

    cout<<"f51  "<<(double)x_loose/100<<" , dEdxSig  "<<(double)y_loose/2<<endl;
    cout << "\t\tExpect\terr\tActual\terr" << endl;

    cout << "Region 7:\t\t\t" << Actual->GetBinContent(1,3)<< '\t'<<Actual->GetBinError(1,3)<<endl;
    cout << "Region 4:\t\t\t" << Actual->GetBinContent(1,2)<< '\t'<<Actual->GetBinError(1,2)<<endl;
    cout << "Region 1:\t\t\t" << Actual->GetBinContent(1,1)<< '\t'<<Actual->GetBinError(1,1)<<endl;
    cout << "Region 2:\t\t\t" << Actual->GetBinContent(2,1)<< '\t'<<Actual->GetBinError(2,1)<<endl;
    cout << "Region 3:\t\t\t" << Actual->GetBinContent(3,1)<< '\t'<<Actual->GetBinError(3,1)<<endl;

    cout << "Region 8:\t" << Expected->GetBinContent(2,3) << '\t' << Expected->GetBinError(2,3)<<'\t'<<Actual->GetBinContent(2,3)<<'\t'<<Actual->GetBinError(2,3)<< endl;
    cout << "Region 5:\t" << Expected->GetBinContent(2,2) << '\t' << Expected->GetBinError(2,2)<<'\t'<<Actual->GetBinContent(2,2)<<'\t'<<Actual->GetBinError(2,2)<< endl;
    cout << "Region 6:\t" << Expected->GetBinContent(3,2) << '\t' << Expected->GetBinError(3,2)<<'\t'<<Actual->GetBinContent(3,2)<<'\t'<<Actual->GetBinError(3,2)<< endl;

    if(Unblind==0)
        cout << "Region 9:\t" << a[3]*a[7]/a[1] << '\t' << a[9] << endl;
    else{
    
	cout << "Region 9:\t" << Expected->GetBinContent(3,3) << '\t' << Expected->GetBinError(3,3)<<'\t'<<Actual->GetBinContent(3,3)<<'\t'<<Actual->GetBinError(3,3)<< endl;
//        Expecteded->SetBinContent(3, 3, (a[3]+a[6])*(a[7]+a[8])/(a[1]+a[2]+a[4]+a[5]));//????????
    }

    newPlot->SetStats(0);
    newPlot->SetLineColor(51);
    newPlot->SetFillColor(51);
    Actual->SetMarkerColor(kBlack);//the color of words for actual
    Actual->SetLineColor(kBlack);
    Actual->SetMarkerSize(1.5);
    Expected->SetMarkerColor(kRed);
    Expected->SetLineColor(kRed);
    Expected->SetMarkerSize(1.5);

    newPlot->Draw("col");
    Actual->Draw("same text");
    Expected->Draw("same text e");

    TLegend *leg = new TLegend(0.15,0.67,0.36,0.88);
    leg->AddEntry(Actual,"Actual count","l");
    leg->AddEntry(Expected,"Expected count","l");
    leg->Draw();

    TLine *l1=new TLine(xbins[0], ybins[1], xbins[3], ybins[1]);
    TLine *l2=new TLine(xbins[0], ybins[2], xbins[3], ybins[2]);
    TLine *l3=new TLine(xbins[1], ybins[0], xbins[1], ybins[3]);
    TLine *l4=new TLine(xbins[2], ybins[0], xbins[2], ybins[3]);

    l1->SetLineStyle(2); l2->SetLineStyle(2); l3->SetLineStyle(2); l4->SetLineStyle(2);
    l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();

    ChiSquare.push_back(pow(Actual->GetBinContent(2,2)-Expected->GetBinContent(2,2),2)/(pow(Actual->GetBinError(2,2),2)+pow(Expected->GetBinError(2,2),2))
				+pow(Actual->GetBinContent(2,3)-Expected->GetBinContent(2,3),2)/(pow(Actual->GetBinError(2,3),2)+pow(Expected->GetBinError(2,3),2))
				+pow(Actual->GetBinContent(3,2)-Expected->GetBinContent(3,2),2)/(pow(Actual->GetBinError(3,2),2)+pow(Expected->GetBinError(3,2),2)));

	    ofstream fout("output/csv_file/ABCD_"+year+"_"+LooseX+"_"+LooseY+".csv"); 
	    fout<<year<<",f51,dEdXSig"<<endl;
	    fout<<","<<(double)x_loose/100<<","<<(double)y_loose/2<<endl;
	    fout<<"Region,actual,error,expect,error,uncertainty"<<endl;
	    fout<<"1,"<<Actual->GetBinContent(1,1)<<","<<Actual->GetBinError(1,1)<<endl;
	    fout<<"2,"<<Actual->GetBinContent(2,1)<<","<<Actual->GetBinError(2,1)<<endl;
	    fout<<"3,"<<Actual->GetBinContent(3,1)<<","<<Actual->GetBinError(3,1)<<endl;
	    fout<<"4,"<<Actual->GetBinContent(1,2)<<","<<Actual->GetBinError(1,2)<<endl;
	    fout<<"7,"<<Actual->GetBinContent(1,3)<<","<<Actual->GetBinError(1,3)<<endl;
	    fout<<"5,"<<Actual->GetBinContent(2,2)<<","<<Actual->GetBinError(2,2)<<","<<Expected->GetBinContent(2,2)<<","<<Expected->GetBinError(2,2)<<","<<Expected->GetBinError(2,2)/Expected->GetBinContent(2,2)<<endl;
	    fout<<"6,"<<Actual->GetBinContent(3,2)<<","<<Actual->GetBinError(3,2)<<","<<Expected->GetBinContent(3,2)<<","<<Expected->GetBinError(3,2)<<","<<Expected->GetBinError(3,2)/Expected->GetBinContent(3,2)<<endl;
	    fout<<"8,"<<Actual->GetBinContent(2,3)<<","<<Actual->GetBinError(2,3)<<","<<Expected->GetBinContent(2,3)<<","<<Expected->GetBinError(2,3)<<","<<Expected->GetBinError(2,3)/Expected->GetBinContent(2,3)<<endl;
	    fout<<"9,"<<Actual->GetBinContent(3,3)<<","<<Actual->GetBinError(3,3)<<","<<Expected->GetBinContent(3,3)<<","<<Expected->GetBinError(3,3)<<","<<Expected->GetBinError(3,3)/Expected->GetBinContent(3,3)<<endl;
	    fout.close();

    x_[i] = x_loose;
    y_[i] = y_loose;
    i++;

    y_loose = y_loose+1;
    }
    x_loose = x_loose+5;
  } 
	int nloop = 0;
	nloop = i;	
	cout<<"Loop "<<nloop<<endl;
	ChiSquare_scan(x_,y_,nloop,ChiSquare,year);
}
void ChiSquare_scan(double *x_,double *y_, int nloop, vector<double> ChiSquare,string year){

        ofstream Scan_result(("output/csv_file/ABCD_scan_"+year+".csv").c_str()); 
	Scan_result<<"f51,dEdxSig,Chi-square"<<endl;
	for(int i =0;i<nloop;i++){
		Scan_result<<x_[i]<<","<<y_[i]<<","<<ChiSquare[i]<<endl;
	}

}
	
