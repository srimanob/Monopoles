//////////////////////////////////
//
//	MakeBlindPlot.cc
//	Shih Lin
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
    vector<double> key;

    cout<<"NEvents "<<Mono->GetEntries()<<endl;
    for(int ev=0; ev<Mono->GetEntries(); ev++){
	//one event is one track/candidate(MC)
	Mono->GetEntry(ev);
	if (ev%10000000==0) cout<<(float)ev<<"/"<<Mono->GetEntries()<<endl;	
		// you need to apply these preselection and trigger cut
		// so that you can "estimate" the background !!( those cutted entries are not relate with our bkg!!)
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
	//note that same event can also push back into the point.
	//so the point is a vector array(like candidates in one event but
	//pass the preselection)
	//After enter the other event(so we enter the Event!=LastEvent),
	//this point(candidate) is sorted. and then fill into the plot. 
	//Besides, since we need to get new candidate info from new even
	//so we "clear" the Point, then put new candidate into Point.
   }
}


void MakeBlindPlot(int Unblind, string year){

    TCanvas *c = new TCanvas("c","",800,600);
    TChain *Mono  = new TChain("monopoles");
    
    SetMonoAddress(Mono,year);

    double f51_loose = 0; //0.6
    double dEdx_loose = 0; // 6.5
	
    if( year == "16"){
	
	f51_loose = 60; //0.6
	dEdx_loose = 13; // 6.5
    }
    else if( year == "1718"){

	f51_loose = 75; // 0.75
	dEdx_loose = 14; // 7
	
    }
	
	cout<<f51_loose<<" "<<dEdx_loose<<endl;

    double x[4]={0,f51_loose, 85, 105};
    double y[4]={0,dEdx_loose, 18, 60};
    double VetoX = 0.0;
    double VetoY = 0.0;

    double xbins[4] = {x[0]/100, x[1]/100, x[2]/100, x[3]/100};//0, 0.60, 0.85,1.05 x-axis
    double ybins[4] = {y[0]/2, y[1]/2, y[2]/2, y[3]/2};//0, 7, 9, 30 y-axis
    
 
    TH2F *Plot = new TH2F("Plot", "Data", 105, 0, 1.05, 60, 0, 30);
    Plot->SetXTitle("Frac51");
    Plot->SetYTitle("HighDeDx Significance");

    doAnalysis(Mono,Plot,year);

    TH2F *Actual = new TH2F("Actual","",3,xbins,3,ybins);
    TH2F *Expected = new TH2F("Expected","",3,xbins,3,ybins);

    int LowEdgeX, LowEdgeY;
    //Unblind is 0 to show 5 regions, 1 to show 8 regions, 2 to show full signal, 3 for MC
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
            Plot->SetBinContent(i, j, 0);//to let the Plot have nothing. The data will be show on the Actual and Expected.
        }
    }


    float a[10];//a[0] nothing

    a[7] = Plot->Integral(x[0]+1,x[1],y[2]+1,y[3]);
    a[4] = Plot->Integral(x[0]+1,x[1],y[1]+1,y[2]);
    a[1] = Plot->Integral(x[0]+1,x[1],y[0]+1,y[1]);
    a[2] = Plot->Integral(x[1]+1,x[2],y[0]+1,y[1]);
    a[3] = Plot->Integral(x[2]+1,x[3],y[0]+1,y[1]);
    a[5] = Plot->Integral(x[1]+1,x[2],y[1]+1,y[2]);
    a[8] = Plot->Integral(x[1]+1,x[2],y[2]+1,y[3]);
    a[6] = Plot->Integral(x[2]+1,x[3],y[1]+1,y[2]);
    a[9] = Plot->Integral(x[2]+1,x[3],y[2]+1,y[3]);

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

    cout << "\t\tExpect\tActual" << endl;

    cout << "Region 7:\t\t" << a[7] << endl;
    cout << "Region 4:\t\t" << a[4] << endl;
    cout << "Region 1:\t\t" << a[1] << endl;
    cout << "Region 2:\t\t" << a[2] << endl;
    cout << "Region 3:\t\t" << a[3] << endl;

    cout << "Region 8:\t" << a[2]*a[7]/a[1] << '\t' << a[8] << endl;
    cout << "Region 5:\t" << a[2]*a[4]/a[1] << '\t' << a[5] << endl;
    cout << "Region 6:\t" << a[3]*a[4]/a[1] << '\t' << a[6] << endl;

    if(Unblind==0)
        cout << "Region 9:\t" << a[3]*a[7]/a[1] << '\t' << a[9] << endl;
    else{
        cout << "Region 9:\t" << (a[3]+a[6])*(a[7]+a[8])/(a[1]+a[2]+a[4]+a[5]) << '\t' << a[9] << endl;
        Expected->SetBinContent(3, 3, (a[3]+a[6])*(a[7]+a[8])/(a[1]+a[2]+a[4]+a[5]));
    }



    Plot->SetStats(0);
    Plot->SetLineColor(51);
    Plot->SetFillColor(51);
    Actual->SetMarkerColor(kBlack);
    Actual->SetLineColor(kBlack);
    Actual->SetMarkerSize(1.5);
    Expected->SetMarkerColor(kRed);
    Expected->SetLineColor(kRed);
    Expected->SetMarkerSize(1.5);

    ofstream fout("output/csv_file/BKG_"+year+".csv"); 
    fout<<year<<", bkg uncertainty"<<endl;
    fout<<"exp,"<<Expected->GetBinContent(3,3)<<endl;
    fout<<"err,"<<Expected->GetBinError(3,3)<<endl;
    fout<<"unc,"<<Expected->GetBinError(3,3)/Expected->GetBinContent(3,3)<<endl;
    fout.close();

    Plot->Draw("col");
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

    c->SaveAs(("ABCD_plot"+year+".pdf").c_str()); 

}	
