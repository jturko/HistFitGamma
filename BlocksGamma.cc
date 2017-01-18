
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TColor.h"
#include "TApplication.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom3.h"

#include "TMath.h"

#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

class vec 
{
public:
    vec() { 
        set(0,0,0); 
    }
    ~vec() {} 
    
    vec(double e1, double e2, double e3) {
        set(e1,e2,e3);
    }

    double at(int i) { 
        if(i<0||i>=int(fVector.size())) { 
            std::cout << "out of range!" << std::endl;
            return 0.;
        }
        else return fVector.at(i); 
    }

    int size() { return fVector.size(); }
 
    void set(double e1, double e2, double e3) {
        fVector.resize(0);
        fVector.push_back(e1); fArray[0] = e1;
        fVector.push_back(e2); fArray[1] = e2;
        fVector.push_back(e3); fArray[2] = e3;
    }
    void set(int i, double val) { 
        if(i<0||i>=int(fVector.size())) { 
            std::cout << "out of range!" << std::endl;
            return;
        }
        fVector.at(i) = val; fArray[i] = val;
    }
    
    void add(vec v) { 
        if(v.size() != int(fVector.size())) { std::cout << "error, different vector sizes" << std::endl; return; }
        for(int i=0; i<int(fVector.size()); i++) fVector.at(i) = fVector.at(i) + v.at(i);   
    }
    void subtract(vec v) { 
        if(v.size() != int(fVector.size())) { std::cout << "error, different vector sizes" << std::endl; return; }
        for(int i=0; i<int(fVector.size()); i++) fVector.at(i) = fVector.at(i) - v.at(i);   
    }
    
    vec scalar_multiply(double num) {
        vec return_vec;
        for(int i=0; i<int(fVector.size()); i++) return_vec.set(i,num*fVector.at(i));
        return return_vec;
    }
    
    vec midpoint(vec v) { 
        vec return_vec;
        if(v.size() != int(fVector.size())) {
            std::cout << "error, different vector sizes" << std::endl;
            return return_vec;
        }
        for(int i=0; i<int(fVector.size()); i++) { return_vec.set(i,0.5*(fVector.at(i) + v.at(i))); }
        return return_vec;
    }
    
    double * par_array() {
        return fArray;
    }

private:
    std::vector<double> fVector;
    double fArray[3];
};

/////////////////////////////////////////////////////////////////////////////////////

class HistFit 
{

public:
    HistFit(std::string run_id = "0_0");
    ~HistFit();

    void PrintParameters() {
        std::cout << "      A = " << fParameters[0] << std::endl;
        std::cout << "      B = " << fParameters[1] << std::endl;
        std::cout << "      C = " << fParameters[2] << std::endl;
        std::cout << " offset = " << fParameters[8] << " keVee " << std::endl;
    }

    void Sort(double * par);
    void Sort(double A=0.1, double B=0.05, double C=1e-4, double offset=0.);

    void SetParameters(double * par);

    double LightOutput(double E, double * par) {
        return ( par[0]*E-par[1]*(1.0-TMath::Exp(-par[2]*TMath::Power(E,par[3]))) );
        //return ( fLightOffset/1000. + par[0]*E-par[1]*(1.0-TMath::Exp(-par[2]*TMath::Power(E,par[3]))) );
    }
    double Resolution(double E, double * par) {
        return (E*TMath::Sqrt(TMath::Power(par[0],2)+TMath::Power(par[1],2)/E+TMath::Power(par[2]/E,2)))/(2.*TMath::Sqrt(2.*TMath::Log(2.))) ;
    }

    void ApplyCutoffLow(double cutoff, std::string str) {
        if(str == "sim")      for(int i=0; i<fSimHist->FindBin(cutoff); i++) fSimHist->SetBinContent(i,0.);
        else if(str == "exp") for(int i=0; i<fExpHist->FindBin(cutoff); i++) fExpHist->SetBinContent(i,0.);
        else                  std::cout << "--->  Error with ApplyCutoffLow!" << std::endl;
    }
    void ApplyCutoffHigh(double cutoff, std::string str) {
        if(str == "sim")      for(int i=fSimHist->FindBin(cutoff); i<fSimHist->GetNbinsX(); i++) fSimHist->SetBinContent(i,0.);
        else if(str == "exp") for(int i=fExpHist->FindBin(cutoff); i<fExpHist->GetNbinsX(); i++) fExpHist->SetBinContent(i,0.);
        else                  std::cout << "--->  Error with ApplyCutoffHigh!" << std::endl;
    }
    
    double DoChi2() {
        TH1F * h_e = (TH1F*)fExpHist->Clone();
        TH1F * h_s = (TH1F*)fSimHist->Clone();
        for(int i=h_e->FindBin(fCutoffHigh); i<h_e->GetNbinsX(); i++) h_e->SetBinContent(i,0.);
        for(int i=h_s->FindBin(fCutoffHigh); i<h_s->GetNbinsX(); i++) h_s->SetBinContent(i,0.);
        for(int i=0; i<h_e->FindBin(fCutoffLow); i++) h_e->SetBinContent(i,0.);
        for(int i=0; i<h_s->FindBin(fCutoffLow); i++) h_s->SetBinContent(i,0.);
        fChi2 = h_e->Chi2Test(h_s,"CHI2/NDF");
        delete h_e;
        delete h_s;

        return fChi2;
    }

    void Draw() {
        fExpHist->SetLineColor(kRed);
        if(fSimHist) fSimHist->SetLineColor(kBlack);
        else { std::cout << "no sim hist sorted yet!" << std::endl; return; }
        fExpHist->GetXaxis()->SetRangeUser(fCutoffLow-30,fCutoffHigh+150);
        double ymax = 2*fExpHist->GetBinContent(fExpHist->GetMaximumBin());
        fExpHist->GetYaxis()->SetRangeUser(0.1,ymax);

        fExpHist->Draw();
        fSimHist->Draw("same");   
    }

    std::string GetRunId() { return fRunId; }

    void SetCutoffHigh(double cut) { fCutoffHigh = cut; }
    void SetCutoffLow(double cut) { fCutoffLow = cut; }
    double GetCutoffHigh() { return fCutoffHigh; }    
    double GetCutoffLow() { return fCutoffLow; }    

    TH1F * GetSimHist() { return fSimHist; }
    TH1F * GetExpHist() { return fExpHist; }
    
    double GetSimEntries() { return fNumEntries; }

    std::string fRunId;

    double fCutoffHigh;
    double fCutoffLow;
    
    TFile * fSimFile;
    TFile * fExpFile;

    TH1F * fSimHist;
    TH1F * fExpHist;

    TTree * fSimTree;
    TBranch * fEdepBranch;
    TBranch * fEkinBranch;
    TBranch * fPtypeBranch;
    std::vector<double> * fEdepVector;
    std::vector<double> * fEkinVector;
    std::vector<int> * fPtypeVector;

    double fProtonCoeff[4];
    double fDeuteronCoeff[4];
    double fCarbonCoeff[4];
    double fAlphaCoeff[4];
    double fBeCoeff[4];
    double fBCoeff[4];
    double fSmearingCoeff[3];
    double fParameters[9];

    int fNumEntries;
    TRandom3 fRandom;

    double fLightOffset;

    double fChi2;
    int fExpBinNum;

};

HistFit::HistFit(std::string run_id) :
    fRunId(run_id),
    fSimFile(NULL),
    fExpFile(NULL),
    fSimHist(NULL),
    fExpHist(NULL),
    fSimTree(NULL),
    fEdepBranch(NULL),
    fEkinBranch(NULL),
    fPtypeBranch(NULL),
    fEdepVector(NULL),
    fEkinVector(NULL),
    fPtypeVector(NULL)
{

    //fCutoffHigh = 80.;
    //fCutoffLow = 3.0;
    
    std::string title;
    if(run_id=="0_2"||run_id=="3_1"||run_id=="0_3"||run_id=="3_2") { fCutoffHigh = 3000.; fCutoffLow = 100.; title="24Na";} 
    else if(run_id=="0_4"||run_id=="3_3") { fCutoffHigh = 1500; fCutoffLow = 50; title="60Co";} 
    else if(run_id=="0_0"||run_id=="1_0"||run_id=="3_0") { fCutoffHigh = 700; fCutoffLow = 70; title="137Cs";} 
    else if(run_id=="0_1"||run_id=="1_1"||run_id=="2_0"||run_id=="1_2"||run_id=="2_1") { fCutoffHigh = 80.; fCutoffLow = 3.; title="241Am";} 
    else{ std::cout << "not a valid run id, this is a garbage object!!!" << std::endl; return; }
 
    fExpFile = TFile::Open("~/data/hists2012.root"); 

    std::string hist_name = "ProtonCalibratedSource" + fRunId;
    fExpHist = (TH1F*)(fExpFile->Get(hist_name.c_str())->Clone());
    fExpHist->SetNameTitle(hist_name.c_str(),title.c_str());
    fExpHist->SetLineColor(kBlack);
    fExpHist->GetXaxis()->UnZoom(); 
    fExpHist->GetYaxis()->UnZoom();
    fExpHist->SetStats(false);
    ApplyCutoffLow(fCutoffLow,"exp");
    fExpBinNum = fExpHist->GetNbinsX();
    
    std::string name = "../G4_RAW/Sim" + fRunId + ".root";
    fSimFile = TFile::Open(name.c_str());     

    fSimTree = (TTree*)(fSimFile->Get("ntuple/ntuple")); 
    fNumEntries = fSimTree->GetEntries();
    
    fSimTree->SetBranchAddress("eDepVector",&fEdepVector,&fEdepBranch);
    fSimTree->SetBranchAddress("eKinVector",&fEkinVector,&fEkinBranch);
    fSimTree->SetBranchAddress("particleTypeVector",&fPtypeVector,&fPtypeBranch);
    
    fProtonCoeff[0] = 0.74; fProtonCoeff[1] = 3.2; fProtonCoeff[2] = 0.20; fProtonCoeff[3] = 0.97;
    fDeuteronCoeff[0] = 0.75; fDeuteronCoeff[1] = 2.80; fDeuteronCoeff[2] = 0.25; fDeuteronCoeff[3] = 0.93;
    fCarbonCoeff[0] = 0.05; fCarbonCoeff[1] = 0.0; fCarbonCoeff[2] = 0.0;fCarbonCoeff[3] = 0.0;
    fAlphaCoeff[0] = 0.14; fAlphaCoeff[1] = 0.59; fAlphaCoeff[2] = 0.065; fAlphaCoeff[3] = 1.01;
    fBeCoeff[0] = 0.0821; fBeCoeff[1] = 0.0; fBeCoeff[2] = 0.0; fBeCoeff[3] = 0.0;
    fBCoeff[0] = 0.0375; fBCoeff[1] = 0.0; fBCoeff[2] = 0.0; fBCoeff[3] = 0.0;
    
    fSmearingCoeff[0] = 0.4; fSmearingCoeff[1] = 0.07; fSmearingCoeff[2] = 0.008;

    fParameters[0] = 0.1;
    fParameters[1] = 0.05;
    fParameters[2] = 1e-4;
    fParameters[3] = 0.639;
    fParameters[4] = 1.462;
    fParameters[5] = 0.373;
    fParameters[6] = 0.968;
    fParameters[7] = 0.0;
    fParameters[8] = 0.0;
    SetParameters(fParameters);
   
    fLightOffset = 0.0;

    std::cout << "Run# = " << fRunId << std::endl;

    //if(fExpBinNum == 50100) fExpHist->Rebin(10);
}

HistFit::~HistFit() {}

void HistFit::SetParameters(double * par)
{
    fSmearingCoeff[0] = par[0];
    fSmearingCoeff[1] = par[1];
    fSmearingCoeff[2] = par[2];
    fProtonCoeff[0] = par[3];
    fProtonCoeff[1] = par[4];
    fProtonCoeff[2] = par[5];
    fProtonCoeff[3] = par[6];
    fCarbonCoeff[0] = par[7];
    fLightOffset = par[8];   
 
    fParameters[0] = par[0];
    fParameters[1] = par[1];
    fParameters[2] = par[2];
    fParameters[3] = par[3];
    fParameters[4] = par[4];
    fParameters[5] = par[5];
    fParameters[6] = par[6];
    fParameters[7] = par[7];
    fParameters[8] = par[8];
}

void HistFit::Sort(double A, double B, double C, double offset)
{
    double par[9];
    par[0] = A;
    par[1] = B;
    par[2] = C;
    par[3] = fParameters[3];
    par[4] = fParameters[4];
    par[5] = fParameters[5];
    par[6] = fParameters[6];
    par[7] = fParameters[7];
    par[8] = offset;
    
    Sort(par);
}

void HistFit::Sort(double * par)
{
    fRandom.SetSeed(1);
    gErrorIgnoreLevel = kError;    
    
    SetParameters(par);

    fExpBinNum = fExpHist->GetNbinsX();
    if(fSimHist) { delete fSimHist; fSimHist = NULL; }
    fSimHist = new TH1F("fSimHist","fSimHist",fExpBinNum,0,4000); 
    int nHits = 0;
    double light = 0.;
    double centroidEkin = 0.;    
    double centroidEres = 0.;    
    
    int counter = 0;

    for(int i=0; i<fNumEntries; i++)
    {
        counter++;
        if( counter%50000==0 ) std::cout << "sorting evt " << counter << "/" << fNumEntries << "; " << double(counter)/double(fNumEntries)*100 << "% complete \r"  << std::flush; 
     
        fEdepBranch->GetEntry(i);   
        fEkinBranch->GetEntry(i);   
        fPtypeBranch->GetEntry(i);   
        nHits = fEdepVector->size();
        light = 0.;
        for(int j=0; j<nHits; j++)
        {
            if(fPtypeVector->at(j) == 2 || fPtypeVector->at(j) == 3) {
                centroidEkin = fEkinVector->at(j);
                centroidEres = fEkinVector->at(j)-fEdepVector->at(j);
            }
            else if(fPtypeVector->at(j) == 4) {
                centroidEkin = LightOutput(fEkinVector->at(j), fProtonCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fProtonCoeff);
            }
            else if(fPtypeVector->at(j) == 6) {
                centroidEkin = LightOutput(fEkinVector->at(j), fDeuteronCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fDeuteronCoeff);
            }
            else if(fPtypeVector->at(j) == 7) {
                centroidEkin = LightOutput(fEkinVector->at(j), fCarbonCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fCarbonCoeff);
            }
            else if(fPtypeVector->at(j) == 8) {
                centroidEkin = LightOutput(fEkinVector->at(j), fAlphaCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fAlphaCoeff);
            }
            else if(fPtypeVector->at(j) == 9) {
                centroidEkin = LightOutput(fEkinVector->at(j), fBeCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fBeCoeff);
            }
            else if(fPtypeVector->at(j) == 10) {
                centroidEkin = LightOutput(fEkinVector->at(j), fBCoeff );
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fBCoeff);
            }
            else { 
                centroidEkin = 0.; 
                centroidEres = 0.; 
            } 

            if(centroidEkin>0.){
                light += 1000.*fRandom.Gaus(centroidEkin, Resolution(centroidEkin,fSmearingCoeff));
            }
            if(centroidEres>0.){
                light -= 1000.*fRandom.Gaus(centroidEres, Resolution(centroidEres,fSmearingCoeff));
            } 
        }//end scatters loop       
        if(light>0.) fSimHist->Fill(light+fLightOffset);
        //if(light>0.) fSimHist->Fill(light);
    }//end event loop

    ApplyCutoffLow(fCutoffLow,"sim");    
    fSimHist->Scale(fExpHist->Integral(fExpHist->FindBin(fCutoffLow),fExpHist->FindBin(fCutoffHigh),"width")/fSimHist->Integral(fSimHist->FindBin(fCutoffLow),fSimHist->FindBin(fCutoffHigh),"width"));
    fSimHist->SetStats(false);
    std::cout << "sorting... done!                                                   " << std::endl;

}

/////////////////////////////////////////////////////////////////////////////////////

class Fitter 
{

public:
    Fitter();
    ~Fitter();

    Fitter(std::string);
    Fitter(std::string, std::string);
    Fitter(std::string, std::string, std::string);
    Fitter(std::string, std::string, std::string, std::string);
    Fitter(std::string, std::string, std::string, std::string, std::string, std::string);
    Fitter(std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string);
    
    void Draw();
    void Run(double A=0.10, double B=0.05, double C=1e-4, double offset = 0);

    bool Check(int i) { if(i<=-1||i>=GetNumberOfHistFits()) return false; else return true; }
    
    HistFit * GetHistFit(int i) { if(Check(i)) return &fHistFitVector[i]; else{ std::cout << "out of bounds!" << std::endl; return NULL;} }
    
    void SetNextHistFit(HistFit hfit) { fHistFitVector.push_back(hfit); fRunIdVector.push_back(hfit.GetRunId());  }
    void SetNextHistFit(std::string val) {
        HistFit * hfit = new HistFit(val);
        SetNextHistFit(*hfit);
    }    

    int GetNumberOfHistFits() { return fHistFitVector.size(); }
    
    void DrawHistFit(int i) { 
        if(!Check(i)) { std::cout << "out of bounds!" << std::endl; return; }
        if(fCanvas) delete fCanvas; 
        fHistFitVector.at(i).Draw(); 
    }
    
    void SortRun(int num) { 
        if(Check(num)) { 
            fHistFitVector.at(num).Sort(fParameters); 
        }
        else{ std::cout << "out of bounds!" << std::endl;} 
    }
    
    void SortAllRuns() { 
        //PrintParameters();
        for(int num=0; num<GetNumberOfHistFits(); num++) SortRun(num); 
        DoChi2();
    }
    
    void Print() { for(int num=0; num<GetNumberOfHistFits(); num++) std::cout << "Run# = " << fRunIdVector.at(num) << std::endl; }

    void SetParameters(double A=0.10, double B=0.05, double C=1e-4, double a1=0.639, double a2=1.462, double a3=0.373, double a4=0.968, double carbon=0,double offset=0) {   
        fParameters[0]=A;fParameters[1]=B;fParameters[2]=C;
        fParameters[3]=a1; fParameters[4]=a2; fParameters[5]=a3; fParameters[6]=a4; 
        fParameters[7]=carbon; 
        fParameters[8] = offset;
        for(int i=0; i<GetNumberOfHistFits(); i++) fHistFitVector.at(i).SetParameters(fParameters);
    }
    void SetParameters(double * par) { // expects a par array w/ 8 elements
        for(int i=0; i<8; i++) fParameters[i] = par[i];
        for(int i=0; i<GetNumberOfHistFits(); i++) fHistFitVector.at(i).SetParameters(fParameters);
    }
    void SetOffset(double offset) {
        fParameters[8] = offset;
        for(int i=0; i<GetNumberOfHistFits(); i++) fHistFitVector.at(i).SetParameters(fParameters);
    }
    
    void PrintParameters() { 
        std::cout << "      A = " << fParameters[0] << std::endl;
        std::cout << "      B = " << fParameters[1] << std::endl;
        std::cout << "      C = " << fParameters[2] << std::endl;
        std::cout << "     a1 = " << fParameters[3] << std::endl;
        std::cout << "     a2 = " << fParameters[4] << std::endl;
        std::cout << "     a3 = " << fParameters[5] << std::endl;
        std::cout << "     a4 = " << fParameters[6] << std::endl;
        std::cout << " carbon = " << fParameters[7] << std::endl;
        std::cout << " offset = " << fParameters[8] << " keVee" << std::endl;
    }

    void InitializeParameters();

    void NelderMead3(double A=0.10, double B=0.05, double C=0.0001, int itermax=50);

    std::vector<HistFit> fHistFitVector;   
    std::vector<std::string> fRunIdVector;

    double fParameters[9];   
 
    TCanvas * fCanvas;
    
    double fSum;
    double fSum2;
    double DoChi2() { 
        fSum = 0.;
        fSum2 = 0.;
        for(int i=0;i<GetNumberOfHistFits();i++) fSum += fHistFitVector.at(i).DoChi2();
        for(int i=0;i<GetNumberOfHistFits();i++) fSum2 += fHistFitVector.at(i).DoChi2() * fHistFitVector.at(i).DoChi2();
        fSum /= double(GetNumberOfHistFits());
        fSum2 /= double(GetNumberOfHistFits());    
    
        return fSum2; // CHANGE THIS TO SET WHAT WE WILL MINIMIZE [ chi2 or (chi2)^2 ]
    }
        
    double nm_val(double * par) {
        SetParameters(par);
        SortAllRuns();
        return DoChi2();
    }
    double nm_val(vec v) {
        SetParameters(v.par_array());
        SortAllRuns();
        return DoChi2();
    }

    void DrawToFile(std::string name);


};

void Fitter::InitializeParameters()
{
    fParameters[0] = 0.1;
    fParameters[1] = 0.05;
    fParameters[2] = 1e-4;
    fParameters[3] = 0.639;
    fParameters[4] = 1.462;
    fParameters[5] = 0.373;
    fParameters[6] = 0.968;
    fParameters[7] = 0.;
    fParameters[8] = 0;

    fSum = 0;
    fSum2 = 0;
}

void Fitter::Draw()
{
    if(!fCanvas) {
        fCanvas = new TCanvas();
        if(GetNumberOfHistFits() == 1) fCanvas->Divide(1);
        else if(GetNumberOfHistFits() == 2) fCanvas->Divide(2);
        else if(GetNumberOfHistFits() == 3) fCanvas->Divide(3);
        else if(GetNumberOfHistFits() == 4) fCanvas->Divide(2,2);
        else if(GetNumberOfHistFits() == 6) fCanvas->Divide(2,3);
        else if(GetNumberOfHistFits() == 8) fCanvas->Divide(2,4);
        else if(GetNumberOfHistFits() == 10) fCanvas->Divide(2,5);
        else { std::cout << "more/less HistFits that I know how to draw!" << std::endl; return; }
    }
    for(int i=0; i<GetNumberOfHistFits(); i++) {
        fCanvas->cd(i+1);
        gPad->Clear();
        fHistFitVector.at(i).Draw();
        gPad->Update();
    }
}

Fitter::Fitter() { InitializeParameters(); } 

Fitter::~Fitter() {}

Fitter::Fitter(std::string one)
{
    fCanvas = NULL;
    SetNextHistFit(one);
    InitializeParameters();
}

Fitter::Fitter(std::string one, std::string two)
{
    fCanvas = NULL;
    SetNextHistFit(one);
    SetNextHistFit(two);
    InitializeParameters();
}

Fitter::Fitter(std::string one, std::string two, std::string three)
{
    fCanvas = NULL;
    SetNextHistFit(one);
    SetNextHistFit(two);
    SetNextHistFit(three);
    InitializeParameters();
}

Fitter::Fitter(std::string one, std::string two, std::string three, std::string four)
{
    fCanvas = NULL;
    SetNextHistFit(one);
    SetNextHistFit(two);
    SetNextHistFit(three);
    SetNextHistFit(four);
    InitializeParameters();
}

Fitter::Fitter(std::string one, std::string two, std::string three, std::string four, std::string five,  std::string six)
{
    fCanvas = NULL;
    SetNextHistFit(one);
    SetNextHistFit(two);
    SetNextHistFit(three);
    SetNextHistFit(four);
    SetNextHistFit(five);
    SetNextHistFit(six);
    InitializeParameters();
}

Fitter::Fitter(std::string one, std::string two, std::string three, std::string four, std::string five, std::string six, std::string seven, std::string eight)
{
    fCanvas = NULL;
    SetNextHistFit(one);
    SetNextHistFit(two);
    SetNextHistFit(three);
    SetNextHistFit(four);
    SetNextHistFit(five);
    SetNextHistFit(six);
    SetNextHistFit(seven);
    SetNextHistFit(eight);
    InitializeParameters();
}

void Fitter::Run(double A, double B, double C, double offset) 
{
    double a1 = fParameters[3];
    double a2 = fParameters[4];
    double a3 = fParameters[5];
    double a4 = fParameters[6];
    double carbon = fParameters[7];
    offset = fParameters[8];
    
    SetParameters(A,B,C,a1,a2,a3,a4,carbon,offset);
    PrintParameters();
    if(!fCanvas) {
        fCanvas = new TCanvas();
        if(GetNumberOfHistFits() == 1) fCanvas->Divide(1);
        else if(GetNumberOfHistFits() == 2) fCanvas->Divide(2);
        else if(GetNumberOfHistFits() == 3) fCanvas->Divide(3);
        else if(GetNumberOfHistFits() == 4) fCanvas->Divide(2,2);
        else if(GetNumberOfHistFits() == 6) fCanvas->Divide(2,3);
        else if(GetNumberOfHistFits() == 8) fCanvas->Divide(2,4);
        else if(GetNumberOfHistFits() == 10) fCanvas->Divide(2,5);
        else { std::cout << "more/less HistFits that I know how to draw!" << std::endl; return; }
    }
    for(int i=0; i<GetNumberOfHistFits(); i++) {
        SortRun(i);
        fCanvas->cd(i+1);
        gPad->Clear();
        fHistFitVector.at(i).Draw();
        gPad->Update();
    }
    fSum = 0.;
    fSum2 = 0.;
    for(int i=0;i<GetNumberOfHistFits();i++) fSum += fHistFitVector.at(i).DoChi2();
    for(int i=0;i<GetNumberOfHistFits();i++) fSum2 += fHistFitVector.at(i).DoChi2() * fHistFitVector.at(i).DoChi2();
    fSum /= double(GetNumberOfHistFits());
    fSum2 /= double(GetNumberOfHistFits());
    std::cout << "sum(chi2)/nfits = " << fSum << " | sum((chi2)^2)/nfits = " << fSum2 << std::endl;
}

void Fitter::DrawToFile(std::string input)
{
    std::cout << "drawing all HistFits to output file \"" << input << "\" ... " << std::flush;
    TCanvas * out = new TCanvas();
    for(int i=0; i<GetNumberOfHistFits(); i++) 
    {
        std::string name = input;
        fHistFitVector.at(i).Draw();
        if(i==0) {
            name += "(";
            out->Print(name.c_str(),"pdf");
        }
        else if(i==GetNumberOfHistFits()-1) {
            name += ")";
            out->Print(name.c_str(),"pdf");
        }
        else out->Print(name.c_str(),"pdf");
    }
    delete out;
    std::cout << " done!" << std::endl;

}

void Fitter::NelderMead3(double A, double B, double C, int itermax)
{
    std::cout << "starting Nelder Mead method... " << std::endl;
    std::cout << "!!! only fitting the 3 Gaussian parameters : A, B, C" << std::endl;
    
    
    double inc2 = 0.01;   // A
    double inc3 = 0.01;  // B
    double inc4 = 1e-4; // C

    //    ( A   , B    , C     , a1   , a2  , a3  , a4   , carbon)
    vec v1(A,B,C);
    vec v2(v1); v2.set(0,v2.at(0)+inc2);
    vec v3(v1); v3.set(1,v3.at(1)+inc3);
    vec v4(v1); v4.set(2,v4.at(2)+inc4);

    std::vector<vec> nmvec;
    nmvec.push_back(v1);
    nmvec.push_back(v2);
    nmvec.push_back(v3);
    nmvec.push_back(v4);

    std::cout << "calculating chi2's for the initial simplex..." << std::endl;
    std::vector<double> chi2vec;
    for(int i=0; i<int(nmvec.size()); i++) {
        std::cout << "simplex " << i+1 << std::endl;
        SetParameters(nmvec.at(i).par_array());         
        SortAllRuns();
        chi2vec.push_back(DoChi2());
    } 
        
    std::cout << "starting the Nelder-Mead iterations..." << std::endl;
    //////////////////////////////////////////////////////////////////
    for(int iter=1; iter<=itermax; iter++) {

        std::vector<vec> temp_par;
        std::vector<double> temp_chi2;
        temp_par.resize(4);
        temp_chi2.resize(4);

        // reordering...
        double test = 1e100;
        int val = 0;
        for(int i=0; i<4; i++) {
            for(int j=0; j<4; j++) {
                if(chi2vec.at(j) < test) {
                    test = chi2vec.at(j);
                    temp_chi2.at(i) = test;
                    temp_par.at(i) = nmvec.at(j);
                    val = j;
                }
            }
            chi2vec.at(val) = 1e100;
            test = 1e100;
        }
        nmvec = temp_par;
        chi2vec = temp_chi2;

        std::cout << "printing the reordered variables..." << std::endl;
        for(int i=0; i<4; i++) {
            std::cout << " chi2 = " << chi2vec.at(i);
            std::cout << " pars = ";
            for(int j=0; j<2; j++) std::cout << nmvec.at(i).at(j) << " , "; std::cout << nmvec.at(i).at(2);
            std::cout << std::endl;
        }
        
    
        vec B(nmvec.at(0)); double B_chi2 = chi2vec.at(0);
        vec G(nmvec.at(1)); double G_chi2 = chi2vec.at(1);
        vec W(nmvec.at(3)); double W_chi2 = chi2vec.at(3);
        vec M = B.midpoint(G); double M_chi2 = nm_val(M);
        vec R = M.scalar_multiply(2.); R.subtract(W); double R_chi2 = nm_val(R);
        vec E = R.scalar_multiply(2.); E.subtract(M); double E_chi2 = 0;
        vec C; double C_chi2 = 0;        
        vec S; double S_chi2 = 0;    

        // now with the logical decisions....
        if(R_chi2 < G_chi2) {  // case 1
            if(B_chi2 < R_chi2) {
                W = R; W_chi2 = R_chi2;
            }
            else {
                E_chi2 = nm_val(E);
                if(E_chi2 < B_chi2) {
                    W = E; W_chi2 = E_chi2;   
                }
                else {
                    W = R; W_chi2 = R_chi2;
                }        
            }
        }
        else {  // case 2
            if(R_chi2 < W_chi2) {
                W = R; W_chi2 = R_chi2;
            }
            vec C1 = W.midpoint(M); double C1_chi2 = nm_val(C1);
            vec C2 = M.midpoint(R); double C2_chi2 = nm_val(C2);
            if(C1_chi2 < C2_chi2) { C = C1; C_chi2 = C1_chi2; }
            else                  { C = C2; C_chi2 = C2_chi2; }
        
            if(C_chi2 < W_chi2) {
                W = C; W_chi2 = C_chi2;
            }
            else {
                S = B.midpoint(W); S_chi2 = nm_val(S);
                W = S; W_chi2 = S_chi2;
                G = M; G_chi2 = M_chi2;
            }
        }
        nmvec.at(0) = B; chi2vec.at(0) = B_chi2;
        nmvec.at(1) = G; chi2vec.at(1) = G_chi2;
        nmvec.at(3) = W; chi2vec.at(3) = W_chi2;
    
        std::cout << std::endl << "finished iteration # " << iter << "/" << itermax << std::endl << std::endl;
        
        // end of logical loop
    }
    //////////////////////////////////////////////////////////////////
}
