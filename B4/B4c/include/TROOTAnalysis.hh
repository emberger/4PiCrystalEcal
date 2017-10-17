
#ifndef TROOTAnalysis_hh
#define TROOTAnalysis_hh


#include "B4ROOTEvent.hh"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TBrowser.h"
#include "TColor.h"
#include "TH3D.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraph2D.h"
#include "TAttMarker.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TImage.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MinimumParameters.h"
#include "Minuit2/MnUserCovariance.h"
#include "Minimizer.hh"
#include "Minuit2/Minuit2Minimizer.h"

#include <chrono>
#include <string>
#include <cstring>

using namespace ROOT::Math;
using namespace ROOT::Minuit2;



class TROOTAnalysis {

public:
TROOTAnalysis(TChain* ch);
~TROOTAnalysis();
Int_t GetNofEntries();

void PrintERes();
void PlotRMSx();
void plotEvent(Int_t pev);
void plotEventPion(Int_t pev);
//void plotCOGs();
void PlotProjection(Int_t distance);
void PlotChimap(Int_t event);

//void AnalyzePions

Bool_t CheckEvent(Int_t event);

void CalcCOGPion(Int_t event);

std::tuple<Double_t, Double_t, Double_t> TransformCoordinates(Double_t x, Double_t y, Double_t z, Double_t seg);

void FitCOGsPion(Int_t event);
void FitCOGsPion2(Int_t event);


std::pair<TVector3, TVector3> FindClosestApproach(Int_t event);

void GetInvariantMass(Int_t event);

//void MinimizeClosestApproach(Int_t event, std::pair<TVector3, TVector3> p);

void PionLocator(Int_t event, std::pair<TVector3, TVector3> p);


void SampleFromHIst();



void DrawHists();

/////////////////////////////////////////////////////////////////////////////////////////////
// void CalcCOG(Int_t minevent, Int_t maxevent);
// void FitCOGs( Int_t minevent, Int_t maxevent);
//
// void CalcCOG(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent);
//////////////////////////////////////////////////////////////////////////////////////////////






//void CalcCOGwithFit(Int_t minlayer, Int_t maxlayer);
//void CleanCOGs(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent);


void PrintFitHists(Int_t minevent, Int_t maxevent);



void SetPath(std::string path){
        savepath=path;
        pathset=true;
}


private:
Int_t nofRawEntries;
Int_t nofEntries;     // number of events in Tree

Double_t Eges;  //Energy in event

Double_t EPhot1;
Double_t EPhot2;



Double_t tiledimX;
Double_t tiledimY;
Double_t calsizeXY;
Double_t AbsoThickness;
Double_t GapThickness;
Double_t nofLayers;

Double_t crystalThickness;
Double_t pcbThickness;
Double_t crystalsizeXY;

TVector3 gunposition;

Double_t histsizeX;
Double_t histsizeY;
Double_t histsizeZ;

std::string savepath;

Bool_t pathset=false;

std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > coglist;

std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > > COGCollectionPH1;
std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > > COGCollectionPH2;

std::vector<std::tuple<Double_t, Double_t, Double_t> > showerCOGPhoton1;
std::vector<std::tuple<Double_t, Double_t, Double_t> > showerCOGPhoton2;


std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > FitParamsGamma;
std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > FitParams;//to be deleted, needed for old analysis
std::vector<std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > > FitParamsPions;

std::vector<std::pair< Double_t, Int_t> > EnergyPhoton1;
std::vector<std::pair<TVector3, TVector3 > > DirectionPhoton1;

std::vector<std::pair< Double_t, Int_t> > EnergyPhoton2;
std::vector<std::pair<TVector3, TVector3 > >  DirectionPhoton2;


std::vector<Double_t> InvariantMass;

std::vector<std::tuple<Double_t,Double_t,Double_t, Double_t> > ClosestApproach;

std::vector<std::tuple<Double_t, Double_t, Double_t> > DeviationFromGun;

std::vector<std::tuple<Double_t, Double_t, Double_t> > ClusteredHits;
std::vector<Double_t> showerCenters;

TTree* EcalTree;


B4ROOTEvent * Cevent;


//TCanvas * c2 = new TCanvas("COGs", "COGs");

TH1D * fitX = new TH1D("X Fit","X Fit", 4000,-600,600);
TH1D * fitY = new TH1D("Y Fit","Y Fit", 4000,-600,600);

TH1D * fitSX = new TH1D("X Slope Fit","X Slope Fit", 800,-2,2);
TH1D * fitSY = new TH1D("Y Slope Fit","Y Slope Fit", 800,-2,2);


TH1D * er1= new TH1D("errx", "errx", 50,0,50);
TH1D * er2= new TH1D("erry", "erry", 50,0,50);



TH3D * h;


TH3D * hA;
TH3D * hB;
TH2D * h1;
TH2D * h2;
TH3D * h3;
//TH3D * h4;



//TCanvas * check= new TCanvas("CHECK", "CHECK");

//Hists for ClosestApproach difference



//Direction difference
TCanvas * D = new TCanvas("DC","Direction vector difference");
TH1D * dx1= new TH1D("DeltaX_photon1","DeltaX_photon1", 1000,-2,2);
TH1D * dy1= new TH1D("DeltaY_photon1","DeltaY_photon1", 1000,-2,2);
TH1D * dz1= new TH1D("DeltaZ_photon1","DeltaZ_photon1", 1000,-2,2);
TH1D * dx2= new TH1D("DeltaX_photon2","DeltaX_photon2", 1000,-2,2);
TH1D * dy2= new TH1D("DeltaY_photon2","DeltaY_photon2", 1000,-2,2);
TH1D * dz2= new TH1D("DeltaZ_photon2","DeltaZ_photon2", 1000,-2,2);

TCanvas * Delta1 = new TCanvas("Delta1", "difference of gun positoin to ClosestApproach");
// TCanvas * Delta2 = new TCanvas("Delta2", "difference of gun positoin to ClosestApproach");
// TCanvas * Delta3 = new TCanvas("Delta3", "difference of gun positoin to ClosestApproach");
// TCanvas * Delta4 = new TCanvas("Delta4", "difference of gun positoin to ClosestApproach");

TH1D * delx = new TH1D("DeltaX_VertexReconstruction", "DeltaX_VertexReconstruction", 1000,-1000,1000);
TH1D * dely = new TH1D("DeltaY_VertexReconstruction", "DeltaY_VertexReconstruction", 1000,-1000,1000);
TH1D * delz = new TH1D("DeltaZ_VertexReconstruction", "DeltaZ_VertexReconstruction", 1000,-1000,1000);
TH1D * appdist1=new TH1D("Closest Approach Distance", "Closest Approach Distance", 500,0,500);

//Hists for InvariantMassCalculation

TCanvas * IM1 = new TCanvas("InvariantMassC");


TH1D * InvMassReco = new TH1D("InvariantMassH1", "MC Energy, reconstructed direction",500, 0,1000);
TH1D * InvMassSim = new TH1D("InvariantMassH2", "MC Direction, reconstructed Energy", 500,0,1000);
TH1D * InvMassAllreco = new TH1D("InvariantMassH3", "reconstructed Direction, reconstructed Energy", 500,0,1000);
TH2D * Angle_vs_massReco = new TH2D("Angle_vs_massRecoH", "Angle_vs_InvariantMass", 500,0,1000,500,0,200);
TH2D * Angle_vs_massSim = new TH2D("Angle_vs_massSimH", "Angle_vs_InvariantMass", 500,0,1000,500,0,200);
TH2D * Angle_vs_massAllreco = new TH2D("Angle_vs_massAllrecoH", "Angle_vs_InvariantMass", 5000,0,1000,500,0,200);

TCanvas * PiLoc = new TCanvas("MinimizeClosestApproachC");
// TCanvas * PiLoc1 = new TCanvas("MinimizeClosestApproachC1");
// TCanvas * PiLoc2 = new TCanvas("MinimizeClosestApproachC2");

TH1D * pilocX = new TH1D("MinimizeClosestApproachX", "MinimizeClosestApproachX",1200,-1500,1500);
TH1D * pilocY = new TH1D("MinimizeClosestApproachY", "MinimizeClosestApproachY",1200,-1500,1500);
TH1D * pilocZ = new TH1D("MinimizeClosestApproachZ", "MinimizeClosestApproachZ",1200,-1500,1500);

TH1D * inputX = new TH1D("InputLocationX", "InputLocationX",1200,-1500,1500);
TH1D * inputY = new TH1D("InputLocationY", "InputLocationY",1200,-1500,1500);
TH1D * inputZ = new TH1D("InputLocationZ", "InputLocationZ",1200,-1500,1500);

TH1D * invmasFitX= new TH1D("FitInvariantMassX","FitInvariantMassX", 1200,-1500,1500);
TH1D * invmasFitY= new TH1D("FitInvariantMassY","FitInvariantMassY", 1200,-1500,1500);
TH1D * invmasFitZ= new TH1D("FitInvariantMassZ","FitInvariantMassZ", 1200,-1500,1500);


TCanvas * DeltaR= new TCanvas("DeltaR", "3D deviation from gun position");
TH1D * delR= new TH1D("3D deviation from gun position","3D deviation from gun position", 500,0,700);

TCanvas * DeltaR_fit= new TCanvas("DeltaR_fit", "3D deviation from gun position_fit");
TH1D * delR_fit= new TH1D("3D deviation from gun position_fit","3D deviation from gun position_fit", 500,0,700);

TCanvas * ang1=new TCanvas("angle", "anglecheck");

TH1D * angh_rec= new TH1D("angle", "angle", 100, 0,180);
TH1D * angh_true= new TH1D("angle_true", "angle_true", 100, 0,180);

TCanvas * totalEnergyC = new TCanvas("total Energy", "total Energy");
TH1D * totalEnergyH = new TH1D("totalEnergyH", "total Energy",1000,0,1);


};
#endif
