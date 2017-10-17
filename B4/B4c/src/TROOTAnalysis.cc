#include "TROOTAnalysis.hh"




TROOTAnalysis::TROOTAnalysis(TChain* ch){
        this->EcalTree=ch->GetTree();
        //std::cout << "assigned Tree" << '\n';
        nofEntries=EcalTree->GetEntries();
        //std::cout<<nofEntries<<std::endl;
        FitParamsPions.resize(nofEntries);

        Cevent = new B4ROOTEvent();
        EcalTree->SetBranchAddress("EventBranch", &Cevent);

        EcalTree->GetEntry(0);

        GapThickness=Cevent->GapThickness();
        AbsoThickness=Cevent->AbsoThickness();
        tiledimX=Cevent->TilesizeX();
        tiledimY=Cevent->TilesizeY();
        calsizeXY=Cevent->calsizeXY();
        nofLayers=Cevent->NumberOfLayers();
        crystalThickness=Cevent->crystalThickness();
        pcbThickness=Cevent->pcbThickness();
        crystalsizeXY=Cevent->crystalsizeXY();

        //Double_t EcalSizeXYZ=calsizeXY/2+(GapThickness+AbsoThickness)*nofLayers;

        histsizeX=calsizeXY/tiledimX;
        histsizeY=calsizeXY/tiledimY;
        histsizeZ=nofLayers;

        h = new TH3D("ECalEvent","ECalEvent",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);


        hA= new TH3D("ECalEvent1","ECalEvent1",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);
        hB= new TH3D("ECalEvent2","ECalEvent2",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);

        h1= new TH2D("h1", "h1", histsizeX+1,-0.5,histsizeX+0.5,histsizeY+1,-0.5,histsizeY+0.5);
        h2= new TH2D("h2", "h2", histsizeX+1,-0.5,histsizeX+0.5,histsizeY+1,-0.5,histsizeY+0.5);

        h3 = new TH3D("h3", "h3",histsizeX, 0, histsizeX, histsizeY,0, histsizeY, histsizeZ,0,histsizeZ);
        //
        // h3->GetXaxis()->SetTitle("X");
        // h3->GetYaxis()->SetTitle("Y");
        // h3->GetZaxis()->SetTitle("Z");

        //  std::cout<<"GT: "<<GapThickness<<" AT: "<<AbsoThickness<<" TDX: "<<tiledimX<<" TDY: "<<tiledimY<<" CXY: "<<calsizeXY <<" NL: "<<nofLayers<<std::endl;
        //  std::cout<<" HX: "<<histsizeX<<" HY: "<<histsizeY<<" HZ: "<<histsizeZ<<std::endl;
        //EcalTree->Print();
}

//----------------------------------------------------------------------------------------------------------------------------------------------------

TROOTAnalysis::~TROOTAnalysis(){
}

Int_t TROOTAnalysis::GetNofEntries(){
        return nofEntries;
}


void TROOTAnalysis::plotEvent(Int_t pev){     //plot 3DHisto of selected event
        TCanvas * plotcanvas1 = new TCanvas("eventplotter", "eventplotter");

        if(nofEntries<pev+1) {
                std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nofEntries<<"."<<std::endl;
                return;
        }



        //for(Int_t i=0;i<nent; i++){

        EcalTree->GetEntry(pev);
        Int_t pnh=Cevent->NHits();

        for(Int_t j=0; j<pnh; j++) {
                h->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
        }

        h->GetXaxis()->SetTitle("X");
        h->GetYaxis()->SetTitle("Y");
        h->GetZaxis()->SetTitle("Z");

        plotcanvas1->cd();
        h->Draw("BOX");

        std::string Path="/home/iwsatlas1/emberger/FuerFrank/PhotonEventDisplay/Event";
        std::string nr = std::to_string(pev);

        std::string extension = ".C";
        std::string PlotPath=Path+nr+extension;
        plotcanvas1->Print(PlotPath.c_str());

        extension = ".pdf";
        PlotPath=Path+nr+extension;
        plotcanvas1->Print(PlotPath.c_str());
        h->Reset();
}

void TROOTAnalysis::plotEventPion(Int_t pev){     //plot 3DHisto of selected event
        TCanvas * plotcanvas1 = new TCanvas("eventplotter", "eventplotter");
        if(nofEntries<pev+1) {
                std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nofEntries<<"."<<std::endl;
                return;
        }



        //for(Int_t i=0;i<nent; i++){

        EcalTree->GetEntry(pev);
        Int_t pnh=Cevent->NHits();

        for(Int_t j=0; j<pnh; j++) {
                if(Cevent->Hit(j)->PhotNr()==1) {
                        hA->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
                }
        }

        //TColor * col1= new TColor(1,0,0,"red");

        hA->GetXaxis()->SetTitle("X");
        hA->GetYaxis()->SetTitle("Y");
        hA->GetZaxis()->SetTitle("Z");
        hA->SetMarkerColor(kRed);
        hA->SetMarkerStyle(20);
        hA->SetMarkerSize(0.3);
        hA->Draw("");

        for(Int_t j=0; j<pnh; j++) {
                if(Cevent->Hit(j)->PhotNr()==2) {
                        hB->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
                }
        }

        //TColor * col2 = new TColor(0,0,1,"blue");
        hB->SetMarkerColor(kBlue);
        hB->SetMarkerStyle(20);
        hB->SetMarkerSize(0.3);
        hB->GetXaxis()->SetTitle("X");
        hB->GetYaxis()->SetTitle("Y");
        hB->GetZaxis()->SetTitle("Z");
        hB->Draw("SAME");

        std::string Path="/home/iwsatlas1/emberger/PionPlots/Event";
        std::string nr = std::to_string(pev);

        std::string extension = ".C";
        std::string PlotPath=Path+nr+extension;
        plotcanvas1->Print(PlotPath.c_str());

        // extension = ".png";
        // PlotPath=Path+nr+extension;
        // plotcanvas1->Print(PlotPath.c_str());

        extension = ".pdf";
        PlotPath=Path+nr+extension;
        plotcanvas1->Print(PlotPath.c_str());

        hB->Reset();
        hA->Reset();
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::PlotRMSx(){
        TCanvas * rms=new TCanvas("rmsx_over_E");


        // Double_t y[10]={0.5298,0.3507,0.3624,0.2524,0.222,0.1861,0.228,0.1751,0.1425,0.1481}; // 0.0 slope
        // Double_t x[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        //
        Double_t y[10]={0.5398,0.3428,0.3087,0.2319,0.2117,0.1933,0.2243,0.1716,0.1411,0.1588};                                                                      //0.4slope
        Double_t x[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        TF1 * rms_x = new TF1("rms", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");

        TGraph *gRMS1=new TGraph(10, x, y );
        gRMS1->Fit(rms_x);
        gStyle->SetOptFit(11111111);
        gRMS1->GetXaxis()->SetTitle("Energy[GeV]");
        gRMS1->GetYaxis()->SetTitle("#frac{RMS_x}{100mm}");
        gRMS1->GetYaxis()->SetTitleOffset(1.3);
        rms->cd(0);
        gRMS1->Draw("A*");
}


//---------------------------------------------------------------------------------------------------------------------------------------------------------

void TROOTAnalysis::SampleFromHIst(){

        TCanvas * dummy = new TCanvas();
        dummy->Divide(2,2,0.01,0.01);
        TFile *f = new TFile("/home/iwsatlas1/emberger/Geant4/Current/4piSensitiveDetector_mod/B4-build/pi0_kinematics.root");
        //f.ls();
        Double_t energy;
        Double_t cos_theta;

        TH2D * hist1 = (TH2D*)f->Get("pi0_energyAngle");

        TH1D * sample1= new TH1D("sample1", "sample1", 100, 0, 1);
        TH1D * sample2= new TH1D("sample2", "sample2", 200, -1, 1);

        TRandom3 * rndGen = new TRandom3();

        TVector3 BeamDirection(0.,0.,1.);
        BeamDirection=BeamDirection.Unit();

        TH1D * check = new TH1D("Check", "check", 200, -1, 1);
        TH2D * fd = new TH2D("fd", "fd", 200,-2,2,200,-2,2);
        for(Int_t i=0; i<1000000; i++) {
                hist1->GetRandom2(energy, cos_theta);
                sample1->Fill(energy);
                sample2->Fill(cos_theta);


                if(cos_theta>0) {
                        Double_t theta = TMath::ACos(cos_theta);

                        Double_t radius = TMath::Tan(theta);

                        Double_t a=rndGen->Uniform(0,TMath::TwoPi());

                        Double_t x=TMath::Cos(a)*TMath::Abs(radius);
                        Double_t y=TMath::Sin(a)*TMath::Abs(radius);
                        fd->Fill(x,y);
                        TVector3 MomentumDirection(x+BeamDirection.x(), y+BeamDirection.y(), BeamDirection.z());
                        MomentumDirection=MomentumDirection.Unit();
                        check->Fill(TMath::Cos(MomentumDirection.Angle(BeamDirection)));
                }
                else if(cos_theta==0) {
                        //Double_t theta = TMath::ACos(cos_theta);

                        Double_t a=rndGen->Uniform(0,TMath::TwoPi());

                        Double_t x=TMath::Cos(a)*1;
                        Double_t y=TMath::Sin(a)*1;
                        fd->Fill(x,y);

                        TVector3 MomentumDirection(x+BeamDirection.x(), y+BeamDirection.y(), 0);
                        MomentumDirection=MomentumDirection.Unit();
                        check->Fill(TMath::Cos(MomentumDirection.Angle(BeamDirection)));
                }
                else if(cos_theta<0) {
                        Double_t theta = TMath::ACos(cos_theta);

                        Double_t radius = TMath::Tan(theta);

                        Double_t a=rndGen->Uniform(0,TMath::TwoPi());

                        Double_t x=TMath::Cos(a)*TMath::Abs(radius);
                        Double_t y=TMath::Sin(a)*TMath::Abs(radius);
                        fd->Fill(x,y);

                        TVector3 MomentumDirection(x+BeamDirection.x(), y+BeamDirection.y(), -BeamDirection.z());
                        MomentumDirection=MomentumDirection.Unit();
                        check->Fill(TMath::Cos(MomentumDirection.Angle(BeamDirection)));
                }

        }

        dummy->cd(1);
        sample1->Draw();
        dummy->cd(2);
        sample2->Draw();
        dummy->cd(3);
        check->Draw();
        dummy->cd(4);
        fd->Draw();


}



void TROOTAnalysis::PrintERes(){

        TCanvas * res = new TCanvas("Energy Resolution", "ERes",2000,1000);
        TCanvas * gap = new TCanvas("GapEnergy", "GapEnergy");

        Double_t y[12]={4.55/19.25, 5.942/38.47, 9.056/77.51, 12.06/115.9, 11.59/154.2, 13.71/192.9, 17.98/268.8, 20.25/384, 53.25/956.2, 69.23/1905, 94.6/2836, 146.7/3779};
        Double_t x[12]={0.05,0.100,0.200,0.300,0.400,0.500,0.700,1.000,2.500,5.000,7.500,10.000};

        TGraph * re1 = new TGraph(12, x, y);
        re1->SetTitle("Energyresolution");
        TF1 * EnergyRes = new TF1("EnergyRes", "sqrt(([0] * [0] / x) + ([1]*[1])  + ([2]*[2]/(x*x)))");

        re1->Fit(EnergyRes);



        re1->GetXaxis()->SetTitle("Energy[GeV]");
        re1->GetYaxis()->SetTitle("#frac{#sigma}{E}");
        re1->GetYaxis()->SetTitleOffset(1.3);
        //re1->GetYaxis()->LabelsOption("v");
        //re1->SetMarkerStyle(23);
        gStyle->SetOptFit();
        res->cd(0);
        re1->Draw("A*");
        TImage * img1 = TImage::Create();
        img1->FromPad(res);
        //img1->Scale(1000, 1000);
        img1->WriteImage("/home/iwsatlas1/emberger/Geant4/Current/SensitiveDetector/B4-build/B4c/GammaEnergyandSlopeScan_Analysis/ERes.png");
        delete img1;

        Double_t Egun[12]={50,100,200,300,400,500,700,1000,2500,5000,7500,10000};
        Double_t Egap[12]={19.25,38.47,77.51,115.9,154.2,192.9,268.8,384,956.2,1905,2836,3779};

        Double_t errx[12]={0,0,0,0,0,0,0,0,0,0,0,0};
        Double_t erry[12]={4.55, 5.942, 9.056, 12.06, 11.59, 13.71, 17.98, 20.25, 53.25, 69.23, 94.6, 146.7};

        TF1 * Efficiency = new TF1("Efficiency", "[0]*x+[1]");

        TGraphErrors * gap1 = new TGraphErrors(12, Egun, Egap, errx, erry );
        gap1->SetTitle("Efficiency");
        gap1->Fit(Efficiency);
        gap1->GetXaxis()->SetTitle("Gun Energy[MeV]");
        gap1->GetYaxis()->SetTitle("Gap Energy[MeV]");
        gap1->GetYaxis()->SetTitleOffset(1.3);
        gStyle->SetOptFit();

        gap->cd(0);
        gap1->Draw("A*");
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::PlotProjection(Int_t distance){

        TCanvas * dist1= new TCanvas("Projection", "Projection", 1500,1000);
        TCanvas * dist2 = new TCanvas("ProjectionXYC", "ProjectionXYC", 1000,1000);


        dist1->Divide(2,1,0.01,0.01);

        TH1D * dx = new TH1D("Projection X", "Projection X", 1000, -600,600);
        dx->GetXaxis()->SetTitle("X[mm]");
        TH1D * dy = new TH1D("Projection Y", "Projection Y", 1000, -600,600);
        dy->GetXaxis()->SetTitle("Y[mm]");

        TH2D * cont1 = new TH2D("ProjectionXY", "ProjectionXY", 2000, -1000,1000,2000,-1000,1000);
        cont1->GetXaxis()->SetTitle("X[mm]");
        cont1->GetYaxis()->SetTitle("Y[mm]");
        std::cout<<"Created Hists"<<std::endl;

        for(Int_t i=0; i<nofEntries; i++) {
                Double_t proX=std::get<0>(FitParamsGamma[i])*(-275-distance)+std::get<1>(FitParamsGamma[i]);
                Double_t proY=std::get<2>(FitParamsGamma[i])*(-275-distance)+std::get<3>(FitParamsGamma[i]);

                dx->Fill(proX);
                dy->Fill(proY);
                cont1->Fill(proX,proY);
                std::cout<<"filled: "<<i<<std::endl;
        }

        gStyle->SetStatY(0.9);
        // Set y-position (fraction of pad size)
        gStyle->SetStatX(0.9);
        // Set x-position (fraction of pad size)

        dist1->cd(1);
        dx->Draw("");

        dist1->cd(2);
        dy->Draw("");

        dist2->cd(0);
        gStyle->SetStatY(0.9);
        // Set y-position (fraction of pad size)
        gStyle->SetStatX(0.9);
        // Set x-position (fraction of pad size)
        cont1->Draw("colz");


        if(pathset) {

                // std::string imgname1= "Projection100mm.png";
                // std::string imgpath1 = savepath + '/' + imgname1;

                std::string imgname2= "ProjectionColour100mm";
                std::string imgpath2c = savepath + '/' + imgname2;


                std::string imgpath2=imgpath2c + ".C";
                dist2->Print(imgpath2.c_str());


                imgpath2=imgpath2c + ".pdf";
                dist2->Print(imgpath2.c_str());




                // TImage * img3 = TImage::Create();
                // TImage * img4 = TImage::Create();
                //
                // img3->FromPad(dist1);
                // img4->FromPad(dist2);
                //
                // img3->WriteImage(imgpath1.c_str());
                // img4->WriteImage(imgpath2.c_str());
                // delete img3;
                // delete img4;
        }

}


//--------------------------------------------------------------------------------------------------------------------------------


// void TROOTAnalysis::plotCOGs(){
//
//         h3->GetXaxis()->SetTitle("X");
//         h3->GetYaxis()->SetTitle("Y");
//         h3->GetZaxis()->SetTitle("Z");
//         h3->SetMarkerStyle(4);
//         c2->cd();
//         h3->Draw("");
//         //h3->Reset();
//
// }
//-------------------------------------------------------------------------------------------------

Bool_t TROOTAnalysis::CheckEvent(Int_t event){
        EcalTree->GetEntry(event);
        gunposition=Cevent->GunPos();

        Int_t hitnr = Cevent->NHits();
        Double_t tmpE1=0;
        Double_t tmpE2=0;
        Double_t seglistph1[6]={0,0,0,0,0,0};
        Double_t seglistph2[6]={0,0,0,0,0,0};

        for(Int_t j=0; j<hitnr; j++) {

                if(Cevent->Hit(j)->PhotNr()==1) {

                        Double_t e1=Cevent->Hit(j)->EnergyDeposit();
                        tmpE1+=e1;

                        Int_t seg1=Cevent->Hit(j)->CalorimeterSegment()-1;
                        seglistph1[seg1]+=e1;
                }
                else if(Cevent->Hit(j)->PhotNr()==2) {
                        Double_t e2=Cevent->Hit(j)->EnergyDeposit();
                        tmpE2+=e2;

                        Int_t seg2=Cevent->Hit(j)->CalorimeterSegment()-1;
                        seglistph2[seg2]+=e2;
                }
        }

        //std::cout<<"energy: "<<seglistph1[0]<<":"<<tmpE1<<"--"<<seglistph2[0]<<":"<<tmpE2<<std::endl;

        Double_t eseg1=0;
        Double_t eseg2=0;
        Int_t segph1=0;
        Int_t segph2=0;

        for(Int_t k =0; k<6; k++) {

                if(seglistph1[k]>eseg1) {
                        segph1=k;
                        eseg1=seglistph1[k];
                }
                if(seglistph2[k]>eseg2) {
                        segph2=k;
                        eseg2=seglistph2[k];
                }
        }
        //std::cout<<segph1+1<<":"<<segph2+1<<std::endl;

        if(seglistph1[segph1]==tmpE1 && seglistph2[segph2]==tmpE2) {

                EnergyPhoton1.push_back(std::make_pair(tmpE1/0.50, segph1+1));
                EnergyPhoton2.push_back(std::make_pair(tmpE2/0.50, segph2+1));

                return true;
                //std::cout<<"true"<<std::endl;
        }
        else{
                return false;
                //std::cout<<"false"<<std::endl;

        }
        //std::cout<<segph1+1<<":"<<segph2+1<<std::endl;
}



//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::CalcCOGPion(Int_t event){            //calculate vector of (X,Y) tuples containing layerwise center of gravity

        // //variables for clustering
        // Int_t xdir=1; // integers for direction
        // Int_t ydir=0;
        // Int_t buf;
        // Int_t stepstodo=1; // how many steps to go before rotation
        // Int_t stepsctr=0; // number of steps since last rotation
        //
        // Int_t currX=0; // starting the spiral on the histogram maximum
        // Int_t currY=0;
        //
        // Double_t esum=0;
        // Double_t curre=0;

        //Int_t maxbin=0;
        //Double_t maxdep=0;
        //Int_t binx, biny, binz;

        // Variables for fit
        Double_t xerr=0;
        Double_t yerr=0;
        Double_t cgx=0;
        Double_t cgy=0;
        Double_t cgz=0;
        Double_t Eweight=0;

        Int_t segment;

        Double_t EGesPH1=0;
        Double_t EGesPH2=0;

        //  for(Int_t eventstodo=minevent; eventstodo<maxevent; eventstodo++) { //loop over all simulated events

        EcalTree->GetEntry(event);         //grab event from tree
        Eges = Cevent->GapEnergy();
        Int_t cnh = Cevent->NHits();
        //std::cout<<"NHits: "<<cnh<<std::endl;
        Double_t integral;
        for(Int_t phnr=1; phnr<3; phnr++) {
                if(phnr==1) {segment=EnergyPhoton1[event].second; }
                if(phnr==2) {segment=EnergyPhoton2[event].second; }

                std::cout<<"segment: "<<segment<<std::endl;

                for(Int_t i=0; i<nofLayers; i++) { //loop over all layers in event

                        for(Int_t j =0; j<cnh; j++) { //loop over all hits in laver i

                                //std::cout<<Cevent->Hit(j)->PhotNr()<<std::endl;
                                Double_t edep= Cevent->Hit(j)->EnergyDeposit();
                                if(Cevent->Hit(j)->Z()==i && Cevent->Hit(j)->PhotNr()==phnr && edep > 0.1) {

                                        //std::cout<<"fill: "<<Cevent->Hit(j)->PhotNr()<<"from photon mode:"<<phnr<<std::endl;
                                        h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), edep);

                                } //fill histogram with hits of layer i
                        }
                        //std::cout<<photonNR<<" : "<<std::endl;
                        integral=h1->Integral();
                        //std::cout<<"integral: "<<integral<<std::endl;
                        //std::cout<<"doing event"<<event<<", photon "<<phnr<<std::endl;
                        //maxbin=h1->GetMaximumBin();
                        //  maxdep=h1->GetBinContent(maxbin); // get coordinates of bin with maximum energy deposition
                        //  h1->GetBinXYZ(maxbin, binx, biny, binz);

                        //Do spiral for clustering only if energy in Layer

                        if(integral!=0) { //check if layer containes energy
                                //std::cout<<"srg"<<std::endl;
                                // xdir=1;  // integers for direction
                                // ydir=0;
                                // stepstodo=1; // how many steps to go before rotation
                                // stepsctr=0; // number of steps since last rotation
                                //
                                // currX=binx; // starting the spiral on the histogram maximum
                                // currY=biny;
                                //
                                // esum=0;
                                // curre=0;
                                // esum+=maxdep;
                                // h2->SetBinContent(currX, currY, maxdep);
                                // auto tp=std::make_tuple(currX, currY, maxdep);
                                // ClusteredHits.push_back(tp);
                                // h1->SetBinContent(binx, biny, 0);

                                // while(esum<integral*0.9 || currX<binx+4) {
                                //         //do spiral until desired energyfraction is reached
                                //         currX+=xdir;
                                //         currY+=ydir;
                                //         curre=h1->GetBinContent(currX, currY);
                                //
                                //         if(curre>0) { //save tile only if it containes energy
                                //                 tp=std::make_tuple(currX, currY, curre);
                                //                 ClusteredHits.push_back(tp);
                                //                 h2->SetBinContent(currX, currY, curre);
                                //         }
                                //         h1->SetBinContent(currX, currY, 0);
                                //         esum+=curre;
                                //         stepsctr++;
                                //         if(stepsctr==stepstodo) {
                                //                 stepsctr=0;
                                //                 buf=xdir; //rotate 90 deg
                                //                 xdir= -ydir;
                                //                 ydir=buf;
                                //
                                //                 if(ydir==0) { //incremant steps at every iteration along x
                                //                         stepstodo++;
                                //                 }
                                //         }
                                // }

                                Double_t e=h1->Integral(); //get energy contained in cluster


                                if(phnr==1) {
                                        Eweight=e/Cevent->EnergyPhoton1();
                                        //calculate weight
                                        EGesPH1+=e;
                                }
                                else if(phnr==2) {
                                        Eweight=e/Cevent->EnergyPhoton2();
                                        EGesPH2+=e;
                                }

                                cgx=h1->GetMean(1);
                                xerr=h1->GetMeanError(1);
                                //extract center of gravity and error
                                cgy=h1->GetMean(2);
                                yerr=h1->GetMeanError(2);

                                if(yerr<0.00001) {yerr=1/TMath::Sqrt(12); }
                                if(xerr<0.00001) {xerr=1/TMath::Sqrt(12); }

                                cgz=i;

                                er1->Fill(i, xerr/nofEntries);
                                er2->Fill(i, yerr/nofEntries);
                                //  std::cout<<"segment: "<<segment<<std::endl;

                                std::tuple<Double_t, Double_t, Double_t> COGtransformed;

                                h3->Fill(cgx,cgy,cgz);

                                COGtransformed=TransformCoordinates(cgx, cgy, cgz, segment);

                                //  std::cout<<std::get<0>(COGtransformed)<<":"<< std::get<1>(COGtransformed)<<":"<< std::get<2>(COGtransformed)<<std::endl;

                                //  std::cout<<"------------------------"<<std::endl;
                                //std::cout<<cgx<<":"<<cgy<<":"<<cgz<<std::endl;

                                auto cg=std::make_tuple(std::get<0>(COGtransformed), std::get<1>(COGtransformed), std::get<2>(COGtransformed), xerr, yerr, Eweight);

                                //  std::cout<<"------------------------"<<std::endl;

                                coglist.push_back(cg);

                        } // end of clustering


                        ClusteredHits.clear();

                        h1->Reset();

                        h2->Reset();

                }         //end of event

                if(phnr==1) {

                        showerCOGPhoton1.push_back(TransformCoordinates(h3->GetMean(1), h3->GetMean(2), h3->GetMean(3), segment));
                        COGCollectionPH1.push_back(coglist);
                        EnergyPhoton1[event]=(std::make_pair(EGesPH1/0.52, segment));
                        coglist.clear();
                        //EnergyPhoton1.push_back(EGesPH1/0.50);
                        std::cout<<"111"<<std::endl;
                }

                else if(phnr==2) {
                        showerCOGPhoton2.push_back(TransformCoordinates(h3->GetMean(1), h3->GetMean(2), h3->GetMean(3), segment));
                        COGCollectionPH2.push_back(coglist);

                        EnergyPhoton2[event]=(std::make_pair(EGesPH2/0.52, segment));
                        coglist.clear();
                        //EnergyPhoton2.push_back(EGesPH2/0.50);
                        std::cout<<"222"<<std::endl;
                }
                h3->Reset();

        }

        totalEnergyH->Fill((EGesPH1+EGesPH2)/1000);

        //check->cd(0);
        //h3->GetYaxis()->SetRange(40,60);
        //h3->Draw("*");
}

//-------------------------------------------------------------------------------------------------------------------------------


std::tuple<Double_t, Double_t, Double_t > TROOTAnalysis::TransformCoordinates(Double_t x, Double_t y, Double_t z, Double_t seg){

        std::tuple<Double_t,Double_t, Double_t> transformedCOGs;


        // if(seg==1) {
        Double_t offsetX=-1*calsizeXY/2+tiledimX/2;
        Double_t offsetY=-1*calsizeXY/2+tiledimY/2;
        Double_t offsetZ=calsizeXY/2+AbsoThickness+GapThickness/2;
        //  std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;

        TVector3 coordinates((offsetX + x * tiledimX),                      // X
                             (offsetY + y * tiledimY),                      // Y
                             (offsetZ + z * (AbsoThickness+GapThickness))); // Z

        //coordinates.Print();

        //rotate into correct segment

        if (seg==1) {
                transformedCOGs=std::make_tuple(coordinates.X(), coordinates.Y(), coordinates.Z());
        }
        if (seg==2) {
                coordinates.RotateX(TMath::Pi());
                transformedCOGs=std::make_tuple(coordinates.X(), coordinates.Y(), coordinates.Z());
        }
        if (seg==3) {
                coordinates.RotateX(-1*TMath::Pi()/2);
                transformedCOGs=std::make_tuple(coordinates.X(), coordinates.Y(), coordinates.Z());
        }
        if (seg==4) {
                coordinates.RotateX(-1*TMath::Pi()*(3/2));
                transformedCOGs=std::make_tuple(coordinates.X(), coordinates.Y(), coordinates.Z());
        }
        if (seg==5) {
                coordinates.RotateY(TMath::Pi()/2);
                transformedCOGs=std::make_tuple(coordinates.X(), coordinates.Y(), coordinates.Z());
        }
        if (seg==6) {
                coordinates.RotateY(-1*TMath::Pi()/2);
                transformedCOGs=std::make_tuple(coordinates.X(), coordinates.Y(), coordinates.Z());
        }





        // }

        // if(seg==2) {
        //         Double_t offsetX=-1*calsizeXY/2+tiledimX/2;
        //         Double_t offsetY=calsizeXY/2-tiledimY/2;
        //         Double_t offsetZ=-1*(calsizeXY/2+AbsoThickness+GapThickness/2);
        //         //  std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;
        //
        //
        //
        //         auto tc = std::make_tuple((offsetX + x * tiledimX),                  // X
        //                                   (offsetY - y * tiledimY),                   // Y
        //                                   (offsetZ - z * (AbsoThickness+GapThickness)));                                 // Z
        //
        //         transformedCOGs=tc;
        //
        // }
        //
        // if(seg==3) {
        //         Double_t offsetX=-1*calsizeXY/2+tiledimX/2;
        //         Double_t offsetY=-1*(calsizeXY/2+AbsoThickness+GapThickness/2);
        //         Double_t offsetZ=calsizeXY/2-tiledimX/2;
        //         //  std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;
        //
        //
        //
        //         auto tc = std::make_tuple((offsetX + x * tiledimX),                  // X
        //                                   (offsetY - z * (AbsoThickness+GapThickness)),                   // Y
        //                                   (offsetZ - y * tiledimX));                                 // Z
        //         transformedCOGs=tc;
        //
        // }
        //
        // if(seg==4) {
        //         Double_t offsetX=-1*calsizeXY/2+tiledimX/2;
        //         Double_t offsetY=-1*(calsizeXY/2+AbsoThickness+GapThickness/2);
        //         Double_t offsetZ=-1*calsizeXY/2+tiledimY/2;
        //         //  std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;
        //
        //
        //         auto tc = std::make_tuple((offsetX + x * tiledimX),                  // X
        //                                   (offsetY - z * (AbsoThickness+GapThickness)),                   // Y
        //                                   (offsetZ + y * tiledimY));                                 // Z
        //         transformedCOGs=tc;
        //
        //
        // }
        //
        // if(seg==5) {
        //         Double_t offsetX=(calsizeXY/2+AbsoThickness+GapThickness/2);
        //         Double_t offsetY=-1*calsizeXY/2+tiledimX/2;
        //         Double_t offsetZ=calsizeXY/2-tiledimX/2;
        //         //  std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;
        //
        //
        //
        //         auto tc = std::make_tuple((offsetX + z * (AbsoThickness+GapThickness)),                  // X
        //                                   (offsetY + y * tiledimY),                   // Y
        //                                   (offsetZ - x * tiledimX));                                // Z
        //
        //         transformedCOGs=tc;
        //
        //
        // }
        //
        // if(seg==6) {
        //         Double_t offsetX=-1*(calsizeXY/2+AbsoThickness+GapThickness/2);
        //         Double_t offsetY=-1*calsizeXY/2+tiledimX/2;
        //         Double_t offsetZ=-1*calsizeXY/2+tiledimX/2;
        //         //  std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;
        //
        //
        //
        //         auto tc = std::make_tuple((offsetX - z * (AbsoThickness+GapThickness)),                  // X
        //                                   (offsetY + y * tiledimX),                   // Y
        //                                   (offsetZ + x * tiledimY));                                 // Z
        //
        //         transformedCOGs=tc;
        //
        //
        // }

        return transformedCOGs;
}


//-------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::FitCOGsPion(Int_t event){
        //gErrorIgnoreLevel = 1001;
        EcalTree->GetEntry(event);

        TVector3 ph1_orig=Cevent->MomentumPh1().Unit();
        TVector3 ph2_orig=Cevent->MomentumPh2().Unit();

        //Int_t seg;


        Fcn myfcn;
        myfcn.SetMode("photon3D");

        MnUserParameters upar;
        double error_minimizer_parameters = 1e-4;


        myfcn.SetCOGs(COGCollectionPH1[event]);

        upar.Add("a_1", std::get<0>(showerCOGPhoton1[event]), error_minimizer_parameters);
        upar.Add("a_2", std::get<1>(showerCOGPhoton1[event]), error_minimizer_parameters);
        upar.Add("a_3", std::get<2>(showerCOGPhoton1[event]), error_minimizer_parameters);

        // upar.Fix("a_1");
        // upar.Fix("a_2");
        // upar.Fix("a_3");

        //std::cout<<"set COGPH1"<<std::endl;
        upar.Add("v_1", ph1_orig.X(), error_minimizer_parameters);
        upar.Add("v_2", ph1_orig.Y(), error_minimizer_parameters);
        upar.Add("v_3", ph1_orig.Z(), error_minimizer_parameters);

        //  upar.SetLimits("v_3", 0.6, 1.0);
        //cout << upar << endl;

        MnMigrad migrad(myfcn, upar, 2);

        myfcn.SetCurrentEvent(event);

        //std::cout<<"event: "<<events<<std::endl;
        FunctionMinimum min = migrad();

        //std::cout<<min<<std::endl;

        MnUserParameterState uParState = min.UserState();

        //std::cout<<"attempt adding params for photon 1"<<std::endl;

        TVector3 A_1(uParState.Value("a_1"), uParState.Value("a_2"), uParState.Value("a_3"));
        TVector3 V_1(uParState.Value("v_1"), uParState.Value("v_2"), uParState.Value("v_3"));
        V_1=V_1.Unit();



        auto tp1 = std::make_pair(A_1, V_1);
        //std::cout<<"errga"<<std::endl;
        DirectionPhoton1.push_back(tp1);
        //std::cout<<"adding params for photon 1"<<std::endl;


}

void TROOTAnalysis::FitCOGsPion2(Int_t event){

        EcalTree->GetEntry(event);

        TVector3 ph1_orig=Cevent->MomentumPh1().Unit();
        TVector3 ph2_orig=Cevent->MomentumPh2().Unit();

        //photon2

        Fcn myfcn2;
        myfcn2.SetMode("photon3D");
        MnUserParameters upar2;
        double error_minimizer_parameters = 1e-4;

        myfcn2.SetCOGs(COGCollectionPH2[event]);

        upar2.Add("a_1", std::get<0>(showerCOGPhoton2[event]), error_minimizer_parameters);
        upar2.Add("a_2", std::get<1>(showerCOGPhoton2[event]), error_minimizer_parameters);
        upar2.Add("a_3", std::get<2>(showerCOGPhoton2[event]), error_minimizer_parameters);

        //
        // upar2.Fix("a_1");
        // upar2.Fix("a_2");
        // upar2.Fix("a_3");

        upar2.Add("v_1", ph2_orig.X(), error_minimizer_parameters);
        upar2.Add("v_2", ph2_orig.Y(), error_minimizer_parameters);
        upar2.Add("v_3", ph2_orig.Z(), error_minimizer_parameters);

        //  upar2.SetLimits("v_3", 0.6, 1.0);



        MnMigrad migrad2(myfcn2, upar2, 2);

        myfcn2.SetCurrentEvent(event);

        //std::cout<<"event: "<<events<<std::endl;
        FunctionMinimum min2 = migrad2();

        //std::cout<<min<<std::endl;

        MnUserParameterState uParState2 = min2.UserState();


        //std::cout<<"attempt adding params for photon 2"<<std::endl;
        TVector3 A_2(uParState2.Value("a_1"), uParState2.Value("a_2"), uParState2.Value("a_3"));

        TVector3 V_2(uParState2.Value("v_1"), uParState2.Value("v_2"), uParState2.Value("v_3"));
        V_2=V_2.Unit();


        auto tp2 = std::make_pair(A_2, V_2);
        //std::cout<<"errga"<<std::endl;
        DirectionPhoton2.push_back(tp2);
        //std::cout<<"adding params for photon 2"<<std::endl





        Double_t angle_rec=DirectionPhoton1[event].second.Angle(DirectionPhoton2[event].second);
        angle_rec=(360/(2*TMath::Pi()))*angle_rec;
        angh_rec->Fill(angle_rec);

        Double_t angle_true= ph1_orig.Angle(ph2_orig);
        angle_true=(360/(2*TMath::Pi()))*angle_true;
        angh_true->Fill(angle_true);


        dx1->Fill(ph1_orig.X()-DirectionPhoton1[event].second.X());
        dy1->Fill(ph1_orig.Y()-DirectionPhoton1[event].second.Y());
        dz1->Fill(ph1_orig.Z()-DirectionPhoton1[event].second.Z());

        dx2->Fill(ph2_orig.X()-DirectionPhoton2[event].second.X());
        dy2->Fill(ph2_orig.Y()-DirectionPhoton2[event].second.Y());
        dz2->Fill(ph2_orig.Z()-DirectionPhoton2[event].second.Z());




}

//--------------------------------------------------------------------------------------------------------------------

// void TROOTAnalysis::AnalyzePions(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent){
//         GetInvariantMass();
//         std::cout<<"PlottingChimap"<<std::endl;
//         //PlotChimap(3);
//         std::cout<<"DoneChimap"<<std::endl;
//
//
//         //MinimizeClosestApproach();
//         PionLocator();
//
//         FitParamsPions.clear();
//         ClosestApproach.clear();
// }


//--------------------------------------------------------------------------------------------------------------------------------


std::pair<TVector3, TVector3> TROOTAnalysis::FindClosestApproach(Int_t event){
        gErrorIgnoreLevel = 1001;
        Fcn myfcnPi;

        //myfcnPi.SetParamsPions(FitParamsPions[event]);
        myfcnPi.SetTrack(DirectionPhoton1[event], DirectionPhoton2[event]);
        myfcnPi.SetMode("pion");

        MnUserParameters uparPi;

        double error_minimizer_parameters = 1e-4;

        uparPi.Add("p1", 0., error_minimizer_parameters);
        uparPi.Add("p2", 0., error_minimizer_parameters);

        MnMigrad migradPi(myfcnPi, uparPi, 2);



        uparPi.SetValue("p1", 0.);
        uparPi.SetValue("p2", 0.);

        //myfcnPi.SetCurrentEvent(event);

        FunctionMinimum minPi = migradPi();

        //std::cout<<minPi<<std::endl;

        MnUserParameterState userParameterStatePi = minPi.UserState();


        TVector3 app1=DirectionPhoton1[event].first +userParameterStatePi.Value("p1")*DirectionPhoton1[event].second;
        TVector3 app2=DirectionPhoton2[event].first +userParameterStatePi.Value("p2")*DirectionPhoton2[event].second;
        //app1.Print();
        //app2.Print();

        TVector3 ca = app1+app2;
        TVector3 d=app1-app2;

        Double_t dist=d.Mag();

        // Double_t mx1=(std::get<0>(FitParamsPions[event][0])*userParameterStatePi.Value("a") + std::get<1>(FitParamsPions[event][0]));
        // Double_t mx2=(std::get<0>(FitParamsPions[event][1])*userParameterStatePi.Value("a") + std::get<1>(FitParamsPions[event][1]));
        //
        // Double_t my1=(std::get<2>(FitParamsPions[event][0])*userParameterStatePi.Value("a") + std::get<3>(FitParamsPions[event][0]));
        // Double_t my2=(std::get<2>(FitParamsPions[event][1])*userParameterStatePi.Value("a") + std::get<3>(FitParamsPions[event][1]));
        //
        // Double_t ReconstX= (mx1 + mx2)/2;
        // Double_t ReconstY= (my1 + my2)/2;
        //
        // Double_t X_=(mx1-mx2);
        // Double_t Y_=(my1-my2);
        //
        // Double_t approachdist= TMath::Sqrt(X_*X_ + Y_*Y_);
        //
        ClosestApproach.push_back(std::make_tuple(ca.X()/2,ca.Y()/2,ca.Z()/2, dist));

        TVector3 ClApp(ca.X()/2,ca.Y()/2,ca.Z()/2);
        TVector3 cadiff=ClApp - gunposition;
        delR->Fill(cadiff.Mag());

        Double_t DeltaX=gunposition.X() - ca.X()/2;
        Double_t DeltaY=gunposition.Y() - ca.Y()/2;
        Double_t DeltaZ=gunposition.Z() - ca.Z()/2;

        delx->Fill(DeltaX);
        dely->Fill(DeltaY);
        delz->Fill(DeltaZ);
        appdist1->Fill(dist);

        return std::make_pair(app1, app2);

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void TROOTAnalysis::GetInvariantMass(Int_t event){
        gErrorIgnoreLevel = 1001;
        EcalTree->GetEntry(event);

        TVector3 Ph1_direction;         //reconstructed direction unitvectors
        TVector3 Ph2_direction;

        TVector3 Ph1_dir_true;
        TVector3 Ph2_dir_true;

        TLorentzVector L_ph1;
        TLorentzVector L_ph2;
        TLorentzVector L_pi;

        TLorentzVector L_ph1_reco;
        TLorentzVector L_ph2_reco;
        TLorentzVector L_pi_reco;

        TLorentzVector L_ph1_true;
        TLorentzVector L_ph2_true;
        TLorentzVector L_pi_true;

        //True Direction reconstructed energy
        Ph1_dir_true = Cevent->MomentumPh1().Unit();            //get truedirectionfrom Tree
        Ph2_dir_true = Cevent->MomentumPh2().Unit();
        Ph1_direction = DirectionPhoton1[event].second;
        Ph2_direction = DirectionPhoton2[event].second;


        L_ph1_true.SetPxPyPzE(Ph1_dir_true.X()*EnergyPhoton1[event].first,
                              Ph1_dir_true.Y()*EnergyPhoton1[event].first,
                              Ph1_dir_true.Z()*EnergyPhoton1[event].first,
                              EnergyPhoton1[event].first);

        L_ph2_true.SetPxPyPzE(Ph2_dir_true.X()*EnergyPhoton2[event].first,
                              Ph2_dir_true.Y()*EnergyPhoton2[event].first,
                              Ph2_dir_true.Z()*EnergyPhoton2[event].first,
                              EnergyPhoton2[event].first);

        //calculate Invariant Mass
        L_pi_true=L_ph1_true + L_ph2_true;
        Double_t invMassPion_true =  L_pi_true.M();
        InvMassSim->Fill(invMassPion_true);
        //  std::cout<<"calculated invariant mass1"<<std::endl;
        //true energy reconstructed direction


        L_ph1.SetPxPyPzE(Ph1_direction.X() * Cevent->EnergyPhoton1(),
                         Ph1_direction.Y() * Cevent->EnergyPhoton1(),
                         Ph1_direction.Z() * Cevent->EnergyPhoton1(),
                         Cevent->EnergyPhoton1());
        L_ph2.SetPxPyPzE(Ph2_direction.X() * Cevent->EnergyPhoton2(),
                         Ph2_direction.Y() * Cevent->EnergyPhoton2(),
                         Ph2_direction.Z() * Cevent->EnergyPhoton2(),
                         Cevent->EnergyPhoton2());

        L_pi=L_ph1 + L_ph2;
        Double_t invMassPion=L_pi.M();
        //std::cout<<"calculated invariant mass2"<<std::endl;
        InvMassReco->Fill(invMassPion);

        InvariantMass.push_back(invMassPion);

        //reconstruchted energy reconstruchted direction

        L_ph1_reco.SetPxPyPzE(Ph1_direction.X() * EnergyPhoton1[event].first,
                              Ph1_direction.Y() * EnergyPhoton1[event].first,
                              Ph1_direction.Z() * EnergyPhoton1[event].first,
                              EnergyPhoton1[event].first);
        L_ph2_reco.SetPxPyPzE(Ph2_direction.X() * EnergyPhoton2[event].first,
                              Ph2_direction.Y() * EnergyPhoton2[event].first,
                              Ph2_direction.Z() * EnergyPhoton2[event].first,
                              EnergyPhoton2[event].first);

        L_pi_reco=L_ph1_reco+L_ph2_reco;

        Double_t invMassPion_reco=L_pi_reco.M();

        //std::cout<<"calculated invariant mass3"<<std::endl;

        InvMassAllreco->Fill(invMassPion_reco);

        //reconstruchted direction

        Double_t anglereco=(360/(2*TMath::Pi()))*Ph1_direction.Angle(Ph2_direction);
        Angle_vs_massReco->Fill(invMassPion, anglereco);

        //simulated direction

        Double_t anglesim=(360/(2*TMath::Pi()))*Ph1_dir_true.Angle(Ph2_dir_true);
        Angle_vs_massSim->Fill(invMassPion_true, anglesim);

        //all reconstructed

        Angle_vs_massAllreco->Fill(invMassPion_reco, anglereco);



}

//------------------------------------------------------------------------------------------------------------------------------------------

// void TROOTAnalysis::MinimizeClosestApproach(Int_t event, std::pair<TVector3, TVector3> p){
//
//         EcalTree->GetEntry(event);
//
//         Fcn myfcnLoc;
//         myfcnLoc.SetMode("approach");
//
//         TVector3 V_ph1_0_loc;
//         TVector3 V_ph1_1_loc;
//         TVector3 V_ph2_0_loc;
//         TVector3 V_ph2_1_loc;
//
//         MnUserParameters uparLoc;
//
//         double error_minimizer_parameters = 1e-4;
//
//         uparLoc.Add("x", std::get<0>(ClosestApproach[event]), error_minimizer_parameters);
//         uparLoc.Add("y", std::get<1>(ClosestApproach[event]), error_minimizer_parameters);
//         uparLoc.Add("z", std::get<2>(ClosestApproach[event]), error_minimizer_parameters);
//
//         TVector3 Ph1_dir_true_loc=Cevent->MomentumPh1().Unit();
//         TVector3 Ph2_dir_true_loc=Cevent->MomentumPh2().Unit();
//
//         TVector3 Ph1_direction_loc =DirectionPhoton1[event].second;
//         TVector3 Ph2_direction_loc =DirectionPhoton2[event].second;
//
//         TVector3 startvalues1(std::get<0>(ClosestApproach[event]) - std::get<0>(showerCOGPhoton1[event]),
//                               std::get<1>(ClosestApproach[event]) - std::get<1>(showerCOGPhoton1[event]),
//                               std::get<2>(ClosestApproach[event]) - std::get<2>(showerCOGPhoton1[event]));
//         Double_t deltaT1=startvalues1.Mag();
//
//         TVector3 startvalues2(std::get<0>(ClosestApproach[event]) - std::get<0>(showerCOGPhoton2[event]),
//                               std::get<1>(ClosestApproach[event]) - std::get<1>(showerCOGPhoton2[event]),
//                               std::get<2>(ClosestApproach[event]) - std::get<2>(showerCOGPhoton2[event]));
//         Double_t deltaT2=startvalues2.Mag();
//
//         myfcnLoc.SetPhotonPoints(p.first.X(), p.first.Y(), p.first.Z(), EnergyPhoton1[event].first,
//                                  p.second.X(), p.second.Y(), p.second.Z(), EnergyPhoton2[event].first);
//         myfcnLoc.SetDeltas(deltaT1, deltaT2);
//
//         myfcnLoc.SetshowerCOG(showerCOGPhoton1[event], showerCOGPhoton2[event]);
//
//         MnMigrad migradLoc(myfcnLoc, uparLoc, 2);
//
//         FunctionMinimum minLoc = migradLoc();
//
//         //std::cout<<minLoc<<std::endl;
//
//         MnUserParameterState userParameterStateLoc = minLoc.UserState();
//
//         pilocX->Fill(userParameterStateLoc.Value("x"));
//         pilocY->Fill(userParameterStateLoc.Value("y"));
//         pilocZ->Fill(userParameterStateLoc.Value("z"));
//
//         truelocX->Fill(std::get<0>(ClosestApproach[event]));
//         truelocY->Fill(std::get<1>(ClosestApproach[event]));
//         truelocZ->Fill(std::get<2>(ClosestApproach[event]));
//
// }


void TROOTAnalysis::PionLocator(Int_t event, std::pair<TVector3, TVector3> p){
        gErrorIgnoreLevel = 1001;
        EcalTree->GetEntry(event);

        Fcn myfcnLoc;
        Fcn myfcnEn;

        myfcnLoc.SetMode("approach");


        myfcnEn.SetMode("InvMass");



        MnUserParameters uparLoc;
        MnUserParameters uparEn;

        double error_minimizer_parameters = 1e-4;

        uparLoc.Add("x", std::get<0>(ClosestApproach[event]), error_minimizer_parameters);
        uparLoc.Add("y", std::get<1>(ClosestApproach[event]), error_minimizer_parameters);
        uparLoc.Add("z", std::get<2>(ClosestApproach[event]), error_minimizer_parameters);


        TVector3 Ph1_dir_true_loc=Cevent->MomentumPh1().Unit();
        TVector3 Ph2_dir_true_loc=Cevent->MomentumPh2().Unit();

        TVector3 Ph1_direction_loc = DirectionPhoton1[event].second;

        TLorentzVector L_ph1_loc;

        L_ph1_loc.SetPxPyPzE(Ph1_direction_loc.X() * EnergyPhoton1[event].first,
                             Ph1_direction_loc.Y() * EnergyPhoton1[event].first,
                             Ph1_direction_loc.Z() * EnergyPhoton1[event].first,
                             EnergyPhoton1[event].first);

        TVector3 Ph2_direction_loc = DirectionPhoton2[event].second;

        TLorentzVector L_ph2_loc;

        L_ph2_loc.SetPxPyPzE(Ph2_direction_loc.X() * EnergyPhoton2[event].first,
                             Ph2_direction_loc.Y() * EnergyPhoton2[event].first,
                             Ph2_direction_loc.Z() * EnergyPhoton2[event].first,
                             EnergyPhoton2[event].first);


        TVector3 startvalues1(std::get<0>(ClosestApproach[event]) - std::get<0>(showerCOGPhoton1[event]),
                              std::get<1>(ClosestApproach[event]) - std::get<1>(showerCOGPhoton1[event]),
                              std::get<2>(ClosestApproach[event]) - std::get<2>(showerCOGPhoton1[event]));
        Double_t deltaT1=startvalues1.Mag();

        TVector3 startvalues2(std::get<0>(ClosestApproach[event]) - std::get<0>(showerCOGPhoton2[event]),
                              std::get<1>(ClosestApproach[event]) - std::get<1>(showerCOGPhoton2[event]),
                              std::get<2>(ClosestApproach[event]) - std::get<2>(showerCOGPhoton2[event]));
        Double_t deltaT2=startvalues2.Mag();




        myfcnLoc.SetPhotonPoints(p.first.X(), p.first.Y(), p.first.Z(), EnergyPhoton1[event].first,
                                 p.second.X(), p.second.Y(), p.second.Z(), EnergyPhoton2[event].first);

        myfcnLoc.SetDeltas(deltaT1, deltaT2);

        myfcnLoc.SetshowerCOG(showerCOGPhoton1[event], showerCOGPhoton2[event]);

        MnMigrad migradLoc(myfcnLoc, uparLoc, 2);

        FunctionMinimum minLoc = migradLoc();

        //std::cout<<minLoc<<std::endl;

        MnUserParameterState userParameterStateLoc = minLoc.UserState();





        inputX->Fill(std::get<0>(ClosestApproach[event]));
        inputY->Fill(std::get<1>(ClosestApproach[event]));
        inputZ->Fill(std::get<2>(ClosestApproach[event]));

        pilocX->Fill(userParameterStateLoc.Value("x"));
        pilocY->Fill(userParameterStateLoc.Value("y"));
        pilocZ->Fill(userParameterStateLoc.Value("z"));





        TVector3 fittedCA(userParameterStateLoc.Value("x"),userParameterStateLoc.Value("y"),userParameterStateLoc.Value("z"));
        TVector3 diff=fittedCA-Cevent->GunPos();
        Double_t deltaR= diff.Mag();
        delR_fit->Fill(deltaR);

        Double_t error_minimizer_parameters_En=0.0001;

        uparEn.Add("x", userParameterStateLoc.Value("x"), error_minimizer_parameters_En);
        uparEn.Add("y", userParameterStateLoc.Value("y"), error_minimizer_parameters_En);
        uparEn.Add("z", userParameterStateLoc.Value("z"), error_minimizer_parameters_En);

        //uparEn.SetLimits("x",userParameterStateLoc.Value("x")-20, userParameterStateLoc.Value("x")+20);
        //uparEn.SetLimits("x",userParameterStateLoc.Value("y")-30, userParameterStateLoc.Value("y")+30);

        myfcnEn.SetDeltas(deltaT1, deltaT2);

        myfcnLoc.SetPhotonPoints(p.first.X(), p.first.Y(), p.first.Z(), EnergyPhoton1[event].first,
                                 p.second.X(), p.second.Y(), p.second.Z(), EnergyPhoton2[event].first);

        myfcnEn.SetshowerCOG(showerCOGPhoton1[event], showerCOGPhoton2[event]);

        myfcnEn.SetClosestApproach(userParameterStateLoc.Value("x"),userParameterStateLoc.Value("y"),userParameterStateLoc.Value("z"));

        MnMigrad migradEn(myfcnEn, uparEn, 2);

        FunctionMinimum minEn = migradEn();

        MnUserParameterState userParameterStateEn=minEn.UserState();


        invmasFitX->Fill(userParameterStateEn.Value("x"));
        invmasFitY->Fill(userParameterStateEn.Value("y"));
        invmasFitZ->Fill(userParameterStateEn.Value("z"));










        // TCanvas * sgd = new TCanvas("sgdC", "sgdC");
        // TH1D * lkl = new TH1D("lkl", "lkl", 1000,0,1000);
        //
        // Double_t x[1000];
        // Double_t y[1000];
        // for(Int_t j=0; j<nofEntries; j++) {
        //
        //         EcalTree->GetEntry(j);
        //
        //         TVector3 origdir1=Cevent->MomentumPh1().Unit();
        //         TVector3 origdir2=Cevent->MomentumPh2().Unit();
        //
        //
        //
        //         TVector3 dir_ph1(std::get<0>(showerCOGPhoton1[j])-std::get<0>(ClosestApproach[j]),
        //                          std::get<1>(showerCOGPhoton1[j])-std::get<1>(ClosestApproach[j]),
        //                          std::get<2>(showerCOGPhoton1[j])-std::get<2>(ClosestApproach[j]));
        //         dir_ph1=dir_ph1.Unit();
        //
        //         TLorentzVector lv_ph1;
        //         lv_ph1.SetPxPyPzE(dir_ph1.X()*Cevent->EnergyPhoton1(),dir_ph1.Y()*Cevent->EnergyPhoton1(),dir_ph1.Z()*Cevent->EnergyPhoton1(),Cevent->EnergyPhoton1());
        //
        //         TVector3 dir_ph2(std::get<0>(showerCOGPhoton2[j])-std::get<0>(ClosestApproach[j]),
        //                          std::get<1>(showerCOGPhoton2[j])-std::get<1>(ClosestApproach[j]),
        //                          std::get<2>(showerCOGPhoton2[j])-std::get<2>(ClosestApproach[j]));
        //         dir_ph2=dir_ph2.Unit();
        //
        //         TLorentzVector lv_ph2;
        //         lv_ph2.SetPxPyPzE(dir_ph2.X()*Cevent->EnergyPhoton2(),dir_ph2.Y()*Cevent->EnergyPhoton2(),dir_ph2.Z()*Cevent->EnergyPhoton2(),Cevent->EnergyPhoton2());
        //
        //         lv_ph1.Print();
        //         lv_ph2.Print();
        //         TLorentzVector lv_pi;
        //         lv_pi  = lv_ph1 + lv_ph2;
        //
        //         lv_pi.Print();
        //
        //         Double_t curr_invmass=lv_pi.M();
        //         std::cout<<curr_invmass<<std::endl;
        //         lkl->Fill(curr_invmass);
        //
        //
        //         // origdir1.Print();
        //         // dir_ph1.Print();
        //         //
        //         // std::cout<<"PH2____________________________"<<std::endl;
        //         //
        //         // origdir2.Print();
        //         // dir_ph2.Print();
        //
        //         std::cout<<"PH1____________________________"<<std::endl;
        //
        // }
        // // TGraph* g1= new TGraph(1000, x,y);
        // sgd->cd(0);
        // lkl->Draw();

}

void TROOTAnalysis::DrawHists(){


        dx1->GetXaxis()->SetTitle("Delta X[mm]");
        dx1->GetYaxis()->SetTitle("Entries");

        dx2->GetXaxis()->SetTitle("Delta X[mm]");
        dx2->GetYaxis()->SetTitle("Entries");

        dy1->GetXaxis()->SetTitle("Delta Y[mm]");
        dy1->GetYaxis()->SetTitle("Entries");

        dy2->GetXaxis()->SetTitle("Delta Y[mm]");
        dy2->GetYaxis()->SetTitle("Entries");

        dz1->GetXaxis()->SetTitle("Delta Z[mm]");
        dz1->GetYaxis()->SetTitle("Entries");

        dz2->GetXaxis()->SetTitle("Delta Z[mm]");
        dz2->GetYaxis()->SetTitle("Entries");

        D->Divide(2,3,0.01,0.01);
        gStyle->SetOptStat(2222);
        D->cd(1);
        dx1->Draw();
        D->cd(3);
        dy1->Draw();
        D->cd(5);
        dz1->Draw();
        D->cd(2);
        dx2->Draw();
        D->cd(4);
        dy2->Draw();
        D->cd(6);
        dz2->Draw();

        Delta1->Divide(2,2,0.01,0.01);

        gStyle->SetOptStat(2222);

        delx->GetXaxis()->SetTitle("DeltaX[mm]");
        dely->GetXaxis()->SetTitle("DeltaY[mm]");
        delz->GetXaxis()->SetTitle("DeltaZ[mm]");

        delx->GetYaxis()->SetTitle("Entries");
        dely->GetYaxis()->SetTitle("Entries");
        delz->GetYaxis()->SetTitle("Entries");
        appdist1->GetXaxis()->SetTitle("Distance of closest Approach [mm]");
        appdist1->GetYaxis()->SetTitle("Entries");

        Delta1->cd(1);
        delx->Draw();
        Delta1->cd(2);
        dely->Draw();
        Delta1->cd(3);
        delz->Draw();
        Delta1->cd(4);
        appdist1->Draw();

        // Delta2->cd(0);
        // delx->Draw();
        //
        // Delta3->cd(0);
        // appdist1->Draw();
        //
        // Delta4->cd(0);
        // delz->Draw();


        // histograms for different invariaont mass calculations

        InvMassReco->GetXaxis()->SetTitle("InvariantMass[MeV]");
        InvMassReco->GetYaxis()->SetTitle("Entries");

        InvMassSim->GetXaxis()->SetTitle("InvariantMass[MeV]");
        InvMassSim->GetYaxis()->SetTitle("Entries");

        InvMassAllreco->GetXaxis()->SetTitle("InvariantMass[MeV]");
        InvMassAllreco->GetYaxis()->SetTitle("Entries");

        //histograms for Angle_vs_InvariantMass plots

        Angle_vs_massReco->GetXaxis()->SetTitle("Invariant Mass [MeV]");
        Angle_vs_massReco->GetYaxis()->SetTitle("Angle[deg]");
        Angle_vs_massReco->GetYaxis()->SetTitleOffset(1.4);

        Angle_vs_massSim->GetXaxis()->SetTitle("Invariant Mass [MeV]");
        Angle_vs_massSim->GetYaxis()->SetTitle("Angle[deg]");
        Angle_vs_massSim->GetYaxis()->SetTitleOffset(1.4);


        Angle_vs_massAllreco->GetXaxis()->SetTitle("Invariant Mass [MeV]");
        Angle_vs_massAllreco->GetYaxis()->SetTitle("Angle[deg]");
        Angle_vs_massAllreco->GetYaxis()->SetTitleOffset(1.4);

        IM1->Divide(2,2,0.01,0.01);

        gStyle->SetOptStat(2222);
        // IM1->cd(1);
        // InvMassReco->Draw();

        IM1->cd(1);
        InvMassSim->Draw();

        IM1->cd(2);
        InvMassAllreco->Draw();



        Int_t palette[5];
        palette[0] = 3;
        palette[1] = 4;
        palette[2] = 5;
        palette[3] = 6;
        palette[4] = 2;
        gStyle->SetPalette(5, palette);
        gStyle->SetOptStat(0);

        // IM1->cd(4);
        // Angle_vs_massReco->GetZaxis()->SetNdivisions(5);
        // Angle_vs_massReco->Draw("colz");
        // Angle_vs_massReco->SetMaximum(5);
        // Angle_vs_massReco->SetMinimum(0);

        IM1->cd(3);
        Angle_vs_massSim->GetZaxis()->SetNdivisions(5);
        Angle_vs_massSim->Draw("colz");
        Angle_vs_massSim->SetMaximum(5);
        Angle_vs_massSim->SetMinimum(0);

        IM1->cd(4);
        Angle_vs_massAllreco->GetZaxis()->SetNdivisions(5);
        Angle_vs_massAllreco->Draw("colz");
        Angle_vs_massAllreco->SetMaximum(5);
        Angle_vs_massAllreco->SetMinimum(0);
        gStyle->SetOptStat(2222);
        PiLoc->Divide(3,3,0.01,0.01);

        pilocX->GetXaxis()->SetTitle("Closest Approach Position X [mm]");
        pilocX->GetYaxis()->SetTitle("Entries");
        PiLoc->cd(4);
        pilocX->Draw();

        pilocY->GetXaxis()->SetTitle("Closest Approach Position Y [mm]");
        pilocY->GetYaxis()->SetTitle("Entries");
        PiLoc->cd(5);
        pilocY->Draw();

        pilocZ->GetXaxis()->SetTitle("Closest Approach Position Z [mm]");
        pilocZ->GetYaxis()->SetTitle("Entries");
        PiLoc->cd(6);
        pilocZ->Draw();

        inputX->GetXaxis()->SetTitle("Closest Approach Position X [mm]");
        inputX->GetYaxis()->SetTitle("Entries");
        PiLoc->cd(1);
        inputX->Draw();

        inputY->GetXaxis()->SetTitle("Closest Approach Position Y [mm]");
        inputY->GetYaxis()->SetTitle("Entries");
        PiLoc->cd(2);
        inputY->Draw();

        inputZ->GetXaxis()->SetTitle("Closest Approach Position Z [mm]");
        inputZ->GetYaxis()->SetTitle("Entries");

        PiLoc->cd(3);
        inputZ->Draw();

        // PiLoc1->cd(0);
        // inputZ->Draw();
        //
        // PiLoc2->cd(0);
        // pilocZ->Draw();

        PiLoc->cd(7);
        invmasFitX->Draw();

        PiLoc->cd(8);
        invmasFitY->Draw();

        PiLoc->cd(9);
        invmasFitZ->Draw();

        delR->GetXaxis()->SetTitle("Distance of true to reconstructed vertex [mm]");
        delR->GetYaxis()->SetTitle("Entries");

        DeltaR->cd(0);
        delR->Draw();
        Double_t x, q;
        q = 0.5; // 0.5 for "median"
        delR->ComputeIntegral(); // just a precaution
        delR->GetQuantiles(1, &x, &q);
        std::cout << "median = " << x << std::endl;
        DeltaR->Update();

        TLine *line = new TLine(x,0,x,16);
        line->SetLineColor(kRed);
        line->SetLineWidth(3);
        line->Draw();

        DeltaR_fit->cd(0);
        delR_fit->Draw();
        delR_fit->ComputeIntegral(); // just a precaution
        delR_fit->GetQuantiles(1, &x, &q);
        std::cout << "median = " << x << std::endl;
        DeltaR_fit->Update();

        TLine *line_fit = new TLine(x,0,x,16);
        line_fit->SetLineColor(kRed);
        line_fit->SetLineWidth(3);
        line_fit->Draw();

        std::string path="/home/iwsatlas1/emberger/Presentations/CALICE_japan/4Pi/X_700mm_400MeV/";

        std::string path1=path+"DeltasDirection.pdf";
        D->Print(path1.c_str());
        path1=path+"DeltasDirection.C";
        D->Print(path1.c_str());

        std::string path2=path+"Invariant_mass.pdf";
        IM1->Print(path2.c_str());
        path2=path+"Invariant_mass.C";
        IM1->Print(path2.c_str());

        std::string path3=path+"Pi0Vertex.pdf";
        Delta1->Print(path3.c_str());
        path3=path+"Pi0Vertex.C";
        Delta1->Print(path3.c_str());

        std::string path4=path+"PionLocator.pdf";
        PiLoc->Print(path4.c_str());
        path4=path+"PionLocator.C";
        PiLoc->Print(path4.c_str());


        // path="/home/iwsatlas1/emberger/Presentations/CALICE_japan/Finals/";
        //
        // std::string path5=path+"X_VertexREco.pdf";
        // Delta2->Print(path5.c_str());
        // std::string path6=path+"Dist_VertexREco.pdf";
        // Delta3->Print(path6.c_str());
        // std::string path7=path+"Z_VertexREco.pdf";
        // Delta4->Print(path7.c_str());
        //
        // std::string path8=path+"PionLocator_inputZ.pdf";
        // PiLoc1->Print(path8.c_str());
        // std::string path9=path+"PionLocatorZ.pdf";
        // PiLoc2->Print(path9.c_str());
        //
        // std::string path10=path+"3DDist.pdf";
        // DeltaR->Print(path10.c_str());
        //
        // std::string path11=path+"3DDist_fit.pdf";
        // DeltaR_fit->Print(path11.c_str());


        ang1->Divide(1,2,0.01,0.01);
        ang1->cd(1);
        angh_true->Draw();

        ang1->cd(2);
        angh_rec->Draw();

        totalEnergyC->cd(0);
        totalEnergyH->GetXaxis()->SetTitle("GeV");
        totalEnergyH->GetYaxis()->SetTitle("Entries");

        totalEnergyH->Draw();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

void TROOTAnalysis::PlotChimap(Int_t event){

        // TCanvas * ChiMap=new TCanvas("ChiMapC","Chimap");
        // ChiMap->Divide(2,1,0.01,0.01);
        // TH3D * chimap = new TH3D("chimap", "chimap",1000,-500,500,1000,-500,500,1000,-1600,400);
        //
        //
        //
        // for(Int_t i=0; i<nofEntries; i++) {
        //
        //         EcalTree->GetEntry(i);
        //
        //         // Double_t e1=Cevent->EnergyPhoton1();
        //         // Double_t e2=Cevent->EnergyPhoton2();
        //         Double_t e1=EnergyPhoton1[i];
        //         Double_t e2=EnergyPhoton2[i];
        //         TVector3 startvalues1(std::get<0>(ClosestApproach[i]) - std::get<0>(showerCOGPhoton1[i]),
        //                               std::get<1>(ClosestApproach[i]) - std::get<1>(showerCOGPhoton1[i]),
        //                               std::get<2>(ClosestApproach[i]) - std::get<2>(showerCOGPhoton1[i]));
        //         Double_t deltaT1=startvalues1.Mag();
        //
        //         TVector3 startvalues2(std::get<0>(ClosestApproach[i]) - std::get<0>(showerCOGPhoton2[i]),
        //                               std::get<1>(ClosestApproach[i]) - std::get<1>(showerCOGPhoton2[i]),
        //                               std::get<2>(ClosestApproach[i]) - std::get<2>(showerCOGPhoton2[i]));
        //         Double_t deltaT2=startvalues2.Mag();
        //
        //         TVector3 PointPh1(std::get<0>(FitParamsPions[i][0]) * std::get<2>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][0]),
        //                           std::get<2>(FitParamsPions[i][0]) * std::get<2>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][0]),
        //                           std::get<2>(ClosestApproach[i]));
        //
        //         TVector3 PointPh2(std::get<0>(FitParamsPions[i][1]) * std::get<2>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][1]),
        //                           std::get<2>(FitParamsPions[i][1]) * std::get<2>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][1]),
        //                           std::get<2>(ClosestApproach[i]));
        //
        //         PointPh1.Print();
        //         PointPh2.Print();
        //         std::cout<<std::get<0>(ClosestApproach[i])<<":"<<std::get<1>(ClosestApproach[i])<<":"<<std::get<2>(ClosestApproach[i])<<std::endl;
        //
        //         Double_t chisq;
        //
        //         for(Int_t x = -500; x<501; x+=2) {
        //                 for(Int_t y = -500; y<501; y+=2) {
        //                         for(Int_t z = -1600; z< 401; z+=2) {
        //
        //                                 Double_t sigPh1=(1/TMath::Sqrt(e1))*deltaT1;
        //                                 Double_t sigPh2=(1/TMath::Sqrt(e2))*deltaT2;
        //
        //
        //                                 Double_t delR1=TMath::Sqrt(((PointPh1.X()-x)*(PointPh1.X()-x)) + ((PointPh1.Y()-y)*(PointPh1.Y()-y)) + ((PointPh1.Z()-z)*(PointPh1.Z()-z)));
        //                                 Double_t delR2=TMath::Sqrt(((PointPh2.X()-x)*(PointPh2.X()-x)) + ((PointPh2.Y()-y)*(PointPh2.Y()-y)) + ((PointPh2.Z()-z)*(PointPh2.Z()-z)));
        //
        //                                 TVector3 Dir_ph1(x-std::get<0>(showerCOGPhoton1[i]),y-std::get<1>(showerCOGPhoton1[i]),z-std::get<2>(showerCOGPhoton1[i]));
        //
        //                                 Dir_ph1=Dir_ph1.Unit();
        //
        //                                 TLorentzVector Lv_ph1(Dir_ph1.X()*e1,Dir_ph1.Y()*e1,Dir_ph1.Z()*e1,e1);
        //
        //
        //                                 TVector3 Dir_ph2(x-std::get<0>(showerCOGPhoton2[i]),y-std::get<1>(showerCOGPhoton2[i]),z-std::get<2>(showerCOGPhoton2[i]));
        //
        //                                 Dir_ph2=Dir_ph2.Unit();
        //
        //                                 TLorentzVector Lv_ph2(Dir_ph2.X()*e2,Dir_ph2.Y()*e2,Dir_ph2.Z()*e2,e2);
        //
        //                                 TLorentzVector Lv_pi=Lv_ph1+Lv_ph2;
        //
        //                                 Double_t Curr_invmass=Lv_pi.M();
        //
        //
        //
        //
        //                                 chisq=((delR1/sigPh1)*(delR1/sigPh1))+((delR2/sigPh2)*(delR2/sigPh2))+(((Curr_invmass-134.9766)/5)*((Curr_invmass-134.9766)/5));
        //
        //                                 chimap->Fill(x,y,z,chisq);
        //
        //
        //                         }
        //                         //std::cout<<y<<std::endl;
        //                 }
        //
        //         }
        //         std::cout<<"projection"<<std::endl;
        //         ChiMap->cd(1);
        //         chimap->Project3D("xz")->Draw("colz");
        //         std::cout<<"projection"<<std::endl;
        //         ChiMap->cd(2);
        //         chimap->Project3D("yz")->Draw("colz");
        //
        //         std::string ChiMapPath="Chi2Mapz/ChiMap" +std::to_string(Cevent->GapEnergy())+"_" +std::to_string(Cevent->GunPos().Z())+"_"+ std::to_string(i) + ".pdf";
        //
        //         ChiMap->Print(ChiMapPath.c_str());
        //
        //         ChiMapPath="Chi2Mapz/ChiMap" +std::to_string(Cevent->GapEnergy()) +"_" +std::to_string(Cevent->GunPos().Z())+"_"+ std::to_string(i) + ".C";
        //
        //         ChiMap->Print(ChiMapPath.c_str());
        //
        //         chimap->Reset();
        //
        //
        //
        //
        // }

}



// void TROOTAnalysis::CalcCOG(Int_t minevent, Int_t maxevent){            //calculate vector of (X,Y) tuples containing layerwise center of gravity
//
//         //variables for clustering
//         Int_t xdir=1; // integers for direction
//         Int_t ydir=0;
//         Int_t buf;
//         Int_t stepstodo=1; // how many steps to go before rotation
//         Int_t stepsctr=0; // number of steps since last rotation
//
//         Int_t currX=0; // starting the spiral on the histogram maximum
//         Int_t currY=0;
//
//         Double_t esum=0;
//         Double_t curre=0;
//
//         Int_t maxbin=0;
//         Double_t maxdep=0;
//         Int_t binx, biny, binz;
//
//         // Variables for fit
//         Double_t xerr=0;
//         Double_t yerr=0;
//         Double_t cgx=0;
//         Double_t cgy=0;
//         Double_t cgz=0;
//         Double_t Eweight=0;
//
//         for(Int_t eventstodo=minevent; eventstodo<maxevent; eventstodo++) {    //loop over all simulated events
//                 std::cout<<"CenterEvent: "<<eventstodo<<std::endl;
//                 std::cout<<eventstodo-minevent<<std::endl;
//                 EcalTree->GetEntry(eventstodo); //grab evetn from tree
//                 Eges = Cevent->GapEnergy();
//                 Int_t cnh = Cevent->NHits();
//                 Double_t integSum=0;
//                 Double_t integral;
//                 Int_t i = showerCenters[eventstodo-minevent]-7;
//                 std::cout<<i<<std::endl;
//                 Int_t collected=0;
//                 if(i<0) {
//                         i=0;
//                         collected=4;
//                 }
//
//
//                 while(collected<7 && i<50) { //loop over layers until some layers after showerstart
//                         if(i==49 && collected==0) {i=0; }
//                         std::cout<<"i: "<<i<<std::endl;
//                         for(Int_t j =0; j<cnh; j++) { //loop over all hits in laver i
//                                 if(Cevent->Hit(j)->Z()==i)
//                                         h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit()); //fill histogram with hits of layer i
//                         }
//
//                         integral=h1->Integral();
//
//                         maxbin=h1->GetMaximumBin();
//                         maxdep=h1->GetBinContent(maxbin); // get coordinates of bin with maximum energy deposition
//                         h1->GetBinXYZ(maxbin, binx, biny, binz);
//
//                         //Do spiral for clustering only if energy in Layer
//
//                         if(integral!=0) {     //check if layer containes energy
//                                 collected++;
//                                 std::cout<<collected<<":"<<i<<std::endl;
//                                 xdir=1;        // integers for direction
//                                 ydir=0;
//                                 stepstodo=1;   // how many steps to go before rotation
//                                 stepsctr=0;   // number of steps since last rotation
//
//                                 currX=binx;    // starting the spiral on the histogram maximum
//                                 currY=biny;
//
//                                 esum=0;
//                                 curre=0;
//                                 esum+=maxdep;
//                                 h2->SetBinContent(currX, currY, maxdep);
//                                 auto tp=std::make_tuple(currX, currY, maxdep);
//                                 ClusteredHits.push_back(tp);
//                                 h1->SetBinContent(binx, biny, 0);
//
//                                 while(esum<integral*0.9 || currX<binx+4) {
//                                         //do spiral until desired energyfraction is reached
//                                         currX+=xdir;
//                                         currY+=ydir;
//                                         curre=h1->GetBinContent(currX, currY);
//
//                                         if(curre>0) {     //save tile only if it containes energy
//                                                 tp=std::make_tuple(currX, currY, curre);
//                                                 ClusteredHits.push_back(tp);
//                                                 h2->SetBinContent(currX, currY, curre);
//                                         }
//                                         h1->SetBinContent(currX, currY, 0);
//                                         esum+=curre;
//                                         stepsctr++;
//                                         if(stepsctr==stepstodo) {
//                                                 stepsctr=0;
//                                                 buf=xdir; //rotate 90 deg
//                                                 xdir= -ydir;
//                                                 ydir=buf;
//
//                                                 if(ydir==0) { //incremant steps at every iteration along x
//                                                         stepstodo++;
//                                                 }
//                                         }
//                                 }
//
//                                 Double_t e=h2->Integral(); //get energy contained in cluster
//
//                                 Eweight=e/Eges; //calculate weight
//
//                                 cgx=h2->GetMean(1);
//                                 xerr=h2->GetMeanError(1);
//                                 //extract center of gravity and error
//                                 cgy=h2->GetMean(2);
//                                 yerr=h2->GetMeanError(2);
//
//                                 if(yerr<0.00001) {yerr=1/TMath::Sqrt(12); }
//                                 if(xerr<0.00001) {xerr=1/TMath::Sqrt(12); }
//
//                                 cgz=i;
//
//                                 er1->Fill(i, xerr/nofEntries);
//                                 er2->Fill(i, yerr/nofEntries);
//
//                                 //h3->Fill(cgx, cgy, cgz);
//
//                                 auto cg=std::make_tuple(cgx,cgy,cgz, xerr, yerr, Eweight);
//                                 coglist.push_back(cg);
//
//                         } // end of clustering
//
//                         integSum += integral;
//
//                         ClusteredHits.clear();
//
//                         h1->Reset();
//
//                         h2->Reset();
//
//                         i++;
//                 }
//                 COGCollection.push_back(coglist);
//                 coglist.clear();
//         }
//
// }


//-----------------------------------------------------------------------------------------------------------------------------------------------------



// void TROOTAnalysis::PrintFitHists(Int_t minevent, Int_t maxevent){
//         Int_t ent=maxevent-minevent;
//         TCanvas * c1 = new TCanvas("Fit", "Fit", 1500,1300);
//
//         c1->Divide(2,2,0.01, 0.01);
//
//         //gStyle->SetOptFit(111111111111);
//
//
//         for(Int_t i=0; i<ent; i++) {
//                 fitX->Fill(std::get<1>(FitParamsGamma[i]));
//         }
//
//         for(Int_t i=0; i<ent; i++) {
//                 fitY->Fill(std::get<3>(FitParamsGamma[i]));
//         }
//
//         for(Int_t i=0; i<ent; i++) {
//                 fitSX->Fill(std::get<0>(FitParamsGamma[i]));
//         }
//
//         for(Int_t i=0; i<ent; i++) {
//                 fitSY->Fill(std::get<2>(FitParamsGamma[i]));
//         }
//
//         c1->cd(1);
//         fitX->GetXaxis()->SetTitle("X[mm]");
//         fitX->GetYaxis()->SetTitle("Counts");
//         gStyle->SetOptStat(2222);
//         // Set stat options
//         gStyle->SetStatY(0.9);
//         // Set y-position (fraction of pad size)
//         gStyle->SetStatX(0.9);
//         // Set x-position (fraction of pad size)
//         gStyle->SetStatW(0.25);
//         // Set width of stat-box (fraction of pad size)
//         gStyle->SetStatH(0.25);
//         // Set height of stat-box (fraction of pad size)
//         fitX->Draw();
//
//         c1->cd(2);
//         fitY->GetXaxis()->SetTitle("Y[mm]");
//         fitY->GetYaxis()->SetTitle("Counts");
//         fitY->Draw();
//
//         c1->cd(3);
//         fitSX->GetXaxis()->SetTitle("X Slope");
//         fitSX->GetYaxis()->SetTitle("Counts");
//         fitSX->Draw();
//
//         c1->cd(4);
//         fitSY->GetXaxis()->SetTitle("Y Slope");
//         fitSY->GetYaxis()->SetTitle("Counts");
//         fitSY->Draw();
//
//         //c1->cd(5);
//         //er1->Draw();
//
//         //c1->cd(6);
//         //er2->Draw();
//
//         if(pathset) {
//                 //
//                 std::string imgname= "Fit.png";
//                 std::string imgpath = savepath + '/' + imgname;
//                 // TImage * img2 =TImage::Create();
//                 // img2->FromPad(c1);
//                 // img2->WriteImage(imgpath.c_str());
//                 // delete img2;
//
//                 c1->Print(imgpath.c_str());
//
//                 imgname= "Fit.C";
//                 imgpath = savepath + '/' + imgname;
//
//                 c1->Print(imgpath.c_str());
//         }
// //  er1->Reset();
// //  er2->Reset();
// }

//-----------------------------------------------------------------------------------------------------------------------------------------------------


// void TROOTAnalysis::CleanCOGs(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent){
//
//
//
//         TCanvas * c3 = new TCanvas("CleanCOGs", "CleanCOGs");
//         TH3D * clCOG = new TH3D("COGs after cleaning", "COGs after celaning",20,40,60,20,40,60,50,0,50);
//
//         Double_t MoliereRaduis=4.87; // in cm
//
//         //CalcCOG(0,10, minevent-1, maxevent-1);
//         //FitCOGs(minevent-1, maxevent-1, 5);
//
//         COGCollection.clear();
//         coglist.clear();
//         ClusteredHits.clear();
//
//         //PrintFitHists();
//         //CalcCOG(minlayer-1,maxlayer-1,minevent-1, maxevent-1 );
//
//         Int_t erasecounter=0;
//         Int_t totalcounter=0;
//         for(int i=0; i<nofEntries; i++) {                               //Delete the COGs, outside of one MoliereRadius
//                 for(int z=0; z<COGCollection[i].size(); z++) {
//
//                         Double_t distx=(std::get<0>(FitParamsGamma[i])*z+std::get<1>(FitParamsGamma[i]))-std::get<0>(COGCollection[i][z]);
//                         Double_t disty=(std::get<2>(FitParamsGamma[i])*z+std::get<3>(FitParamsGamma[i]))-std::get<1>(COGCollection[i][z]);
//                         totalcounter++;
//                         if(TMath::Sqrt(distx*distx)>MoliereRaduis || TMath::Sqrt(disty*disty)>MoliereRaduis) {
//                                 COGCollection[i].erase(COGCollection[i].begin()+z);
//                                 //std::cout<<"X distance: "<<distx<<" Y distance: "<<disty<<std::endl;
//                                 erasecounter++;
//                         }
//
//                 }
//
//         }
//
//         for(int i=0; i<nofEntries; i++) {                               //Delete the COGs, outside of one MoliereRadius
//                 for(int z=0; z<COGCollection[i].size(); z++) {
//                         clCOG->Fill(std::get<0>(COGCollection[i][z]), std::get<1>(COGCollection[i][z]), std::get<2>(COGCollection[i][z]));
//                 }
//         }
//
//         clCOG->GetXaxis()->SetTitle("X");
//         clCOG->GetYaxis()->SetTitle("Y");
//         clCOG->GetZaxis()->SetTitle("Z");
//         clCOG->SetMarkerStyle(4);
//         c3->cd();
//         clCOG->Draw("");
//
//         FitParamsGamma.clear();
//
//
//         //FitCOGs(minevent, maxevent, 5);          //refit the leftover COGs
//         //PrintFitHists(minevent, maxevent);
//         //plotCOGs();
//
//         std::cout<<"Discarded COGs: "<<erasecounter<<" Total COGs:"<<totalcounter<<std::endl;
//
//
//
//
// }

//////////////////////////////Old Section///////////////////////////////////

// void TROOTAnalysis::CalcCOG(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent){            //calculate vector of (X,Y) tuples containing layerwise center of gravity
//
//         //variables for clustering
//         Int_t xdir=1; // integers for direction
//         Int_t ydir=0;
//         Int_t buf;
//         Int_t stepstodo=1; // how many steps to go before rotation
//         Int_t stepsctr=0; // number of steps since last rotation
//
//         Int_t currX=0; // starting the spiral on the histogram maximum
//         Int_t currY=0;
//
//         Double_t esum=0;
//         Double_t curre=0;
//
//         Int_t maxbin=0;
//         Double_t maxdep=0;
//         Int_t binx, biny, binz;
//
//         // Variables for fit
//         Double_t xerr=0;
//         Double_t yerr=0;
//         Double_t cgx=0;
//         Double_t cgy=0;
//         Double_t cgz=0;
//         Double_t Eweight=0;
//
//         for(Int_t eventstodo=minevent; eventstodo<maxevent; eventstodo++) {    //loop over all simulated events
//
//                 EcalTree->GetEntry(eventstodo);       //grab evetn from tree
//                 Eges = Cevent->GapEnergy();
//                 Int_t cnh = Cevent->NHits();
//                 Double_t integral;
//                 Bool_t foundstart=false;
//                 for(Int_t i=minlayer; i<maxlayer; i++) { //loop over all layers in event
//                         for(Int_t j =0; j<cnh; j++) { //loop over all hits in laver i
//                                 if(Cevent->Hit(j)->Z()==i)
//                                         h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit());
//                                 //fill histogram with hits of layer i
//                         }
//
//                         integral=h1->Integral();
//
//                         maxbin=h1->GetMaximumBin();
//                         maxdep=h1->GetBinContent(maxbin); // get coordinates of bin with maximum energy deposition
//                         h1->GetBinXYZ(maxbin, binx, biny, binz);
//
//                         //Do spiral for clustering only if energy in Layer
//
//                         if(integral!=0) {     //check if layer containes energy
//
//                                 xdir=1;        // integers for direction
//                                 ydir=0;
//                                 stepstodo=1;   // how many steps to go before rotation
//                                 stepsctr=0;   // number of steps since last rotation
//
//                                 currX=binx;    // starting the spiral on the histogram maximum
//                                 currY=biny;
//
//                                 esum=0;
//                                 curre=0;
//                                 esum+=maxdep;
//                                 h2->SetBinContent(currX, currY, maxdep);
//                                 auto tp=std::make_tuple(currX, currY, maxdep);
//                                 ClusteredHits.push_back(tp);
//                                 h1->SetBinContent(binx, biny, 0);
//
//                                 while(esum<integral*0.9 || currX<binx+4) {
//                                         //do spiral until desired energyfraction is reached
//                                         currX+=xdir;
//                                         currY+=ydir;
//                                         curre=h1->GetBinContent(currX, currY);
//
//                                         if(curre>0) {     //save tile only if it containes energy
//                                                 tp=std::make_tuple(currX, currY, curre);
//                                                 ClusteredHits.push_back(tp);
//                                                 h2->SetBinContent(currX, currY, curre);
//                                         }
//                                         h1->SetBinContent(currX, currY, 0);
//                                         esum+=curre;
//                                         stepsctr++;
//                                         if(stepsctr==stepstodo) {
//                                                 stepsctr=0;
//                                                 buf=xdir; //rotate 90 deg
//                                                 xdir= -ydir;
//                                                 ydir=buf;
//
//                                                 if(ydir==0) { //incremant steps at every iteration along x
//                                                         stepstodo++;
//                                                 }
//                                         }
//                                 }
//
//                                 Double_t e=h2->Integral(); //get energy contained in cluster
//
//                                 Eweight=e/Eges; //calculate weight
//
//                                 cgx=h2->GetMean(1);
//                                 xerr=h2->GetMeanError(1);
//                                 //extract center of gravity and error
//                                 cgy=h2->GetMean(2);
//                                 yerr=h2->GetMeanError(2);
//
//                                 if(yerr<0.00001) {yerr=1/TMath::Sqrt(12); }
//                                 if(xerr<0.00001) {xerr=1/TMath::Sqrt(12); }
//
//                                 cgz=i;
//
//                                 er1->Fill(i, xerr/nofEntries);
//                                 er2->Fill(i, yerr/nofEntries);
//
//                                 h3->Fill(cgx, cgy, cgz);
//
//                                 auto cg=std::make_tuple(cgx,cgy,cgz, xerr, yerr, Eweight);
//                                 coglist.push_back(cg);
//
//                         } // end of clustering
//
//                         ClusteredHits.clear();
//
//                         h1->Reset();
//
//                         h2->Reset();
//
//                 }
//                 COGCollection.push_back(coglist);
//                 coglist.clear();
//
//         }
//         plotCOGs();
// }
//
//
// void TROOTAnalysis::FitCOGs( Int_t minevent, Int_t maxevent){
//         //histograms for correlation
//         TCanvas * corr1 = new TCanvas("Correlations");
//         corr1->Divide(3,2,0.01,0.01);
//
//         TH1D * co1 = new TH1D("MxTx correllation", "MxTx correllation",100, -1,1 );
//         co1->GetXaxis()->SetTitle("corellation of X slope and X intercept");
//         TH1D * co2 = new TH1D("MxMy correllation", "MxMy correllation",100, -1,1 );
//         co2->GetXaxis()->SetTitle("corellation of X slope and Y slope");
//         TH1D * co3 = new TH1D("MxTy correllation", "MxTy correllation",100, -1,1 );
//         co3->GetXaxis()->SetTitle("correlation of X slope and Y intercept");
//         TH1D * co4 = new TH1D("TxMy correllation", "TxMy correllation",100, -1,1 );
//         co4->GetXaxis()->SetTitle("correlation of X intercept and Y slope");
//         TH1D * co5 = new TH1D("TxTy correllation", "TxTy correllation",100, -1,1 );
//         co5->GetXaxis()->SetTitle("correlation of X intercapt an Y intercept");
//         TH1D * co6 = new TH1D("MyTy correllation", "MyTy correllation",100, -1,1 );
//         co6->GetXaxis()->SetTitle("correlation of Y slope and Y intercept");
//
//
//
//
//         //transform from copynumber coordinates to geant4 coordinates
//         std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > Transcoglist;
//         std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > > TransfomedCOGs;
//
//
//
//
//         Double_t offsetX=((-1*calsizeXY/2)+tiledimX/2);
//         Double_t offsetY=((-1*calsizeXY/2)+tiledimY/2);
//         Double_t offsetZ=((((AbsoThickness+GapThickness)*(-1)*nofLayers)/2)+(AbsoThickness+GapThickness/2));
//         std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;
//
//         for(Int_t i = 0; i<nofEntries; i++) {                                         //transformation to geant4 coordinate system
//                 for(Int_t j = 0; j<COGCollection[i].size(); j++) {
//
//                         auto tc = std::make_tuple((offsetX + std::get<0>(COGCollection[i][j]) * tiledimX),                // X
//                                                   (offsetY + std::get<1>(COGCollection[i][j]) * tiledimY),                 // Y
//                                                   (offsetZ + std::get<2>(COGCollection[i][j]) * (AbsoThickness+GapThickness)),                               // Z
//                                                   std::get<3>(COGCollection[i][j]),
//                                                   std::get<4>(COGCollection[i][j]),
//                                                   std::get<5>(COGCollection[i][j])  );
//
//                         Transcoglist.push_back(tc);
//
//                 }
//                 TransfomedCOGs.push_back(Transcoglist);
//                 Transcoglist.clear();
//         }
//
//
//
//
//
//         Fcn myfcn;
//         myfcn.SetCOGs(TransfomedCOGs);
//         myfcn.SetMode("photon");
//         //myfcn.PrintCOGs();
//         MnUserParameters upar;
//         double error_minimizer_parameters = 1e-4;
//
//         upar.Add("mx", 0., error_minimizer_parameters);
//         upar.Add("tx", 1., error_minimizer_parameters);
//         upar.Add("my", 0., error_minimizer_parameters);
//         upar.Add("ty", 1., error_minimizer_parameters);
//
//         cout << upar << endl;
//
//         MnMigrad migrad(myfcn, upar, 2);
//
//
//
//         for(int events = 0; events < maxevent-minevent; events++) {
//
//                 upar.SetValue("mx", 0.);
//                 upar.SetValue("tx", 0.);
//                 upar.SetValue("my", 0.);
//                 upar.SetValue("ty", 0.);
//
//                 myfcn.SetCurrentEvent(events);
//
//                 std::cout<<"event: "<<events<<std::endl;
//                 FunctionMinimum min = migrad();
//
//                 std::cout<<min<<std::endl;
//
//
//                 //MnUserCovariance cov = min.UserCovariance();
//
//                 //std::cout<<cov(0 ,0)<<" "<<cov(1,1)<<" "<<cov(2,2)<<" "<<cov(3,3)<<std::endl;
//
//
//                 MnUserParameterState userParameterState = min.UserState();
//
//                 auto tp = std::make_tuple(userParameterState.Value("mx"), userParameterState.Value("tx"),
//                                           userParameterState.Value("my"), userParameterState.Value("ty"));
//                 FitParamsGamma.push_back(tp);
//
//
//
//                 // Double_t covar[4][4];
//                 // Double_t error[4];
//                 //
//                 // error[0]= userParameterState.Error("mx");
//                 // error[1]= userParameterState.Error("tx");
//                 // error[2]= userParameterState.Error("my");
//                 // error[3]= userParameterState.Error("ty");
//                 //
//                 // //Fill correlation Histograms
//                 // if(COGCollection[events].size() > 1) {
//                 //
//                 //
//                 //         for(Int_t row=0; row<4; row++) {
//                 //                 for(Int_t col=0; col<4; col++) {
//                 //                         covar[row][col]=cov(row,col);
//                 //
//                 //                 }
//                 //         }
//                 //
//                 //         co1->Fill(covar[0][1]/(error[0]*error[1]));
//                 //         co2->Fill(covar[0][2]/(error[0]*error[2]));
//                 //         co3->Fill(covar[0][3]/(error[0]*error[3]));
//                 //         co4->Fill(covar[1][2]/(error[1]*error[2]));
//                 //         co5->Fill(covar[1][3]/(error[1]*error[3]));
//                 //         co6->Fill(covar[2][3]/(error[2]*error[3]));
//                 // }
//                 // auto tp=std::make_tuple(userParameterState.Value("a1"),userParameterState.Value("a2"),userParameterState.Value("a3"),
//                 //                         userParameterState.Value("x1"),userParameterState.Value("x2"),userParameterState.Value("x3"));
//                 //
//                 // FitParams.push_back(tp);
//         }
//
//
//
//         // corr1->cd(1);
//         // co1->Draw("");
//         //
//         // corr1->cd(2);
//         // co2->Draw("");
//         //
//         // corr1->cd(3);
//         // co3->Draw("");
//         //
//         // corr1->cd(4);
//         // co4->Draw("");
//         //
//         // corr1->cd(5);
//         // co5->Draw("");
//         //
//         // corr1->cd(6);
//         // co6->Draw("");
//
// }
