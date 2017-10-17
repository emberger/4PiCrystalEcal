// Main function for analyzing the TTree

#include <iostream>
#include "TApplication.h"
#include "TStyle.h"
#include "TChain.h"
#include "TROOTAnalysis.hh"
#include <stdlib.h>


int main(int argc, char * argv[]) {

        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

        TApplication* app = new TApplication("app", 0, 0, 0);
//
//    gStyle->SetCanvasBorderMode(0);
//    gStyle->SetPadBorderMode(0);
//
//    gStyle->SetPalette(1);
//    gStyle->SetOptStat(1);
//    gStyle->SetOptFit(11);
//    //gStyle->SetOptTitle(0);
//    gStyle->SetStatBorderSize(1);
//    gStyle->SetStatColor(10);
//    gStyle->SetCanvasColor(10);
//    gStyle->SetPadLeftMargin(0.16);
//    gStyle->SetPadBottomMargin(0.16);
//    gStyle->SetPadTickX(1);
//    gStyle->SetPadTickY(1);
//    gStyle->SetOptTitle(0);
//    gStyle->SetTitleSize(0.048,"xy");
//    gStyle->SetLabelSize(0.04,"xy");
//    gStyle->SetTitleOffset(1.3,"x");
//    gStyle->SetTitleOffset(1.3,"y");

        TChain * ch1 = new TChain("eventTree");
        ch1->Add("AbsoFirst.root");
        ch1->Draw("");
        TROOTAnalysis A(ch1);

        Int_t event=0;
        Int_t progress=0;
        //A.SampleFromHIst();
        for(Int_t i=0; i<A.GetNofEntries(); i++) {

                if(A.CheckEvent(i)) {

                        A.CalcCOGPion(event);

                        A.FitCOGsPion(event);
                        A.FitCOGsPion2(event);

                        std::pair<TVector3, TVector3> ca= A.FindClosestApproach(event);

                        A.GetInvariantMass(event);



                        A.PionLocator(event, ca);

                        event++;
                }


                progress = (i /A.GetNofEntries() ) * 100;
                if (progress % 5 == 0)
                {
                        std::cout << "\r" << std::string(progress/5, '|') << progress << "%";
                        std::cout.flush();
                }
        }



        A.DrawHists();


        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Computing took "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  <<" seconds"<<std::endl;
        app->Run();
        return 0;
}
