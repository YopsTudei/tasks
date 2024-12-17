#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>

void process() {
    TString inputDataFile = "m3pimc.root"; 
    TString outputDataFile = "newroot.root"; 

    TFile *inputData = TFile::Open(inputDataFile); 
    TTree *dataTree = (TTree*)inputData->Get("h10"); 
    if (!dataTree) {
        std::cerr << "Tree not found!" << std::endl;
        return;
    }

    TFile *outputData = new TFile(outputDataFile, "RECREATE"); 
    TTree *selectedTree = dataTree->CopyTree("Isrfilter==1 && chi2_3p < 30"); 

    selectedTree->SetBranchStatus("*", 0); 
    selectedTree->SetBranchStatus("nph", 1);
    selectedTree->SetBranchStatus("eph", 1);
    selectedTree->SetBranchStatus("thetaph", 1);
    selectedTree->SetBranchStatus("phiph", 1);

    TTree *finalDataTree = selectedTree->CloneTree();
    finalDataTree->SetName("MyTree");
    finalDataTree->Write();

    TH1F *energyHistogram = new TH1F("energyHistogram", "Photon-energy; MeV; count", 100, 0, 9);
    TF1 *fitFunction = new TF1("fitFunction", "[0]/x^[1]", 0, 9);
    fitFunction->SetParameter(1, 3);
    TCanvas *canvas = new TCanvas();
    finalDataTree->Draw("eph>>energyHistogram", "", "", finalDataTree->GetEntries(), 0);
    canvas->SetLogy();
    energyHistogram->Fit(fitFunction);
    energyHistogram->Write();
    energyHistogram->SaveAs("photon_energy_distribution.png");

    Float_t maxPhotonEnergy = finalDataTree->GetMaximum("eph"); 
    Float_t minPhotonEnergy = finalDataTree->GetMinimum("eph"); 

    std::cout << "Maximum photon energy: " << maxPhotonEnergy << " MeV" << std::endl;
    std::cout << "Minimum photon energy: " << minPhotonEnergy << " MeV" << std::endl;

    double size_in = inputData->GetSize();
    double size_out = outputData->GetSize();
    if (size_out != 0) {
        std::cout << "File size reduction ratio: " << size_in / size_out << std::endl;
    } else {
        std::cout << "Output file size is zero!" << std::endl;
    }

    
    outputData->Close();
    inputData->Close();
}

