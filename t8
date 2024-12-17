#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>

Float_t invMass(Float_t e1, Float_t e2, Float_t theta1, Float_t phi1, Float_t theta2, Float_t phi2) {
    Float_t p1 = e1 * sin(theta1);
    Float_t p2 = e2 * sin(theta2);
    
    Float_t px1 = p1 * cos(phi1);
    Float_t py1 = p1 * sin(phi1);
    Float_t pz1 = e1 * cos(theta1);
    
    Float_t px2 = p2 * cos(phi2);
    Float_t py2 = p2 * sin(phi2);
    Float_t pz2 = e2 * cos(theta2);
    
    Float_t massSquared = (e1 + e2) * (e1 + e2) - (px1 + px2) * (px1 + px2) - (py1 + py2) * (py1 + py2) - (pz1 + pz2) * (pz1 + pz2);
    return sqrt(massSquared);
}

Float_t angleBetweenPhotons(Float_t theta1, Float_t phi1, Float_t theta2, Float_t phi2) {
    Float_t cosTheta = cos(theta1) * cos(theta2) + sin(theta1) * sin(theta2) * cos(phi1 - phi2);
    return acos(cosTheta);
}

void t8() {
    TFile *inputFile = TFile::Open("newroot.root");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Ошибка при открытии файла newroot.root" << std::endl;
        return;
    }
    TTree *inputTree = (TTree *)inputFile->Get("MyTree");
    if (!inputTree) {
        std::cerr << "Дерево MyTree не найдено" << std::endl;
        inputFile->Close();
        return;
    }

    Int_t nph;
    Float_t eph[1000], thetaph[1000], phiph[1000];
    inputTree->SetBranchAddress("nph", &nph);
    inputTree->SetBranchAddress("eph", eph);
    inputTree->SetBranchAddress("thetaph", thetaph);
    inputTree->SetBranchAddress("phiph", phiph);

    
    TFile *outputFile = new TFile("filtered_results.root", "RECREATE");
    TTree *outputTree = new TTree("MyTree", "Filtered Photon Data");

    
    Int_t nPi0Candidates;
    Float_t pi0Mass[10], pi0Angle[10];  
    outputTree->Branch("nPi0Candidates", &nPi0Candidates, "nPi0Candidates/I");
    outputTree->Branch("pi0Mass", pi0Mass, "pi0Mass[nPi0Candidates]/F");
    outputTree->Branch("pi0Angle", pi0Angle, "pi0Angle[nPi0Candidates]/F");

    Long64_t nEntries = inputTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        inputTree->GetEntry(i);

        std::vector<std::pair<int, int>> pi0Candidates;

       
        for (Int_t j1 = 0; j1 < nph; j1++) {
            for (Int_t j2 = j1 + 1; j2 < nph; j2++) {
                Float_t mass = invMass(eph[j1], eph[j2], thetaph[j1], phiph[j1], thetaph[j2], phiph[j2]);

                if (mass >= 0.1 && mass <= 0.2) {
                    pi0Candidates.push_back(std::make_pair(j1, j2));
                }
            }
        }

        
        if (pi0Candidates.size() == 2) {
            nPi0Candidates = pi0Candidates.size();
            for (int i = 0; i < nPi0Candidates; i++) {
                int j1 = pi0Candidates[i].first;
                int j2 = pi0Candidates[i].second;

                pi0Mass[i] = invMass(eph[j1], eph[j2], thetaph[j1], phiph[j1], thetaph[j2], phiph[j2]);
                pi0Angle[i] = angleBetweenPhotons(thetaph[j1], phiph[j1], thetaph[j2], phiph[j2]);
            }
            outputTree->Fill();
        }
    }

    
    outputFile->Write();
    outputFile->Close();
    inputFile->Close();

    std::cout << "Анализ завершен. Результаты сохранены в filtered_results.root" << std::endl;
}
