#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <cmath>

void task9() {
    TFile *inputFile1 = TFile::Open("filtered_results.root");
    if (!inputFile1 || inputFile1->IsZombie()) {
        std::cerr << "Error opening file filtered_results.root" << std::endl;
        return;
    }
    TFile *inputFile2 = TFile::Open("newroot.root");
    if (!inputFile2 || inputFile2->IsZombie()) {
        std::cerr << "Error opening file newroot.root" << std::endl;
        return;
    }
    TTree *outputTree = (TTree *)inputFile2->Get("MyTree;1");
    if (!outputTree) {
        std::cerr << "Tree MyTree not found in newroot.root" << std::endl;
        inputFile1->Close();
        inputFile2->Close();
        return;
    }

    Int_t nph;
    Float_t eph[1000], thetaph[1000], phiph[1000];
    outputTree->SetBranchAddress("nph", &nph);
    outputTree->SetBranchAddress("eph", eph);
    outputTree->SetBranchAddress("thetaph", thetaph);
    outputTree->SetBranchAddress("phiph", phiph);

    TH1F *hPolarAngle = new TH1F("hPolarAngle", "Number of #pi^{0} candidates vs polar angle", 36, 0, 180);
    TH1F *hAzimuthalAngle = new TH1F("hAzimuthalAngle", "Number of #pi^{0} candidates vs azimuthal angle", 36, -180, 180);
    TH1F *hNumCandidates = new TH1F("hNumCandidates", "Distribution of pion count", 50, 0, 100);

    Long64_t nEntries = outputTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        outputTree->GetEntry(i);
        Int_t numCandidates = 0;
        for (Int_t j1 = 0; j1 < nph; j1++) {
            for (Int_t j2 = j1 + 1; j2 < nph; j2++) {
                Float_t e1 = eph[j1], e2 = eph[j2];
                Float_t p1[3] = {e1 * sin(thetaph[j1]) * cos(phiph[j1]), e1 * sin(thetaph[j1]) * sin(phiph[j1]), e1 * cos(thetaph[j1])};
                Float_t p2[3] = {e2 * sin(thetaph[j2]) * cos(phiph[j2]), e2 * sin(thetaph[j2]) * sin(phiph[j2]), e2 * cos(thetaph[j2])};
                Float_t mass2 = (e1 + e2) * (e1 + e2) - (p1[0] + p2[0]) * (p1[0] + p2[0]) - (p1[1] + p2[1]) * (p1[1] + p2[1]) - (p1[2] + p2[2]) * (p1[2] + p2[2]);
                Float_t mass = sqrt(mass2);
                if (mass > 0.1 && mass < 0.2) {
                    numCandidates++;
                    Float_t deltaPhi = phiph[j1] - phiph[j2];
                    if (deltaPhi > TMath::Pi()) deltaPhi -= 2 * TMath::Pi();
                    if (deltaPhi < -TMath::Pi()) deltaPhi += 2 * TMath::Pi();
                    Float_t deltaTheta = thetaph[j1] - thetaph[j2];
                    hPolarAngle->Fill(fabs(deltaTheta) * 180 / TMath::Pi());
                    hAzimuthalAngle->Fill(deltaPhi * 180 / TMath::Pi());
                }
            }
        }

        hNumCandidates->Fill(numCandidates);
    }

    Float_t varCandidates = hNumCandidates->GetStdDev();
    std::cout << "Variance for number of pions: " << varCandidates << std::endl;

    TCanvas *canvas1 = new TCanvas("canvas1", "Pi0 hist", 800, 600);
    hNumCandidates->Draw();
    canvas1->SaveAs("new_directory/num_candidates_distribution.png");

    Float_t varPolar = hPolarAngle->GetStdDev();
    Float_t varAzimuth = hAzimuthalAngle->GetStdDev();
    std::cout << "Variance for polar angle: " << varPolar << std::endl;
    std::cout << "Variance for azimuthal angle: " << varAzimuth << std::endl;

    TGraphErrors *graphPolar = new TGraphErrors();
    TGraphErrors *graphAzimuth = new TGraphErrors();
    for (Int_t bin = 1; bin <= hPolarAngle->GetNbinsX(); bin++) {
        graphPolar->SetPoint(bin - 1, hPolarAngle->GetBinCenter(bin), hPolarAngle->GetBinContent(bin));
        graphPolar->SetPointError(bin - 1, 0, hPolarAngle->GetBinError(bin));
    }
    for (Int_t bin = 1; bin <= hAzimuthalAngle->GetNbinsX(); bin++) {
        graphAzimuth->SetPoint(bin - 1, hAzimuthalAngle->GetBinCenter(bin), hAzimuthalAngle->GetBinContent(bin));
        graphAzimuth->SetPointError(bin - 1, 0, hAzimuthalAngle->GetBinError(bin));
    }

    TCanvas *canvas = new TCanvas("canvas", "Graph with Errors", 800, 600);
    graphPolar->SetMarkerStyle(21);  
    graphPolar->SetMarkerColor(kBlue);
    graphPolar->SetLineColor(kBlue);
    graphPolar->SetTitle("Angle vs Candidate Count");
    graphPolar->GetXaxis()->SetTitle("Angle (degrees)");
    graphPolar->GetYaxis()->SetTitle("Candidate Count");
    graphPolar->Draw("AP");
    graphAzimuth->SetMarkerStyle(21);  
    graphAzimuth->SetMarkerColor(kRed);
    graphAzimuth->SetLineColor(kRed);
    graphAzimuth->Draw("P");
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(graphPolar, "Polar Angle", "p");
    legend->AddEntry(graphAzimuth, "Azimuthal Angle", "p");
    legend->Draw();
    canvas->SaveAs("new_directory/results_task9.png");

    inputFile1->Close();
    inputFile2->Close();
}

