#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>

void plot_energy_distribution(TH1F *hist_eph) {
    TCanvas *c1 = new TCanvas("c1", "Energy Distribution", 800, 600);
    hist_eph->GetXaxis()->SetTitle("Photon Energy");
    hist_eph->GetYaxis()->SetTitle("Counts");
    c1->SetLogx();
    hist_eph->Draw();
    hist_eph->SaveAs("energy_distribution.png");
    std::cout << "Energy distribution saved as energy_distribution.png" << std::endl;
    delete c1;
}

void fit_energy_spectrum(TH1F *hist_eph) {
    TCanvas *c2 = new TCanvas("c2", "Energy Spectrum Fit", 800, 600);
    TF1 *fit_func = new TF1("fit_func", "[0] * log(x) + [1]", 0.1, 22000);
    int fit_status = hist_eph->Fit(fit_func, "R");
    fit_func->Draw("same");
    std::cout << "Energy spectrum fitted with function: " << fit_func->GetExpFormula() << std::endl;
    hist_eph->SaveAs("energy_spectrum_fitted.png");
    std::cout << "Energy spectrum with fit saved as energy_spectrum_fitted.png" << std::endl;
    if (fit_status != 0) {
        std::cerr << "Fit did not converge successfully!" << std::endl;
    }
    delete fit_func;
    delete c2;
}

void process() {
    gROOT->SetBatch(kTRUE);
    TFile *file_in = TFile::Open("m3pimc.root");
    if (!file_in || file_in->IsZombie()) { std::cerr << "Failed to open input file!" << std::endl; return; }
    TTree *tree_in = (TTree*)file_in->Get("h10");
    if (!tree_in) { std::cerr << "Tree not found!" << std::endl; return; }
    TFile *file_out = new TFile("newroot.root", "RECREATE");
    int nph; float eph[100], thetaph[100], phiph[100], chi2_3p; UChar_t Isrfilter;
    tree_in->SetBranchAddress("nph", &nph);
    tree_in->SetBranchAddress("eph", eph);
    tree_in->SetBranchAddress("thetaph", thetaph);
    tree_in->SetBranchAddress("phiph", phiph);
    tree_in->SetBranchAddress("chi2_3p", &chi2_3p);
    tree_in->SetBranchAddress("Isrfilter", &Isrfilter);
    tree_in->SetBranchStatus("*", 0);
    tree_in->SetBranchStatus("nph", 1);
    tree_in->SetBranchStatus("eph", 1);
    tree_in->SetBranchStatus("thetaph", 1);
    tree_in->SetBranchStatus("phiph", 1);
    tree_in->SetBranchStatus("Isrfilter", 1);
    tree_in->SetBranchStatus("chi2_3p", 1);
    TTree *filtered_tree = tree_in->CopyTree("Isrfilter == 1 && chi2_3p < 30");
    if (!filtered_tree) { std::cerr << "Failed to copy tree with the specified filter!" << std::endl; file_in->Close(); file_out->Close(); return; }
    filtered_tree->SetName("MyTree");
    filtered_tree->Write();
    std::cout << "Number of entries in MyTree: " << filtered_tree->GetEntries() << std::endl;
    TH1F *hist_eph = new TH1F("hist_eph", "Energy distribution of photons", 200, 0, 22000);
    Long64_t nentries = filtered_tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        filtered_tree->GetEntry(i);
        for (int j = 0; j < nph && j < 100; ++j) {
            hist_eph->Fill(eph[j]);
        }
        if (nph > 100) { std::cerr << "Warning: nph (" << nph << ") exceeds array size!" << std::endl; }
    }
    std::cout << "Number of entries in hist_eph: " << hist_eph->GetEntries() << std::endl;
    plot_energy_distribution(hist_eph);
    fit_energy_spectrum(hist_eph);
    double actual_min = filtered_tree->GetMinimum("eph");
    double actual_max = filtered_tree->GetMaximum("eph");
    std::cout << "Actual minimum photon energy: " << actual_min << std::endl;
    std::cout << "Actual maximum photon energy: " << actual_max << std::endl;
    hist_eph->Write();
    file_out->Write();
    double size_in = file_in->GetSize(); 
    double size_out = file_out->GetSize();
    if (size_out != 0) { 
        std::cout << "File size reduction ratio: " << size_in / size_out << std::endl; 
    }
    else { 
        std::cout << "Output file size is zero!" << std::endl; 
    }
    delete hist_eph;
    file_out->Close();
    file_in->Close();
}

int main() {
    process();
    return 0;
}
