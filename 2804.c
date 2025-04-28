#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TF1.h>
#include <TH2F.h>
#include <iostream>
#include <vector>
#include <cmath>

const double PI = 3.14159265358979;
const double H = 0.985;
const int N = 1000; 
const double r = 0.3;

double deltaPhi(double phi1, double phi2) {
    double dphi = phi1 - phi2;
    while (dphi > PI) dphi -= 2 * PI; 
    while (dphi < -PI) dphi += 2 * PI;
    return dphi;
}

void work() {
    TFile *f = TFile::Open("tr_ph_run132028.root");
    TTree *tree = (TTree*)f->Get("tr_ph");

    int nt, nsim;
    int it[2];
    float simphi[10], simtheta[10], tptot[10];
    Int_t tcharge[10];
    float simvtz[10];
    int nmu;
    int much[10];

    tree->SetBranchAddress("nt", &nt);
    tree->SetBranchAddress("nsim", &nsim);
    tree->SetBranchAddress("it", it);
    tree->SetBranchAddress("simphi", simphi);
    tree->SetBranchAddress("simtheta", simtheta);
    tree->SetBranchAddress("simvtz", simvtz);
    tree->SetBranchAddress("tptot", tptot);
    tree->SetBranchAddress("tcharge", tcharge);
    tree->SetBranchAddress("nmu", &nmu);
    tree->SetBranchAddress("much", much);

    Long64_t nentries = tree->GetEntries();
    int selected = 0;

    std::vector<double> z_points;
    std::vector<double> y_points;
    std::vector<double> z_hits;
    std::vector<double> y_hits;
    std::vector<double> z_misses;
    std::vector<double> y_misses;
    std::vector<int> detectors_hit;

    for (Long64_t i = 0; i < nentries && selected < N; i++) {
        tree->GetEntry(i);
        if (nt != 2) continue;

        int i1 = it[0];
        int i2 = it[1];

        float phi1 = simphi[i1];
        float phi2 = simphi[i2];
        float theta1 = simtheta[i1];
        float theta2 = simtheta[i2];
        float p1 = tptot[i1];
        float p2 = tptot[i2];
        float q1 = tcharge[i1];
        float q2 = tcharge[i2];

        int ipos = -1;
        if (q1 > 0) ipos = i1;
        else if (q2 > 0) ipos = i2;
        else continue;


        float phi_pos = simphi[ipos];
        float theta_pos = simtheta[ipos];
        double ptot_pos = tptot[ipos];

        if (ptot_pos < 790 || ptot_pos > 800) continue;

        bool in_phi_range =
            (phi_pos > 0 && phi_pos < 0.1) || 
            (phi_pos > 2*PI  && phi_pos < 2 * PI);
        if (!in_phi_range) continue;

        if (!(theta_pos > PI / 4 && theta_pos < 2*PI-PI/4)) continue;

        double R = ptot_pos * sin(theta_pos) / 390;
        double alpha = (-1)* 2 * asin(r / (2 * R)); 
        double xr = r * cos(phi_pos - alpha / 2);
        double yr = r * sin(phi_pos - alpha / 2);

        double y = yr + tan(phi_pos - alpha) * (H - xr);
        double z = (1 / tan(theta_pos)) * (R * (-1)*alpha + (H - xr) / cos(phi_pos - alpha)) + simvtz[ipos]/100; 

       
        if (fabs(y) <= 1 && fabs(z) <= 1) { 
            detectors_hit.clear();
            bool hit = false;
                int detectorID = much[1]; //число сработавших счётчиков ровно один
                if (detectorID == 1 || detectorID == 2 || detectorID == 35 || detectorID == 36) { //проверяю сработали ли нужные счётчики
                    detectors_hit.push_back(detectorID);
                    hit = true;
                }
            

            if (hit) {
                z_hits.push_back(z);
                y_hits.push_back(y);
            } else {
                z_misses.push_back(z);
                y_misses.push_back(y);
            }

            z_points.push_back(z);
            y_points.push_back(y);

           
            std::cout << "Event " << selected + 1 
            << " | phi = " << simphi[ipos]
            << " | simvtz = " << simvtz[ipos]/100
            << " | theta = " << simtheta[ipos]
            << " | impulse = " << tptot[ipos]
            << " | alpha = " << alpha
            << " | y = " << y
            << " | z = " << z
            << " | detected = ";

            if (!detectors_hit.empty()) {
                for (size_t j = 0; j < detectors_hit.size(); ++j) {
                    std::cout << detectors_hit[j];
                    if (j < detectors_hit.size() - 1) std::cout << ", ";
                }
            } else {
                std::cout << "none";
            }

            std::cout << "\n";
            selected++;
        }
    }

    
    TCanvas *canvas = new TCanvas("canvas", "Detector Visualization", 800, 600);
    canvas->SetGrid();

    TH2F *frame = new TH2F("frame", "Track Map;z [m];y [m]", 10, -1, 1, 10, -1, 1); 
    frame->SetStats(0);
    frame->Draw();

 
    TGraph *hitGraph = new TGraph(z_hits.size(), &z_hits[0], &y_hits[0]);
    hitGraph->SetMarkerStyle(3);
    hitGraph->SetMarkerColor(kBlack);
    hitGraph->SetMarkerSize(1);
    hitGraph->Draw("P same");

  
    TGraph *missGraph = new TGraph(z_misses.size(), &z_misses[0], &y_misses[0]);
    missGraph->SetMarkerStyle(20);
    missGraph->SetMarkerColor(kRed);
    missGraph->SetMarkerSize(1);
    missGraph->Draw("P same");

    TLine *l1 = new TLine(-0.75, -0.4, -0.75,  0.4); 
    TLine *l2 = new TLine(-0.75,  0.4,  0.75,  0.4);
    TLine *l3 = new TLine( 0.75,  0.4,  0.75, -0.4); 
    TLine *l4 = new TLine( 0.75, -0.4, -0.75, -0.4); 
    l1->SetLineStyle(2); l2->SetLineStyle(2);
    l3->SetLineStyle(2); l4->SetLineStyle(2);
    l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();
    canvas->Update();
}
