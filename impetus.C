#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TMultiGraph.h>
#include <iostream>
#include <cmath>
#include <vector>

const double r = 0.3;  // радиус цилиндра, в нём действует магнитное поле
const double H = 1.0;  // расстояние от центра цилиндра до середины октанта (октант - сторона правильного восьмиугольника)
const double p_total = 793.0; // импульс мюона, из закона соххранения энергии
const double detector_size = 0.5; // размер каждого окнтанта я считаю 1х1 метр

void go() {
    TFile *f = TFile::Open("tr_ph_run132028.root");
    TTree *tree = (TTree*)f->Get("tr_ph");
    int nsim, nmu;
    float simphi[10], simtheta[10], simvtz[10], tcharge[10], simmom[10];
    int much[10];
    tree->SetBranchAddress("nsim", &nsim);
    tree->SetBranchAddress("simphi", simphi);
    tree->SetBranchAddress("simtheta", simtheta);
    tree->SetBranchAddress("simvtz", simvtz);
    tree->SetBranchAddress("simmom", simmom);
    tree->SetBranchAddress("tcharge", tcharge);
    tree->SetBranchAddress("nmu", &nmu);
    tree->SetBranchAddress("much", much);

    TCanvas *c1 = new TCanvas("c1", "Muon Hit on Detector", 800, 600);
    TMultiGraph *mg = new TMultiGraph();
    TGraph *gHits = new TGraph();
    TGraph *gMisses = new TGraph();

    int total = 0;
    int hitCount = 0;

    const int maxMuons = 1000;
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries && total < maxMuons; i++) {
        tree->GetEntry(i);
        for (int j = 0; j < nsim && total < maxMuons; j++) {
            double theta = simtheta[j];
            double phi = simphi[j];
            double charge = tcharge[j];
// тут я пытался отсеять те mu+ с определенными углами
            if (charge > 0 || simmom[j] < 10) continue;
            if (theta < M_PI/6 || theta > 5*M_PI/6) continue;
            if (phi < -M_PI/4 || phi > M_PI/4) continue;
// тут я вычисляю куда попадает мюон, на примерах конкретных мюонов считалось правильно
            double p_perp = p_total * sin(theta);
            double R = p_perp / 3900.0; //радиус кривизны в магнитном поле, вычисляется как R=p_perp(MeV)/(0,3* 13Gs)=...(sm)=.../100(m)
            double alpha = 2 * asin(r / (2 * R)); //угол, на который происходит поворот в области действия магнитного . поле вдоль оси z, значит отклонение происходит только в плоскости x,y
            double yr = r * sin(phi-alpha/2); //новая Y координата мюона
            double xr = r*cos(phi-alpha/2);//новая X координата мюона 
            double h = H - yr; //вертикальное расстояние от мюона до октанта
            double Xd = xr + tan(alpha) * h; //координата Х вдоль октанта - на примере: вдоль оси Х расположены детекторы 1, 2, 35, 36
            double L_perpendicular = h * cos(alpha); //расстояние, которое пролетает мюон, выйдя из цилиндра с магнитным полем, до октанта
            double Yd = simvtz[j] + L_perpendicular * 1.0 / tan(theta); //по сути исходная ось Z

            std::cout << "Muon " << j << ": simmom = " << simmom[j] << ", phi = " << phi << ", theta = " << theta << std::endl;
          //проверяю, относится ли частица к задетектированным, не уверен, что правильно, тут и нужна помощь
            bool hit = false;
            for (int k = 0; k < nmu; k++) {
                if (much[k] == j) {
                    hit = true;
                    break;
                }
            }

            if (hit) {
                gHits->SetPoint(hitCount, Xd, Yd);
                hitCount++;
            } else {
                gMisses->SetPoint(total - hitCount, Xd, Yd);
            }
            total++;
        }
    }

    gHits->SetMarkerStyle(2); // крестик
    gHits->SetMarkerColor(kRed);
    gHits->SetTitle("Hits");

    gMisses->SetMarkerStyle(20); // точка
    gMisses->SetMarkerColor(kBlue);
    gMisses->SetTitle("Misses");

    mg->Add(gHits);
    mg->Add(gMisses);
    mg->SetTitle("Muon Hits on Detector;X (m);Y (m)");
    mg->Draw("AP");

    c1->Update();
}
