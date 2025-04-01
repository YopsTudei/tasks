#include <TCanvas.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TMultiGraph.h>
#include <iostream>
#include <vector>
#include <cmath>

const int numMuons = 500;       // Количество мюонов
const double r = 0.3;  // Радиус цилиндра (м)
const double R = 1.44; // Радиус кривизны (м)
const double H = 1.0; // Расстояние до детекторов (м)
const int numDetectors = 8;
const double p_total = 793;
const double x_min = -2 * tan(M_PI / 8);
const double x_max = 2 * tan(M_PI / 8);

void go() {
    TRandom3 randGen;
    TCanvas *c1 = new TCanvas("c1", "Muon Hits on Detectors", 1200, 800);
    c1->Divide(4, 2);

    std::vector<TGraph*> graphs;
    for (int i = 0; i < numDetectors; i++) {
        graphs.push_back(new TGraph());
    }
    
    for (int i = 0; i < numMuons; i++) {
        double theta = randGen.Uniform(0, M_PI / 2);  // Полусфера вверх
        double phi = randGen.Uniform(0, 2 * M_PI); // Вокруг оси Z
        
        // Генерация начального импульса
        double p_perp = p_total * sin(theta);
        double pz = p_total * cos(theta);
        
        
        // Перевод в градусы
        double phi_deg = phi * 180.0 / M_PI;
        int n = static_cast<int>(phi_deg) / 45; // Определение приблизительного детектора
        
        if (phi_deg - 11.7 > 0) {
            n = (n + 1) % numDetectors;
            phi_deg = 11.7 - phi_deg;
        }
        
        // Перевод обратно в радианы
        double phi_final = phi_deg * M_PI / 180.0;
        
        // Вычисление координат
        double Xd = r * sin(phi_final +  r/(2*R)) + 
                    (H - r * cos(phi_final + r /(2*R))) * tan(phi_final + r / R);
        double L_perpendicular = sqrt(H * H + Xd * Xd);
        double Yd = L_perpendicular * 1.0 / tan(theta);
        
        graphs[n]->SetPoint(graphs[n]->GetN(), Xd, Yd);
    }
    
    for (int i = 0; i < numDetectors; i++) {
        c1->cd(i + 1);
        graphs[i]->SetTitle(Form("Detector %d;X (m);Y (m)", i + 1));
        graphs[i]->SetMarkerStyle(7);
        graphs[i]->GetXaxis()->SetLimits(x_min, x_max);
        graphs[i]->Draw("AP");
    }

    c1->Update();
}
