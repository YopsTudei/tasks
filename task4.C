#include <iostream>
#include <cmath>
#include <TRandom3.h>
#include <TMath.h>
#include <TH1D.h>
#include <TCanvas.h>

double calculatePiBuffon(int N, double L, double l) {
    TRandom3 rnd;
    int count_intersections = 0;
    double double_crossing_angle = std::acos(L / l);
    double epsilon = 1e-6;  

    for (int i = 0; i < N; i++) {
        double x = rnd.Uniform(0, L); 
        double phi = rnd.Uniform(-TMath::Pi() / 2.0, TMath::Pi() / 2.0);

        
        if (x <= l * std::cos(phi)) {
            count_intersections++;
        }

        
        if (std::abs(x - L) < epsilon && std::abs(phi) <= double_crossing_angle) {
            count_intersections += 2;
        }
    }

    if (count_intersections == 0) {
        std::cout << "Нет пересечений" << std::endl;
        return 0;
    }

    return (2.0 * N * l) / (count_intersections * L);
}

void task4() {  
    int N = 530;
    double L = 4.0;
    double l = 6.0;

    double pi_estimate = calculatePiBuffon(N, L, l);
    std::cout << "π : " << pi_estimate << std::endl;

    TH1D *hist = new TH1D("hist", " Pi Errors", 50, 3.0, 3.5);
    for (int i = 0; i < 1000; i++) {
        pi_estimate = calculatePiBuffon(N, L, l);
        hist->Fill(pi_estimate);
    }

    TCanvas *canvas = new TCanvas("canvas", "Experiment", 800, 600);
    hist->GetXaxis()->SetTitle("Pi");
    hist->GetYaxis()->SetTitle("Frequency");
    hist->Draw();
}
