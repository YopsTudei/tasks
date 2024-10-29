#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TStyle.h>

void analyze_data() {
    // Загружаем экспериментальные данные
    std::ifstream fdata("fdata.dat");
    std::ifstream fMC("fMC.dat");

    // Создаем двумерные гистограммы для данных
    TH2F *h_exp = new TH2F("h_exp", "Experimental Data", 100, -5, 5, 100, -5, 5);
    TH2F *h_mc = new TH2F("h_mc", "MC Data", 100, -5, 5, 100, -5, 5);

    // Считываем данные и заполняем гистограммы
    double x, y;
    while (fdata >> x >> y) {
        h_exp->Fill(x, y);
    }

    while (fMC >> x >> y) {
        h_mc->Fill(x, y);
    }

    // Создаем холст
    TCanvas *c1 = new TCanvas("c1", "Data Analysis", 1200, 800);
    c1->Divide(2, 2); // Разделяем холст на 4 части

    // Рисуем экспериментальную гистограмму
    c1->cd(1);
    h_exp->Draw("COLZ"); // Двумерная гистограмма для эксперимента
    h_exp->SetTitle("Experimental Data");
    h_exp->GetXaxis()->SetTitle("X");
    h_exp->GetYaxis()->SetTitle("Y");

    // Рисуем гистограмму моделирования
    c1->cd(2);
    h_mc->Draw("COLZ"); // Двумерная гистограмма для моделирования
    h_mc->SetTitle("MC Data");
    h_mc->GetXaxis()->SetTitle("X");
    h_mc->GetYaxis()->SetTitle("Y");

    // Проекции на ось X
    TH1D *h_exp_projX = h_exp->ProjectionX();
    TH1D *h_mc_projX = h_mc->ProjectionX();

    // Нормируем моделирование так, чтобы максимальный бин совпадал
    double scale_factor = h_exp_projX->GetMaximum() / h_mc_projX->GetMaximum();
    h_mc_projX->Scale(scale_factor);

    // Рисуем проекции на одном полотне
    c1->cd(3);
    h_exp_projX->SetLineColor(kRed);
    h_mc_projX->SetLineColor(kBlue);
    h_exp_projX->Draw();
    h_mc_projX->Draw("SAME");
    h_exp_projX->SetTitle("Projections on X-axis");
    h_exp_projX->GetXaxis()->SetTitle("X");
    h_exp_projX->GetYaxis()->SetTitle("Counts");

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(h_exp_projX, "Experimental", "l");
    legend->AddEntry(h_mc_projX, "MC (normalized)", "l");
    legend->Draw();

    // Разница гистограмм
    TH1D *h_diff = (TH1D*)h_exp_projX->Clone("h_diff");
    h_diff->Add(h_mc_projX, -1);

    // Рисуем разницу
    c1->cd(4);
    h_diff->Draw();
    h_diff->SetTitle("Difference between Experimental and MC");
    h_diff->GetXaxis()->SetTitle("X");
    h_diff->GetYaxis()->SetTitle("Counts");

    // Сохраняем холст
    c1->SaveAs("data_analysis.pdf");

    // Подсчет фона
    double background = h_diff->Integral();
    std::cout << "Количество фона: " << background << std::endl;
}
