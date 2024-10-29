void task5() {
    using namespace TMath;
    const int num_points = 10000;
    const int N = 100; 
    const double true_integral = 0.11393;  

    TRandom3 *rnd = new TRandom3();


    TH1D *hist_error_neyman = new TH1D("error_neyman", "Error in Neyman Method", 50, -0.02, 0.02);
    TH1D *hist_error_importance = new TH1D("error_importance", "Error in Importance Sampling", 50, -0.02, 0.02);
    TH1D *hist_error_mean = new TH1D("error_mean", "Error in Mean Method", 50, -0.02, 0.02);

   
    auto f = [](double x) -> double {
        return exp(-x) * pow(x, 3);
    };
    auto g = [](double x) -> double {
        return (log(8) / 56.0) * pow(8, x + 1);
    };

    for (int j = 0; j < N; j++) {
        int count_f = 0;
        int count_fg = 0;
        double sum_f = 0;

        
        for (int i = 0; i < num_points; i++) {
            double r1 = rnd->Uniform(0, 1);
            double r2 = rnd->Uniform(0, 1 / exp(1.0));
            sum_f += f(r1);

            if (r2 <= f(r1)) {
                count_f++;
            }

            double r3 = log(7 * r1 + 1) / log(8);
            double r4 = rnd->Uniform(0, 0.16);

            if (r4 <= f(r3) / g(r3)) {
                count_fg++;
            }
        }

       
        double integral_f = ((double)count_f / num_points) * (1.0 / exp(1.0));
        double integral_fg = ((double)count_fg / num_points) * 0.16;
        double mean_f = sum_f / num_points;

        
        double error_neyman = true_integral - integral_f;
        double error_importance = true_integral - integral_fg;
        double error_mean = true_integral - mean_f;

        hist_error_neyman->Fill(error_neyman);
        hist_error_importance->Fill(error_importance);
        hist_error_mean->Fill(error_mean);
    }

    
    TCanvas *c1 = new TCanvas("c1", "Error Comparison", 900, 600);
    c1->Divide(3, 1);
    c1->cd(1);
    hist_error_neyman->Draw();
    c1->cd(2);
    hist_error_importance->Draw();
    c1->cd(3);
    hist_error_mean->Draw();

    std::cout << "..." << std::endl;
}
