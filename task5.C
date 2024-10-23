void task5() {
    using namespace TMath;
    

    TH1D *hist = new TH1D("hist", "histogram", 50, 0, 4);


    TRandom3 *rnd = new TRandom3();


    const int num_points = 10000; 
    double A = 1.0 / exp(1.0);  
    int count_f = 0; 
    int count_fg = 0; 
    

    auto f = [](double x) -> double {
        return exp(-x) * pow(x, 3);
    };
    auto g = [](double x) -> double {
        return (1.0 / 3.0) * pow(x, 2);
    };


    for (int i = 0; i < num_points; i++) {

        double r1 = rnd->Uniform(0, 1);  
        double r2 = rnd->Uniform(0, 1 / exp(1.0));  


        if (r2 <= f(r1)) {
            count_f++;
        }


        double r3 = pow(r1, 1.0 / 3.0);  
        double r4 = rnd->Uniform(0, 3 / exp(1.0));  


        if (r4 <= f(r3) / g(r3)) {
            count_fg++;
        }

     
        hist->Fill(Sqrt(-Log(1 - r1)));  
    }


    double integral_f = ((double)count_f / num_points) * A;
    double integral_fg = ((double)count_fg / num_points) * 3.0/exp(1.0)) ;  

   
    std::cout << "Приближенное значение интеграла: " << integral_f << std::endl;
    std::cout << "Приближенное значение интеграла: " << integral_fg << std::endl;

 
    hist->Draw();
}
