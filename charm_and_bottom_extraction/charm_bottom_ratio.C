

void charm_bottom_ratio(){


    TFile* fin = new TFile("plotting_bbbar_from_DCA_Fits_with_systematics_corrected.root","READ");
    TGraphErrors* graph_stat = (TGraphErrors*)fin->Get("graph_stat");
    graph_stat->SetName("graph_stat");
    TGraphErrors* graph_sys = (TGraphErrors*)fin->Get("graph_sys");
    graph_sys->SetName("graph_sys");

    TF1* func1 = new TF1("func1","pol0", 2.0, 3.0);
    graph_stat->Fit(func1, "R+");

    cout << "Average Charm yield in the 2.0-3.0 GeV/c Region" << endl;
    cout << func1->GetParameter(0) << "\t" << func1->GetParError(0) << endl;

    std::vector<double> fit_values;

    int nPoints = graph_stat->GetN();

    for (int iter = 0; iter < 2e5; ++iter) {

        std::vector<double> x_vals(nPoints);
        std::vector<double> y_vals(nPoints);
        std::vector<double> y_errs(nPoints);

        double r = gRandom->Gaus(0,1);
        int counter = 0;

        for (int i = 0; i < nPoints; ++i) {

            double x_stat, y_stat, y_sys;
            graph_stat->GetPoint(i, x_stat, y_stat);
            graph_sys->GetPoint(i, x_stat, y_sys);

            if(x_stat < 2.0 || x_stat > 3.0) continue;

            double stat_err = graph_stat->GetErrorY(i);
            double sys_err  = graph_sys->GetErrorY(i);

            double y_shifted = y_stat + r*sys_err;

            x_vals[counter] = x_stat;
            y_vals[counter] = y_shifted;
            y_errs[counter] = stat_err;

            counter++;

        }

    
        TGraphErrors* g_varied = new TGraphErrors(counter);

        for (int i = 0; i < counter; ++i) {

            g_varied->SetPoint(i, x_vals[i], y_vals[i]);
            g_varied->SetPointError(i, 0, y_errs[i]);

        }


        TF1* f_const = new TF1("f_const", "pol0", 2.0, 3.0);
        g_varied->Fit(f_const, "RQ");

        double fit_val = f_const->GetParameter(0);
        fit_values.push_back(fit_val);

        delete f_const;
        delete g_varied;
    }

    double mean = TMath::Mean(fit_values.begin(), fit_values.end());
    double stddev = TMath::StdDev(fit_values.begin(), fit_values.end());

    std::cout << "Mean of fitted pol0 values " << mean << std::endl;
    std::cout << "Systematic uncertainty from coherent sys shift: " << stddev << std::endl;

    TH1D* h_fit_vals = new TH1D("h_fit_vals", "Distribution of Fitted Pol0 Values", 100, mean - 5*stddev, mean + 5*stddev);
    for (double val : fit_values) { h_fit_vals->Fill(val); }
    h_fit_vals->Draw();
}


    
