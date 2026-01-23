
const int mass_regions = 7;
double mass_lo[mass_regions] = { 1.75, 2.40, 1.60, 2.40, 2.60, 3.10, 2.80};
double mass_hi[mass_regions] = { 2.40, 2.95, 2.80, 3.45, 4.00, 4.25, 5.00};

double x[mass_regions] = {-999};
double y[mass_regions] = {-999};
double err_x[mass_regions] = {-999};
double err_y[mass_regions] = {-999};

double x_dca[mass_regions] = {-999};
double y_dca[mass_regions] = {-999};
double err_x_dca[mass_regions] = {0};
double err_y_dca[mass_regions] = {-999};


void normalizing_two_TGraphs_bbbar(const int stat){


    TFile* file = new TFile("/Users/vassu/Desktop/new_files/data_files/producing_templates_without_20percent_reduction.root","READ");
    TH1D* h1d_bbbar = (TH1D*)file->Get("h1d_ee_FG_bbbar");
    h1d_bbbar->SetName("h1d_bbbar");
    double x_graph, y_graph;

    TH1D* h1d_bbbar_plot = (TH1D*)file->Get("h1d_ee_FG_bbbar");
    h1d_bbbar_plot->SetName("h1d_bbbar_plot");

    for(int ibin = 0; ibin < mass_regions; ibin++){

        int bin_lo = h1d_bbbar->FindBin(mass_lo[ibin] + 0.01);
        int bin_hi = h1d_bbbar->FindBin(mass_hi[ibin] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_bbbar->GetBinContent(i)*h1d_bbbar->GetBinWidth(i);
            bin_error += pow(h1d_bbbar->GetBinError(i)*h1d_bbbar->GetBinWidth(i), 2);

        }

        x[ibin] = (mass_lo[ibin] + mass_hi[ibin])/2;
        y[ibin] = (bin_content)/(mass_hi[ibin] - mass_lo[ibin]);
        err_x[ibin] = (mass_hi[ibin] - mass_lo[ibin]);
        err_y[ibin] = (sqrt(bin_error))/(mass_hi[ibin] - mass_lo[ibin]);


    }

    TGraphErrors* graph = new TGraphErrors(mass_regions, x, y, err_x, err_y);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kRed);
    graph->SetName("graph");

    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
    h1d_bbbar->Draw();
    graph->Draw("same");

    double sum_pythia = 0;
    double sum_pythia_weights = 0;

    for(int i = 0; i < (int)graph->GetN(); i++){

        graph->GetPoint(i, x_graph, y_graph);
        double weight = 1.0 / pow(graph->GetErrorY(i)*graph->GetErrorX(i), 2);
        sum_pythia += graph->GetErrorX(i)*y_graph*weight;
        sum_pythia_weights += weight;


    }


    sum_pythia = sum_pythia/sum_pythia_weights;

    TH1D* hist = new TH1D("hist", "Distribution", 100, -10, 10);

    TF1* func = new TF1("func","pol1", 1.0, 4.5);
    func->SetParameters(0.4891, -0.1212); //bbbar systematic parameters

    TFile* f1 = new TFile("plotting_background_bbbar_central.root","READ");
    TGraphErrors* graph_stat = (TGraphErrors*)f1->Get("g1");

    TGraphErrors* graph_corrected = (TGraphErrors*)f1->Get("g2");

    int n = graph_stat->GetN();

    for(int itr = 0; itr < 1e5; itr++){

        double r = gRandom->Gaus(0, 1);
        double sum = 0;
        double sum_weights = 0;
    
        for (int i = 0; i < n; ++i) {

            if(i == 4 || i == 6) continue;

            double x_val, y_val, new_y;
            graph_stat->GetPoint(i, x_val, y_val);

            if(stat) new_y = y_val + gRandom->Gaus(0, 1) * graph_stat->GetErrorY(i);
            else new_y = y_val + r * func->Eval(x_val)*y_val;

            double weight = 1.0/pow(graph->GetErrorX(i)*graph_stat->GetErrorY(i), 2);
            sum += graph->GetErrorX(i)*new_y*weight;
            sum_weights += weight;


        }

        sum = sum/sum_weights;

        hist->Fill(sum/sum_pythia);

    }

    TCanvas* c2 = new TCanvas();
    hist->Draw();

    TF1* func1 = new TF1("func1","gaus", hist->GetMean()-2*hist->GetStdDev(), hist->GetMean()+2*hist->GetStdDev());
    hist->Fit(func1,"R+");

    TH1D* h1d_bbbar_central = (TH1D*)h1d_bbbar_plot->Clone();
    h1d_bbbar_central->SetName("h1d_bbbar_central");
    h1d_bbbar_central->Scale(func1->GetParameter(1));

    TH1D* h1d_bbbar_lower = (TH1D*)h1d_bbbar_plot->Clone();
    h1d_bbbar_lower->SetName("h1d_bbbar_lower");
    h1d_bbbar_lower->Scale((func1->GetParameter(1) - func1->GetParameter(2)));

    TH1D* h1d_bbbar_upper = (TH1D*)h1d_bbbar_plot->Clone();
    h1d_bbbar_upper->SetName("h1d_bbbar_upper");
    h1d_bbbar_upper->Scale((func1->GetParameter(1) + func1->GetParameter(2)));

    TFile* f2 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_ccbar_bbbar_SoftQCD_inelastic_75925_gencut_with_correction.root","READ");
    TH2D* h2d_delta_phi_pair_mass_bbbar = (TH2D*)f2->Get("h2d_delta_phi_pair_mass_unlike_bbbar");
    h2d_delta_phi_pair_mass_bbbar->Scale(42./75925.);
    h2d_delta_phi_pair_mass_bbbar->Scale(4.8572);

    TH2D* h2d_delta_phi_pair_mass_bbbar_central = (TH2D*)h2d_delta_phi_pair_mass_bbbar->Clone();
    h2d_delta_phi_pair_mass_bbbar_central->SetName("h2d_delta_phi_pair_mass_bbbar_central");
    h2d_delta_phi_pair_mass_bbbar_central->Scale(func1->GetParameter(1));

    TH2D* h2d_delta_phi_pair_mass_bbbar_lower = (TH2D*)h2d_delta_phi_pair_mass_bbbar->Clone();
    h2d_delta_phi_pair_mass_bbbar_lower->SetName("h2d_delta_phi_pair_mass_bbbar_lower");
    h2d_delta_phi_pair_mass_bbbar_lower->Scale(func1->GetParameter(1) - func1->GetParameter(2));


    TH2D* h2d_delta_phi_pair_mass_bbbar_upper = (TH2D*)h2d_delta_phi_pair_mass_bbbar->Clone();
    h2d_delta_phi_pair_mass_bbbar_upper->SetName("h2d_delta_phi_pair_mass_bbbar_upper");
    h2d_delta_phi_pair_mass_bbbar_upper->Scale(func1->GetParameter(1) + func1->GetParameter(2));

    TH2D* h2d_ee_FG_bbbar = (TH2D*)f2->Get("h2d_ee_FG_unlike_bbbar");
    h2d_ee_FG_bbbar->Scale(42./75925.);
    h2d_ee_FG_bbbar->Scale(4.8572);

    TH2D* h2d_ee_FG_bbbar_central = (TH2D*)h2d_ee_FG_bbbar->Clone();
    h2d_ee_FG_bbbar_central->SetName("h2d_ee_FG_bbbar_central");
    h2d_ee_FG_bbbar_central->Scale(func1->GetParameter(1));

    TH2D* h2d_ee_FG_bbbar_lower = (TH2D*)h2d_ee_FG_bbbar->Clone();
    h2d_ee_FG_bbbar_lower->SetName("h2d_ee_FG_bbbar_lower");
    h2d_ee_FG_bbbar_lower->Scale(func1->GetParameter(1) - func1->GetParameter(2));

    TH2D* h2d_ee_FG_bbbar_upper = (TH2D*)h2d_ee_FG_bbbar->Clone();
    h2d_ee_FG_bbbar_upper->SetName("h2d_ee_FG_bbbar_upper");
    h2d_ee_FG_bbbar_upper->Scale(func1->GetParameter(1) + func1->GetParameter(2));

    TH2D* hist_bbbar = (TH2D*)f2->Get("hist_unlike_bbbar");
    hist_bbbar->Scale(42./75925.);
    hist_bbbar->Scale(4.8572);

    TH2D* hist_bbbar_central = (TH2D*)hist_bbbar->Clone();
    hist_bbbar_central->SetName("hist_bbbar_central");
    hist_bbbar_central->Scale(func1->GetParameter(1));

    TH2D* hist_bbbar_lower = (TH2D*)hist_bbbar->Clone();
    hist_bbbar_lower->SetName("hist_bbbar_lower");
    hist_bbbar_lower->Scale(func1->GetParameter(1) - func1->GetParameter(2));

    TH2D* hist_bbbar_upper = (TH2D*)hist_bbbar->Clone();
    hist_bbbar_upper->SetName("hist_bbbar_upper");
    hist_bbbar_upper->Scale(func1->GetParameter(1) + func1->GetParameter(2));

    TGraphAsymmErrors* band = new TGraphAsymmErrors((int)h1d_bbbar_central->GetNbinsX());
    band->SetName("band");

    for (int i = 1; i <= (int)h1d_bbbar_central->GetNbinsX(); ++i) {

        double x = h1d_bbbar_central->GetBinCenter(i);
        double y = h1d_bbbar_central->GetBinContent(i);

        double y_up   = h1d_bbbar_upper->GetBinContent(i);
        double y_down = h1d_bbbar_lower->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band->SetPoint(i - 1, x, y);
        band->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }



    TCanvas* c3 = new TCanvas();
    c3->SetLogy();

    band->SetFillColorAlpha(kBlue, 0.3);
    band->SetLineColor(kBlue);
    band->Draw("A3");
    graph_corrected->Draw("same p");


    TFile* fout = new TFile("normalizing_two_TGraphs_bbbar_sys.root", "RECREATE");
    fout->cd();
    graph->Write();
    graph_corrected->Write();
    h1d_bbbar_central->Write();
    h1d_bbbar_lower->Write();
    h1d_bbbar_upper->Write();
    h2d_delta_phi_pair_mass_bbbar_central->Write();
    h2d_delta_phi_pair_mass_bbbar_lower->Write();
    h2d_delta_phi_pair_mass_bbbar_upper->Write();
    h2d_ee_FG_bbbar_central->Write();
    h2d_ee_FG_bbbar_lower->Write();
    h2d_ee_FG_bbbar_upper->Write();
    hist_bbbar_central->Write();
    hist_bbbar_lower->Write();
    hist_bbbar_upper->Write();
    band->Write();
    fout->Close();


    
}