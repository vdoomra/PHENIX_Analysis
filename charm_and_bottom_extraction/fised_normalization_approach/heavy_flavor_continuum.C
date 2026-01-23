
const int npts = 40;
double x_lohi[npts] = {-999};
double x[npts] = {-999};
double y[npts] = {-999};
double stat_err[npts] = {-999};
double sys_err[npts] = {-999};

const int nbins_mass = 26;
const double mass_range[nbins_mass+1] = {0.87, 0.95, 1.05, 1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.05, 3.10, 3.15, 3.20, 3.30, 3.45, 3.60, 3.70, 3.80, 4.00, 4.25, 4.5};

const int nbins_mass_hf = 8;
const double mass_range_hf[nbins_mass_hf+1] = {0.87, 1.05, 1.45, 1.75, 2.05, 2.40, 2.80, 3.80, 4.50};

const double nevents = 4.216682807*pow(10,11)*(0.79);

double combinedFunction(double *x, double *par) {

    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    double fn0_val = par[3] + par[4] * x[0] + par[5]*x[0]*x[0];
    double feoverp_val = par[6] + par[7]*x[0] + par[8]*x[0]*x[0];

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}

void heavy_flavor_continuum(const int stat){

    TCanvas *c1 = new TCanvas();

    TF1* fcombined = new TF1("fcombined", combinedFunction, 0.0, 5.5, 9);
    fcombined->SetParameters(0.6885, -0.3019, 0.03388, 0.2268, -0.05276, 0.00609, 0.3725, -0.02848, -0.01353);

    TFile* fdata = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_FG_phi_redef_n0_2_with_correction.root","READ");
    TH1D* h1d_unlike_data_temp = (TH1D*)fdata->Get("h1d_ee_FG_unlike");
    h1d_unlike_data_temp->Scale(1./nevents);
    h1d_unlike_data_temp->Scale(1.0726*0.545*42);

    TH1D* h1d_unlike_data = new TH1D("h1d_unlike_data","Mass Distribution", nbins_mass, mass_range);
    h1d_unlike_data->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_unlike_data_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_unlike_data_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_unlike_data_temp->GetBinContent(i)*h1d_unlike_data_temp->GetBinWidth(i);
            bin_error += pow(h1d_unlike_data_temp->GetBinError(i)*h1d_unlike_data_temp->GetBinWidth(i), 2);

        }

        h1d_unlike_data->SetBinContent(ibin+1, bin_content/h1d_unlike_data->GetBinWidth(ibin+1));
        h1d_unlike_data->SetBinError(ibin+1, sqrt(bin_error)/h1d_unlike_data->GetBinWidth(ibin+1));


    }

    TFile* fjpsi = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_jpsi_2k_gencut_with_correction.root","READ");
    TH1D* h1d_unlike_jpsi_temp = (TH1D*)fjpsi->Get("h1d_ee_FG_unlike");
    h1d_unlike_jpsi_temp->SetName("h1d_unlike_jpsi_temp");
    h1d_unlike_jpsi_temp->Scale(42.*0.05971);
    h1d_unlike_jpsi_temp->Scale(1./(2000.*25000.));

    TH1D* h1d_unlike_jpsi = new TH1D("h1d_unlike_jpsi","Mass Distribution", nbins_mass, mass_range);
    h1d_unlike_jpsi->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_unlike_jpsi_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_unlike_jpsi_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_unlike_jpsi_temp->GetBinContent(i)*h1d_unlike_jpsi_temp->GetBinWidth(i);
            bin_error += pow(h1d_unlike_jpsi_temp->GetBinError(i)*h1d_unlike_jpsi_temp->GetBinWidth(i), 2);

        }

        h1d_unlike_jpsi->SetBinContent(ibin+1, bin_content/h1d_unlike_jpsi->GetBinWidth(ibin+1));
        h1d_unlike_jpsi->SetBinError(ibin+1, 0.0);


    }

    TFile* fpsip = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_psip_2k_gencut_with_correction.root","READ");
    TH1D* h1d_unlike_psip_temp = (TH1D*)fpsip->Get("h1d_ee_FG_unlike");
    h1d_unlike_psip_temp->SetName("h1d_unlike_psip_temp");
    h1d_unlike_psip_temp->Scale(42.*(7.93e-3));
    h1d_unlike_psip_temp->Scale(1./(2000.*25000.));

    TH1D* h1d_unlike_psip = new TH1D("h1d_unlike_psip","Mass Distribution", nbins_mass, mass_range);
    h1d_unlike_psip->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_unlike_psip_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_unlike_psip_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_unlike_psip_temp->GetBinContent(i)*h1d_unlike_psip_temp->GetBinWidth(i);
            bin_error += pow(h1d_unlike_psip_temp->GetBinError(i)*h1d_unlike_psip_temp->GetBinWidth(i), 2);

        }

        h1d_unlike_psip->SetBinContent(ibin+1, bin_content/h1d_unlike_psip->GetBinWidth(ibin+1));
        h1d_unlike_psip->SetBinError(ibin+1, 0.0);


    }

    TFile* fphi = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_phi_ee_2k_gencut_with_correction.root","READ");
    TH1D* h1d_unlike_phi_temp = (TH1D*)fphi->Get("h1d_ee_FG_unlike");
    h1d_unlike_phi_temp->SetName("h1d_unlike_phi_temp");
    h1d_unlike_phi_temp->Scale(42.*0.000297);
    h1d_unlike_phi_temp->Scale(1./(2000.*25000.));

    TH1D* h1d_unlike_phi = new TH1D("h1d_unlike_phi","Mass Distribution", nbins_mass, mass_range);
    h1d_unlike_phi->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_unlike_phi_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_unlike_phi_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_unlike_phi_temp->GetBinContent(i)*h1d_unlike_phi_temp->GetBinWidth(i);
            bin_error += pow(h1d_unlike_phi_temp->GetBinError(i)*h1d_unlike_phi_temp->GetBinWidth(i), 2);

        }

        h1d_unlike_phi->SetBinContent(ibin+1, bin_content/h1d_unlike_phi->GetBinWidth(ibin+1));
        h1d_unlike_phi->SetBinError(ibin+1, 0.0);


    }

    TFile* fBG = new TFile("plotting_BG_with_uncertainties.root","READ");
    TH1D* h1d_unlike_BG = (TH1D*)fBG->Get("h1d_mlp_dataBG");
    h1d_unlike_BG->SetName("h1d_unlike_BG");

    TH1D* h1d_signal_temp = (TH1D*)h1d_unlike_data->Clone();
    h1d_signal_temp->Reset("ICESM");
    h1d_signal_temp->SetName("h1d_signal_temp");

    int nbinsX = h1d_signal_temp->GetNbinsX();

    for(int ibin=1; ibin < nbinsX+1; ibin++){

        double bin_content = h1d_unlike_data->GetBinContent(ibin) - h1d_unlike_jpsi->GetBinContent(ibin) - h1d_unlike_psip->GetBinContent(ibin) - h1d_unlike_BG->GetBinContent(ibin) - h1d_unlike_phi->GetBinContent(ibin);
        double bin_error = sqrt( pow(h1d_unlike_data->GetBinError(ibin), 2) + pow(h1d_unlike_BG->GetBinError(ibin), 2) );

        h1d_signal_temp->SetBinContent(ibin, bin_content);
        h1d_signal_temp->SetBinError(ibin, bin_error);


    }

    TCanvas* c2 = new TCanvas();
    h1d_signal_temp->Draw();

    TH1D* h1d_signal_stat = new TH1D("h1d_signal_stat","Mass Distribution", nbins_mass_hf, mass_range_hf);
    h1d_signal_stat->Sumw2();

    TH1D* h1d_signal_sys = new TH1D("h1d_signal_sys","Mass Distribution", nbins_mass_hf, mass_range_hf);
    h1d_signal_sys->Sumw2();

    for(int ibin=0; ibin < nbins_mass_hf; ibin++){

        int bin_lo = h1d_signal_temp->FindBin(mass_range_hf[ibin] + 0.01);
        int bin_hi = h1d_signal_temp->FindBin(mass_range_hf[ibin+1] - 0.01);

        double bin_center = (mass_range_hf[ibin] + mass_range_hf[ibin+1])/2;

        double bin_content = 0; 
        double bin_error_stat = 0;  
        double bin_error_sys = 0;
        double bin_content_bg = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_signal_temp->GetBinContent(i)*h1d_signal_temp->GetBinWidth(i);
            bin_error_stat += pow(h1d_signal_temp->GetBinError(i)*h1d_signal_temp->GetBinWidth(i), 2);

            bin_content_bg += h1d_unlike_BG->GetBinContent(i);
            

        }

        double sys_signal = bin_content*fcombined->Eval(bin_center);
        double sys_bg = bin_content_bg*0.05212;
            
        bin_error_sys = sqrt( pow(sys_signal, 2) + pow(sys_bg, 2));

        h1d_signal_stat->SetBinContent(ibin+1, bin_content/h1d_signal_stat->GetBinWidth(ibin+1));
        h1d_signal_sys->SetBinContent(ibin+1, bin_content/h1d_signal_sys->GetBinWidth(ibin+1));
        h1d_signal_stat->SetBinError(ibin+1, sqrt(bin_error_stat)/h1d_signal_stat->GetBinWidth(ibin+1));
        h1d_signal_sys->SetBinError(ibin+1, bin_error_sys);


    }

    double h1d_signal_err_stat, h1d_signal_err_sys;
    double integral_stat = h1d_signal_stat->IntegralAndError(1, h1d_signal_stat->GetNbinsX(), h1d_signal_err_stat, "width");
    double integral_sys = h1d_signal_sys->IntegralAndError(1, h1d_signal_sys->GetNbinsX(), h1d_signal_err_sys, "width");

    cout << integral_stat << "\t" << h1d_signal_err_stat << "\t" << h1d_signal_err_sys << endl;


    int nbins = h1d_signal_stat->GetNbinsX();
    int counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = h1d_signal_stat->GetBinCenter(ibin);
        y[counter] = h1d_signal_stat->GetBinContent(ibin);
        x_lohi[counter] = 0.05;

        stat_err[counter] = h1d_signal_stat->GetBinError(ibin);
        sys_err[counter] = h1d_signal_sys->GetBinError(ibin);

        counter++;


    }

    TGraphErrors* graph_stat = new TGraphErrors(counter, x, y, 0, stat_err);
    graph_stat->SetMarkerStyle(20);
    graph_stat->SetMarkerSize(0.9);
    graph_stat->SetMarkerColor(kBlue-4);
    graph_stat->SetLineColor(kBlue-4);
    graph_stat->SetName("graph_stat");

    TGraphErrors* graph_sys = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    graph_sys->SetFillColorAlpha(kBlue-4, 0.3);
    graph_sys->SetLineColor(kBlue-4);
    graph_sys->SetName("graph_sys");

    c1->cd();
    graph_stat->Draw("AP");
    graph_sys->Draw("e2same");

    TFile* fsim = new TFile("/Users/vassu/Desktop/new_files/data_files/producing_templates_without_20percent_reduction.root","READ");
    TH1D* h1d_drell_yan = (TH1D*)fsim->Get("h1d_ee_FG_drellyan");
    h1d_drell_yan->SetName("h1d_drell_yan");

    TH1D* h1d_ccbar = (TH1D*)fsim->Get("h1d_ee_FG_ccbar");
    h1d_ccbar->SetName("h1d_ccbar");

    TH1D* h1d_bbbar = (TH1D*)fsim->Get("h1d_ee_FG_bbbar");
    h1d_bbbar->SetName("h1d_bbbar");

    TH1D* hist_total_temp = (TH1D*)h1d_ccbar->Clone();
    hist_total_temp->SetName("hist_total_temp");
    hist_total_temp->Add(h1d_bbbar);
    hist_total_temp->Add(h1d_drell_yan);

    TH1D* hist_total = new TH1D("hist_total","Mass Distribution", nbins_mass_hf, mass_range_hf);
    hist_total->Sumw2();

    for(int ibin=0; ibin < nbins_mass_hf; ibin++){

        int bin_lo = hist_total_temp->FindBin(mass_range_hf[ibin] + 0.001);
        int bin_hi = hist_total_temp->FindBin(mass_range_hf[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += hist_total_temp->GetBinContent(i)*hist_total_temp->GetBinWidth(i);
            bin_error += pow(hist_total_temp->GetBinError(i)*hist_total_temp->GetBinWidth(i), 2);

        }

        hist_total->SetBinContent(ibin+1, bin_content/hist_total->GetBinWidth(ibin+1));
        hist_total->SetBinError(ibin+1, 0.0);


    }

    hist_total->SetLineColor(kBlue);
    hist_total->SetLineStyle(3);
    hist_total->Draw("psame");

   
    double scale, scale_error;

    if(stat){

        double h1d_signal_error, hist_total_error;
        double S = h1d_signal_stat->IntegralAndError(1, h1d_signal_stat->GetNbinsX(), h1d_signal_error, "width");
        double B = hist_total->IntegralAndError(1, hist_total->GetNbinsX(), hist_total_error, "width");

        scale = S/B;
        scale_error = scale*sqrt(pow(h1d_signal_error/S, 2) + pow(hist_total_error/B, 2));

    }

    else{

        TH1D* hist = new TH1D("hist", "Distribution", 200, -10, 20);

        for(int itr = 0; itr < 2e5; itr++){

            double r = gRandom->Gaus(0, 1);

            TH1D* h1d_signal_iterative = new TH1D(Form("h1d_signal_iterative_%d", itr),"Mass Distribution", nbins_mass_hf, mass_range_hf);
            h1d_signal_iterative->Sumw2();
        
            for (int i = 1; i < h1d_signal_iterative->GetNbinsX()+1; ++i) {
    
                double new_y = (h1d_signal_sys->GetBinContent(i) + r * h1d_signal_sys->GetBinError(i));
                h1d_signal_iterative->SetBinContent(i, new_y);
                h1d_signal_iterative->SetBinError(i, h1d_signal_stat->GetBinError(i));
    
            }
    
            double scale = h1d_signal_iterative->Integral("width")/hist_total->Integral("width");
    
            hist->Fill(scale);
    
        }

        TCanvas* c2 = new TCanvas();
        c2->cd();

        TF1* func = new TF1("func","gaus", hist->GetMean()-2*hist->GetStdDev(), hist->GetMean()+2*hist->GetStdDev());
        hist->Fit(func,"R+");

        scale = func->GetParameter(1);
        scale_error = func->GetParameter(2);

    }

    cout << scale << "\t" << scale_error << endl;

    TH1D* hist_total_central = (TH1D*)hist_total->Clone();
    hist_total_central->SetName("hist_total_central");
    hist_total_central->Scale(scale);

    TH1D* hist_total_lower = (TH1D*)hist_total->Clone();
    hist_total_lower->SetName("hist_total_lower");
    hist_total_lower->Scale(scale-scale_error);

    TH1D* hist_total_upper = (TH1D*)hist_total->Clone();
    hist_total_upper->SetName("hist_total_upper");
    hist_total_upper->Scale(scale+scale_error);

    TH1D* hist_total_dataBinning_central = (TH1D*)hist_total_temp->Clone();
    hist_total_dataBinning_central->SetName("hist_total_dataBinning_central");
    hist_total_dataBinning_central->Scale(scale);

    TH1D* hist_total_dataBinning_lower = (TH1D*)hist_total_temp->Clone();
    hist_total_dataBinning_lower->SetName("hist_total_dataBinning_lower");
    hist_total_dataBinning_lower->Scale(scale-scale_error);

    TH1D* hist_total_dataBinning_upper = (TH1D*)hist_total_temp->Clone();
    hist_total_dataBinning_upper->SetName("hist_total_dataBinning_upper");
    hist_total_dataBinning_upper->Scale(scale+scale_error);

    TH1D* hist_ccbar_dataBinning_central = (TH1D*)h1d_ccbar->Clone();
    hist_ccbar_dataBinning_central->SetName("hist_ccbar_dataBinning_central");
    hist_ccbar_dataBinning_central->Scale(scale);

    TH1D* hist_ccbar_dataBinning_lower = (TH1D*)h1d_ccbar->Clone();
    hist_ccbar_dataBinning_lower->SetName("hist_ccbar_dataBinning_lower");
    hist_ccbar_dataBinning_lower->Scale(scale-scale_error);

    TH1D* hist_ccbar_dataBinning_upper = (TH1D*)h1d_ccbar->Clone();
    hist_ccbar_dataBinning_upper->SetName("hist_ccbar_dataBinning_upper");
    hist_ccbar_dataBinning_upper->Scale(scale+scale_error);

    TH1D* hist_bbbar_dataBinning_central = (TH1D*)h1d_bbbar->Clone();
    hist_bbbar_dataBinning_central->SetName("hist_bbbar_dataBinning_central");
    hist_bbbar_dataBinning_central->Scale(scale);

    TH1D* hist_bbbar_dataBinning_lower = (TH1D*)h1d_bbbar->Clone();
    hist_bbbar_dataBinning_lower->SetName("hist_bbbar_dataBinning_lower");
    hist_bbbar_dataBinning_lower->Scale(scale-scale_error);

    TH1D* hist_bbbar_dataBinning_upper = (TH1D*)h1d_bbbar->Clone();
    hist_bbbar_dataBinning_upper->SetName("hist_bbbar_dataBinning_upper");
    hist_bbbar_dataBinning_upper->Scale(scale+scale_error);

    TH1D* hist_drell_yan_dataBinning_central = (TH1D*)h1d_drell_yan->Clone();
    hist_drell_yan_dataBinning_central->SetName("hist_drell_yan_dataBinning_central");
    hist_drell_yan_dataBinning_central->Scale(scale);

    TH1D* hist_drell_yan_dataBinning_lower = (TH1D*)h1d_drell_yan->Clone();
    hist_drell_yan_dataBinning_lower->SetName("hist_drell_yan_dataBinning_lower");
    hist_drell_yan_dataBinning_lower->Scale(scale-scale_error);

    TH1D* hist_drell_yan_dataBinning_upper = (TH1D*)h1d_drell_yan->Clone();
    hist_drell_yan_dataBinning_upper->SetName("hist_drell_yan_dataBinning_upper");
    hist_drell_yan_dataBinning_upper->Scale(scale+scale_error);

    TGraphAsymmErrors* band = new TGraphAsymmErrors((int)hist_total_central->GetNbinsX());
    band->SetName("band");

    for (int i = 1; i <= (int)hist_total_central->GetNbinsX(); ++i) {

        double x = hist_total_central->GetBinCenter(i);
        double y = hist_total_central->GetBinContent(i);

        double y_up   = hist_total_upper->GetBinContent(i);
        double y_down = hist_total_lower->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band->SetPoint(i - 1, x, y);
        band->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }


    c1->cd();
    band->SetFillColorAlpha(kBlue, 0.3);
    band->SetLineColor(kBlue);
    band->Draw("3same");

    TFile* f3 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_ccbar_bbbar_SoftQCD_inelastic_75925_gencut_with_correction.root","READ");

    TH2D* h2d_delta_phi_pair_mass_total = (TH2D*)f3->Get("h2d_delta_phi_pair_mass_unlike_total");
    h2d_delta_phi_pair_mass_total->Scale(42./75925.);
    h2d_delta_phi_pair_mass_total->Scale(4.8572);

    TH2D* h2d_delta_phi_pair_mass_total_central = (TH2D*)h2d_delta_phi_pair_mass_total->Clone();
    h2d_delta_phi_pair_mass_total_central->SetName("h2d_delta_phi_pair_mass_total_central");
    h2d_delta_phi_pair_mass_total_central->Scale(scale);

    TH2D* h2d_delta_phi_pair_mass_total_lower = (TH2D*)h2d_delta_phi_pair_mass_total->Clone();
    h2d_delta_phi_pair_mass_total_lower->SetName("h2d_delta_phi_pair_mass_total_lower");
    h2d_delta_phi_pair_mass_total_lower->Scale(scale - scale_error);

    TH2D* h2d_delta_phi_pair_mass_total_upper = (TH2D*)h2d_delta_phi_pair_mass_total->Clone();
    h2d_delta_phi_pair_mass_total_upper->SetName("h2d_delta_phi_pair_mass_total_upper");
    h2d_delta_phi_pair_mass_total_upper->Scale(scale + scale_error);

    TH2D* h2d_delta_phi_pair_mass_ccbar = (TH2D*)f3->Get("h2d_delta_phi_pair_mass_unlike_ccbar");
    h2d_delta_phi_pair_mass_ccbar->Scale(42./75925.);
    h2d_delta_phi_pair_mass_ccbar->Scale(4.8572);

    TH2D* h2d_delta_phi_pair_mass_ccbar_central = (TH2D*)h2d_delta_phi_pair_mass_ccbar->Clone();
    h2d_delta_phi_pair_mass_ccbar_central->SetName("h2d_delta_phi_pair_mass_ccbar_central");
    h2d_delta_phi_pair_mass_ccbar_central->Scale(scale);

    TH2D* h2d_delta_phi_pair_mass_ccbar_lower = (TH2D*)h2d_delta_phi_pair_mass_ccbar->Clone();
    h2d_delta_phi_pair_mass_ccbar_lower->SetName("h2d_delta_phi_pair_mass_ccbar_lower");
    h2d_delta_phi_pair_mass_ccbar_lower->Scale(scale - scale_error);

    TH2D* h2d_delta_phi_pair_mass_ccbar_upper = (TH2D*)h2d_delta_phi_pair_mass_ccbar->Clone();
    h2d_delta_phi_pair_mass_ccbar_upper->SetName("h2d_delta_phi_pair_mass_ccbar_upper");
    h2d_delta_phi_pair_mass_ccbar_upper->Scale(scale + scale_error);

    TH2D* h2d_delta_phi_pair_mass_bbbar = (TH2D*)f3->Get("h2d_delta_phi_pair_mass_unlike_bbbar");
    h2d_delta_phi_pair_mass_bbbar->Scale(42./75925.);
    h2d_delta_phi_pair_mass_bbbar->Scale(4.8572);

    TH2D* h2d_delta_phi_pair_mass_bbbar_central = (TH2D*)h2d_delta_phi_pair_mass_bbbar->Clone();
    h2d_delta_phi_pair_mass_bbbar_central->SetName("h2d_delta_phi_pair_mass_bbbar_central");
    h2d_delta_phi_pair_mass_bbbar_central->Scale(scale);

    TH2D* h2d_delta_phi_pair_mass_bbbar_lower = (TH2D*)h2d_delta_phi_pair_mass_bbbar->Clone();
    h2d_delta_phi_pair_mass_bbbar_lower->SetName("h2d_delta_phi_pair_mass_bbbar_lower");
    h2d_delta_phi_pair_mass_bbbar_lower->Scale(scale - scale_error);

    TH2D* h2d_delta_phi_pair_mass_bbbar_upper = (TH2D*)h2d_delta_phi_pair_mass_bbbar->Clone();
    h2d_delta_phi_pair_mass_bbbar_upper->SetName("h2d_delta_phi_pair_mass_bbbar_upper");
    h2d_delta_phi_pair_mass_bbbar_upper->Scale(scale + scale_error);

    TH2D* h2d_ee_FG_total = (TH2D*)f3->Get("h2d_ee_FG_unlike_total");
    h2d_ee_FG_total->Scale(42./75925.);
    h2d_ee_FG_total->Scale(4.8572);

    TH2D* h2d_ee_FG_total_central = (TH2D*)h2d_ee_FG_total->Clone();
    h2d_ee_FG_total_central->SetName("h2d_ee_FG_total_central");
    h2d_ee_FG_total_central->Scale(scale);

    TH2D* h2d_ee_FG_total_lower = (TH2D*)h2d_ee_FG_total->Clone();
    h2d_ee_FG_total_lower->SetName("h2d_ee_FG_total_lower");
    h2d_ee_FG_total_lower->Scale(scale - scale_error);

    TH2D* h2d_ee_FG_total_upper = (TH2D*)h2d_ee_FG_total->Clone();
    h2d_ee_FG_total_upper->SetName("h2d_ee_FG_total_upper");
    h2d_ee_FG_total_upper->Scale(scale + scale_error);

    TH2D* h2d_ee_FG_ccbar = (TH2D*)f3->Get("h2d_ee_FG_unlike_ccbar");
    h2d_ee_FG_ccbar->Scale(42./75925.);
    h2d_ee_FG_ccbar->Scale(4.8572);

    TH2D* h2d_ee_FG_ccbar_central = (TH2D*)h2d_ee_FG_ccbar->Clone();
    h2d_ee_FG_ccbar_central->SetName("h2d_ee_FG_ccbar_central");
    h2d_ee_FG_ccbar_central->Scale(scale);

    TH2D* h2d_ee_FG_ccbar_lower = (TH2D*)h2d_ee_FG_ccbar->Clone();
    h2d_ee_FG_ccbar_lower->SetName("h2d_ee_FG_ccbar_lower");
    h2d_ee_FG_ccbar_lower->Scale(scale - scale_error);

    TH2D* h2d_ee_FG_ccbar_upper = (TH2D*)h2d_ee_FG_ccbar->Clone();
    h2d_ee_FG_ccbar_upper->SetName("h2d_ee_FG_ccbar_upper");
    h2d_ee_FG_ccbar_upper->Scale(scale + scale_error);

    TH2D* h2d_ee_FG_bbbar = (TH2D*)f3->Get("h2d_ee_FG_unlike_bbbar");
    h2d_ee_FG_bbbar->Scale(42./75925.);
    h2d_ee_FG_bbbar->Scale(4.8572);

    TH2D* h2d_ee_FG_bbbar_central = (TH2D*)h2d_ee_FG_bbbar->Clone();
    h2d_ee_FG_bbbar_central->SetName("h2d_ee_FG_bbbar_central");
    h2d_ee_FG_bbbar_central->Scale(scale);

    TH2D* h2d_ee_FG_bbbar_lower = (TH2D*)h2d_ee_FG_bbbar->Clone();
    h2d_ee_FG_bbbar_lower->SetName("h2d_ee_FG_bbbar_lower");
    h2d_ee_FG_bbbar_lower->Scale(scale - scale_error);

    TH2D* h2d_ee_FG_bbbar_upper = (TH2D*)h2d_ee_FG_bbbar->Clone();
    h2d_ee_FG_bbbar_upper->SetName("h2d_ee_FG_bbbar_upper");
    h2d_ee_FG_bbbar_upper->Scale(scale + scale_error);

    TH2D* hist_dca_total = (TH2D*)f3->Get("hist_unlike_total");
    hist_dca_total->Scale(42./75925.);
    hist_dca_total->Scale(4.8572);

    TH2D* hist_dca_total_central = (TH2D*)hist_dca_total->Clone();
    hist_dca_total_central->SetName("hist_dca_total_central");
    hist_dca_total_central->Scale(scale);

    TH2D* hist_dca_total_lower = (TH2D*)hist_dca_total->Clone();
    hist_dca_total_lower->SetName("hist_dca_total_lower");
    hist_dca_total_lower->Scale(scale - scale_error);

    TH2D* hist_dca_total_upper = (TH2D*)hist_dca_total->Clone();
    hist_dca_total_upper->SetName("hist_dca_total_upper");
    hist_dca_total_upper->Scale(scale + scale_error);

    TH2D* hist_dca_ccbar = (TH2D*)f3->Get("hist_unlike_ccbar");
    hist_dca_ccbar->Scale(42./75925.);
    hist_dca_ccbar->Scale(4.8572);

    TH2D* hist_dca_ccbar_central = (TH2D*)hist_dca_ccbar->Clone();
    hist_dca_ccbar_central->SetName("hist_dca_ccbar_central");
    hist_dca_ccbar_central->Scale(scale);

    TH2D* hist_dca_ccbar_lower = (TH2D*)hist_dca_ccbar->Clone();
    hist_dca_ccbar_lower->SetName("hist_dca_ccbar_lower");
    hist_dca_ccbar_lower->Scale(scale - scale_error);

    TH2D* hist_dca_ccbar_upper = (TH2D*)hist_dca_ccbar->Clone();
    hist_dca_ccbar_upper->SetName("hist_dca_ccbar_upper");
    hist_dca_ccbar_upper->Scale(scale + scale_error);

    TH2D* hist_dca_bbbar = (TH2D*)f3->Get("hist_unlike_bbbar");
    hist_dca_bbbar->Scale(42./75925.);
    hist_dca_bbbar->Scale(4.8572);

    TH2D* hist_dca_bbbar_central = (TH2D*)hist_dca_bbbar->Clone();
    hist_dca_bbbar_central->SetName("hist_dca_bbbar_central");
    hist_dca_bbbar_central->Scale(scale);

    TH2D* hist_dca_bbbar_lower = (TH2D*)hist_dca_bbbar->Clone();
    hist_dca_bbbar_lower->SetName("hist_dca_bbbar_lower");
    hist_dca_bbbar_lower->Scale(scale - scale_error);

    TH2D* hist_dca_bbbar_upper = (TH2D*)hist_dca_bbbar->Clone();
    hist_dca_bbbar_upper->SetName("hist_dca_bbbar_upper");
    hist_dca_bbbar_upper->Scale(scale + scale_error);

    TFile* fout = new TFile("dummy.root","RECREATE");
    fout->cd();
    graph_stat->Write();
    graph_sys->Write();

    hist_total_central->Write();
    hist_total_lower->Write();
    hist_total_upper->Write();

    hist_total_dataBinning_central->Write();
    hist_total_dataBinning_lower->Write();
    hist_total_dataBinning_upper->Write();

    hist_ccbar_dataBinning_central->Write();
    hist_ccbar_dataBinning_lower->Write();
    hist_ccbar_dataBinning_upper->Write();

    hist_bbbar_dataBinning_central->Write();
    hist_bbbar_dataBinning_lower->Write();
    hist_bbbar_dataBinning_upper->Write();

    hist_drell_yan_dataBinning_central->Write();
    hist_drell_yan_dataBinning_lower->Write();
    hist_drell_yan_dataBinning_upper->Write();

    h2d_delta_phi_pair_mass_total_central->Write();
    h2d_delta_phi_pair_mass_total_lower->Write();
    h2d_delta_phi_pair_mass_total_upper->Write();

    h2d_delta_phi_pair_mass_ccbar_central->Write();
    h2d_delta_phi_pair_mass_ccbar_lower->Write();
    h2d_delta_phi_pair_mass_ccbar_upper->Write();

    h2d_delta_phi_pair_mass_bbbar_central->Write();
    h2d_delta_phi_pair_mass_bbbar_lower->Write();
    h2d_delta_phi_pair_mass_bbbar_upper->Write();

    h2d_ee_FG_total_central->Write();
    h2d_ee_FG_total_lower->Write();
    h2d_ee_FG_total_upper->Write();

    h2d_ee_FG_ccbar_central->Write();
    h2d_ee_FG_ccbar_lower->Write();
    h2d_ee_FG_ccbar_upper->Write();

    h2d_ee_FG_bbbar_central->Write();
    h2d_ee_FG_bbbar_lower->Write();
    h2d_ee_FG_bbbar_upper->Write();

    hist_dca_total_central->Write();
    hist_dca_total_lower->Write();
    hist_dca_total_upper->Write();

    hist_dca_ccbar_central->Write();
    hist_dca_ccbar_lower->Write();
    hist_dca_ccbar_upper->Write();

    hist_dca_bbbar_central->Write();
    hist_dca_bbbar_lower->Write();
    hist_dca_bbbar_upper->Write();


    band->Write();
    fout->Close();

}