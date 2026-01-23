
const int npts = 200;
double x_lohi[npts] = {-999};
double x[npts] = {-999};
double y[npts] = {-999};
double stat_err[npts] = {-999};
double sys_err[npts] = {-999};

const int nbins_mass = 26;
const double mass_range[nbins_mass+1] = {0.87, 0.95, 1.05, 1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.05, 3.10, 3.15, 3.20, 3.30, 3.45, 3.60, 3.70, 3.80, 4.00, 4.25, 4.5};

const int nbins_mass_hf = 9;
const double mass_range_hf[nbins_mass_hf+1] = {0.87, 1.15, 1.45, 1.75, 2.05, 2.40, 2.80, 3.15, 3.80, 4.50};

const double nevents = 4.216682807*pow(10,11)*(0.79);

double combinedFunction(double *x, double *par) {
    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    if(x[0] > 3.50) fERT_val = 0.0;
    double fn0_val = par[3];
    double feoverp_val = par[4] + par[5]*x[0] + par[6]*x[0]*x[0];
    if(x[0] > 3.20) feoverp_val = 0.03262;

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}

double combinedFunction_jpsi(double *x, double *par) {

    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    double fn0_val = par[3] + par[4] * x[0] + par[5]*x[0]*x[0];
    double feoverp_val = par[6] + par[7]*x[0] + par[8]*x[0]*x[0];

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}

double combinedFunction_psip(double *x, double *par) {
    
    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    double fn0_val = par[3] + par[4] * x[0] + par[5]*x[0]*x[0];
    double feoverp_val = par[6] + par[7]*x[0] + par[8]*x[0]*x[0];

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}


void heavy_flavor_continuum_v2(){

    TCanvas *c1 = new TCanvas("c", "Unified Two Panel Plot", 800, 800);

    TPad *pad1 = new TPad("pad1", "Top pad", 0.0, 0.35, 1.0, 1.0);
    TPad *pad2 = new TPad("pad2", "Bottom pad", 0.0, 0.0, 1.0, 0.35);


    pad1->SetBottomMargin(0.0); 
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.35);

    pad1->SetTickx(0);

    pad1->Draw();
    pad2->Draw();

    TF1* fcombined = new TF1("fcombined", combinedFunction, 0.0, 5.5, 7);
    fcombined->SetParameters(0.102, -0.0022, -0.00492, 0.06424, 0.06622, 0.08956, -0.03127);

    TF1* fcombined_jpsi = new TF1("fcombined_jpsi", combinedFunction_jpsi, 0.0, 5.5, 9);
    fcombined_jpsi->SetParameters(0.0973, 0.00915, -0.00925, 0.289, -0.1783, 0.0318, 0.412, -0.268, 0.0465);

    TF1* fcombined_psip = new TF1("fcombined_psip", combinedFunction_psip, 0.0, 5.5, 9);
    fcombined_psip->SetParameters(0.091, -0.0108, -0.0025, 0.306, -0.1561, 0.0234, 0.4744, -0.272, 0.0404);

    TFile* fdata = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_FG_phi_redef_n0_2.root","READ");
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

    TFile* fjpsi = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_jpsi_2k_gencut.root","READ");
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

    TFile* fpsip = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_psip_2k_gencut.root","READ");
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

    TFile* fBG = new TFile("/Users/vassu/Desktop/new_files/new_approach_v2/plotting_BG_with_uncertainties.root","READ");
    TH1D* h1d_unlike_BG = (TH1D*)fBG->Get("h1d_mlp_dataBG");
    h1d_unlike_BG->SetName("h1d_unlike_BG");

    TH1D* h1d_signal_temp = (TH1D*)h1d_unlike_data->Clone();
    h1d_signal_temp->Reset("ICESM");
    h1d_signal_temp->SetName("h1d_signal_temp");

    int nbinsX = h1d_signal_temp->GetNbinsX();

    for(int ibin=1; ibin < nbinsX+1; ibin++){

        double bin_content = h1d_unlike_data->GetBinContent(ibin) - h1d_unlike_jpsi->GetBinContent(ibin) - h1d_unlike_psip->GetBinContent(ibin) - h1d_unlike_BG->GetBinContent(ibin);
        double bin_error = sqrt( pow(h1d_unlike_data->GetBinError(ibin), 2) + pow(h1d_unlike_BG->GetBinError(ibin), 2) );

        h1d_signal_temp->SetBinContent(ibin, bin_content);
        h1d_signal_temp->SetBinError(ibin, bin_error);


    }

    TH1D* h1d_signal_stat = new TH1D("h1d_signal_stat","Mass Distribution", nbins_mass_hf, mass_range_hf);
    h1d_signal_stat->Sumw2();

    TH1D* h1d_signal_sys = new TH1D("h1d_signal_sys","Mass Distribution", nbins_mass_hf, mass_range_hf);
    h1d_signal_sys->Sumw2();

    for(int ibin=0; ibin < nbins_mass_hf; ibin++){

        int bin_lo = h1d_signal_temp->FindBin(mass_range_hf[ibin] + 0.01);
        int bin_hi = h1d_signal_temp->FindBin(mass_range_hf[ibin+1] - 0.01);

        double bin_center = (mass_range_hf[ibin] + mass_range_hf[ibin+1])/2;

        double bin_content = 0; double bin_content_fg = 0;  double bin_content_bg = 0; double bin_content_jpsi = 0; 
        double bin_content_psip = 0; double bin_error_stat = 0;  double bin_error_sys = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_signal_temp->GetBinContent(i)*h1d_signal_temp->GetBinWidth(i);
            bin_error_stat += pow(h1d_signal_temp->GetBinError(i)*h1d_signal_temp->GetBinWidth(i), 2);

            bin_content_fg += h1d_unlike_data->GetBinContent(i)*h1d_unlike_data->GetBinWidth(i);
            bin_content_bg += h1d_unlike_BG->GetBinContent(i)*h1d_unlike_BG->GetBinWidth(i);
            bin_content_jpsi += h1d_unlike_jpsi->GetBinContent(i)*h1d_unlike_jpsi->GetBinWidth(i);
            bin_content_psip += h1d_unlike_psip->GetBinContent(i)*h1d_unlike_psip->GetBinWidth(i);

        }

        double sys_fg = bin_content_fg*fcombined->Eval(bin_center);
        double sys_jpsi = fcombined_jpsi->Eval(bin_center)*bin_content_jpsi;
        double sys_psip = fcombined_psip->Eval(bin_center)*bin_content_psip;
        double sys_bg = bin_content_bg*0.05212;
            
        bin_error_sys = sqrt( pow(sys_fg, 2) + pow(sys_bg, 2) + pow(sys_jpsi, 2) + pow(sys_psip, 2));

        h1d_signal_stat->SetBinContent(ibin+1, bin_content/h1d_signal_stat->GetBinWidth(ibin+1));
        h1d_signal_sys->SetBinContent(ibin+1, bin_content/h1d_signal_sys->GetBinWidth(ibin+1));
        h1d_signal_stat->SetBinError(ibin+1, sqrt(bin_error_stat)/h1d_signal_stat->GetBinWidth(ibin+1));
        h1d_signal_sys->SetBinError(ibin+1, bin_error_sys/h1d_signal_sys->GetBinWidth(ibin+1));


    }

    int nbins = h1d_signal_stat->GetNbinsX();
    int counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = h1d_signal_stat->GetBinCenter(ibin);
        y[counter] = h1d_signal_stat->GetBinContent(ibin);
        x_lohi[counter] = h1d_signal_stat->GetBinWidth(ibin)/5;

        stat_err[counter] = h1d_signal_stat->GetBinError(ibin);
        sys_err[counter] = h1d_signal_sys->GetBinError(ibin);

        counter++;


    }

    TGraphErrors* graph_stat = new TGraphErrors(counter, x, y, 0, stat_err);
    graph_stat->SetMarkerStyle(20);
    graph_stat->SetMarkerSize(0.9);
    graph_stat->SetMarkerColor(kViolet);
    graph_stat->SetLineColor(kViolet);

    TGraphErrors* graph_sys = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    graph_sys->SetFillColorAlpha(kViolet, 0.3);
    graph_sys->SetLineColor(kViolet);

    pad1->cd();
    graph_stat->Draw("AP");
    graph_sys->Draw("e2same");

    TFile* fdrellyan = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_drell_yan_10k_gencut.root","READ");
    TH1D* h1d_drell_yan_temp = (TH1D*)fdrellyan->Get("h1d_ee_FG_unlike");
    h1d_drell_yan_temp->SetName("h1d_drell_yan_temp");
    h1d_drell_yan_temp->Scale(0.000246/10000.); // Multiplied by Cross section
    h1d_drell_yan_temp->Scale(4.8572);

    TH1D* h1d_drell_yan = new TH1D("h1d_drell_yan","Mass Distribution", nbins_mass_hf, mass_range_hf);
    h1d_drell_yan->Sumw2();

    for(int ibin=0; ibin < nbins_mass_hf; ibin++){

        int bin_lo = h1d_drell_yan_temp->FindBin(mass_range_hf[ibin] + 0.001);
        int bin_hi = h1d_drell_yan_temp->FindBin(mass_range_hf[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_drell_yan_temp->GetBinContent(i)*h1d_drell_yan_temp->GetBinWidth(i);
            bin_error += pow(h1d_drell_yan_temp->GetBinError(i)*h1d_drell_yan_temp->GetBinWidth(i), 2);

        }

        h1d_drell_yan->SetBinContent(ibin+1, bin_content/h1d_drell_yan->GetBinWidth(ibin+1));
        h1d_drell_yan->SetBinError(ibin+1, sqrt(bin_error)/h1d_drell_yan->GetBinWidth(ibin+1));


    }
    

    TH1D* h1d_ccbar_temp = new TH1D("h1d_ccbar_temp", "DCA Distribution", 200, 1.50, 3.50);

    std::ifstream inFile_data("mlp_full_ccbar.txt");
    double x1, y1, err_y1;

    while(!inFile_data.eof()){

        inFile_data >> x1 >> y1 >> err_y1;

        int binX = h1d_ccbar_temp->FindBin(x1);
        h1d_ccbar_temp->SetBinContent(binX, y1);
        h1d_ccbar_temp->SetBinError(binX, err_y1);

    }

    TH1D* h1d_ccbar = new TH1D("h1d_ccbar","DCA Distribution", nbins_mass_hf, mass_range_hf);

    h1d_ccbar->Sumw2();

    for(int ibin=0; ibin < nbins_mass_hf; ibin++){

        double bin_center = (mass_range_hf[ibin] + mass_range_hf[ibin+1])/2.;

        h1d_ccbar->SetBinContent(ibin+1, h1d_ccbar_temp->GetBinContent(h1d_ccbar_temp->FindBin(bin_center)));
        h1d_ccbar->SetBinError(ibin+1, h1d_ccbar_temp->GetBinError(h1d_ccbar_temp->FindBin(bin_center)));

    }

    nbins = h1d_ccbar_temp->GetNbinsX();
    counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = h1d_ccbar_temp->GetBinCenter(ibin);
        y[counter] = h1d_ccbar_temp->GetBinContent(ibin);
        x_lohi[counter] = h1d_ccbar_temp->GetBinWidth(ibin)/2;

        stat_err[counter] = h1d_ccbar_temp->GetBinError(ibin);
        sys_err[counter] = stat_err[counter];

        counter++;


    }

    TGraphErrors* band_ccbar_stat = new TGraphErrors(counter, x, y, 0, stat_err);
    band_ccbar_stat->SetFillColorAlpha(kBlack, 0.3);
    band_ccbar_stat->SetLineColor(kBlack);
    band_ccbar_stat->SetMarkerSize(0.0);
    band_ccbar_stat->SetLineWidth(0);

    TGraphErrors* band_ccbar_sys = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    band_ccbar_sys->SetFillColorAlpha(kBlack, 0.3);
    band_ccbar_sys->SetFillStyle(3022);
    band_ccbar_sys->SetMarkerSize(0.0);
    band_ccbar_sys->SetLineWidth(0);

    band_ccbar_stat->Draw("e3same");
    band_ccbar_sys->Draw("e3same");


    TFile* file = new TFile("/Users/vassu/Desktop/new_files/data_files/producing_templates.root","READ");
    TH1D* h1d_bbbar_plot = (TH1D*)file->Get("h1d_ee_FG_bbbar");
    h1d_bbbar_plot->Scale(4.8572);
    h1d_bbbar_plot->SetName("h1d_bbbar_plot");

    TH1D* h1d_bbbar_central = (TH1D*)h1d_bbbar_plot->Clone();
    h1d_bbbar_central->SetName("h1d_bbbar_central");
    h1d_bbbar_central->Scale(1.8343);

    TH1D* h1d_bbbar_lower = (TH1D*)h1d_bbbar_plot->Clone();
    h1d_bbbar_lower->SetName("h1d_bbbar_lower");
    h1d_bbbar_lower->Scale(1.8343-0.3316);

    TH1D* h1d_bbbar_upper = (TH1D*)h1d_bbbar_plot->Clone();
    h1d_bbbar_upper->SetName("h1d_bbbar_upper");
    h1d_bbbar_upper->Scale(1.8343+0.3316);

    TH1D* h1d_bbbar_temp = (TH1D*)h1d_bbbar_central->Clone();
    h1d_bbbar_temp->SetName("h1d_bbbar_temp");
    h1d_bbbar_temp->Reset("ICESM");

    for(int ibin = 1; ibin < (int)h1d_bbbar_temp->GetNbinsX()+1; ibin++){
        
        
        h1d_bbbar_temp->SetBinContent(ibin+1, 0.8*h1d_bbbar_central->GetBinContent(ibin));
        h1d_bbbar_temp->SetBinError(ibin+1, h1d_bbbar_central->GetBinContent(ibin) - h1d_bbbar_lower->GetBinContent(ibin));

    }

    TH1D* h1d_bbbar = new TH1D("h1d_bbbar","Mass Distribution", nbins_mass_hf, mass_range_hf);
    h1d_bbbar->Sumw2();

    for(int ibin=0; ibin < nbins_mass_hf; ibin++){

        int bin_lo = h1d_bbbar_temp->FindBin(mass_range_hf[ibin] + 0.001);
        int bin_hi = h1d_bbbar_temp->FindBin(mass_range_hf[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_bbbar_temp->GetBinContent(i)*h1d_bbbar_temp->GetBinWidth(i);
            bin_error += pow(h1d_bbbar_temp->GetBinError(i)*h1d_bbbar_temp->GetBinWidth(i), 2);

        }

        h1d_bbbar->SetBinContent(ibin+1, bin_content/h1d_bbbar->GetBinWidth(ibin+1));
        h1d_bbbar->SetBinError(ibin+1, sqrt(bin_error)/h1d_bbbar->GetBinWidth(ibin+1));

    }


    TGraphAsymmErrors* band_bbbar_stat = new TGraphAsymmErrors((int)h1d_bbbar_central->GetNbinsX());
    band_bbbar_stat->SetName("band_bbbar_stat");

    for (int i = 1; i <= (int)h1d_bbbar_central->GetNbinsX(); ++i) {

        double x = h1d_bbbar_central->GetBinCenter(i);
        double y = h1d_bbbar_central->GetBinContent(i);

        double y_up   = h1d_bbbar_upper->GetBinContent(i);
        double y_down = h1d_bbbar_lower->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band_bbbar_stat->SetPoint(i - 1, x, y);
        band_bbbar_stat->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }

    TH1D* h1d_bbbar_central_sys = (TH1D*)h1d_bbbar_plot->Clone();
    h1d_bbbar_central_sys->SetName("h1d_bbbar_central_sys");
    h1d_bbbar_central_sys->Scale(1.8343);

    TH1D* h1d_bbbar_lower_sys = (TH1D*)h1d_bbbar_plot->Clone();
    h1d_bbbar_lower_sys->SetName("h1d_bbbar_lower_sys");
    h1d_bbbar_lower_sys->Scale(1.8343-0.37572);

    TH1D* h1d_bbbar_upper_sys = (TH1D*)h1d_bbbar_plot->Clone();
    h1d_bbbar_upper_sys->SetName("h1d_bbbar_upper_sys");
    h1d_bbbar_upper_sys->Scale(1.8343+0.37572);

    TGraphAsymmErrors* band_bbbar_sys = new TGraphAsymmErrors((int)h1d_bbbar_central_sys->GetNbinsX());
    band_bbbar_sys->SetName("band_bbbar_sys");

    for (int i = 1; i <= (int)h1d_bbbar_central_sys->GetNbinsX(); ++i) {

        double x = h1d_bbbar_central_sys->GetBinCenter(i);
        double y = h1d_bbbar_central_sys->GetBinContent(i);

        double y_up   = h1d_bbbar_upper_sys->GetBinContent(i);
        double y_down = h1d_bbbar_lower_sys->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band_bbbar_sys->SetPoint(i - 1, x, y);
        band_bbbar_sys->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }

    band_bbbar_stat->SetFillColorAlpha(kBlue-4, 0.3);
    band_bbbar_stat->SetLineColor(kBlue-4);
    band_bbbar_stat->SetMarkerSize(0.0);
    band_bbbar_stat->SetLineWidth(0);

    band_bbbar_sys->SetFillColorAlpha(kBlue-4, 0.3);
    band_bbbar_sys->SetFillStyle(3022);
    band_bbbar_sys->SetMarkerSize(0.0);
    band_bbbar_sys->SetLineWidth(0);

    band_bbbar_stat->Draw("e3same");
    band_bbbar_sys->Draw("e3same");


    TH1D* hist_total = (TH1D*)h1d_ccbar->Clone();
    hist_total->SetName("hist_total");
    hist_total->Add(h1d_bbbar);

    hist_total->SetLineColor(kBlue);
    hist_total->SetLineStyle(3);
    hist_total->Draw("psame");

    nbins = hist_total->GetNbinsX();
    counter = 0;


    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = hist_total->GetBinCenter(ibin);
        y[counter] = h1d_signal_stat->GetBinContent(ibin)/hist_total->GetBinContent(ibin);

        x_lohi[counter] = hist_total->GetBinWidth(ibin)/2;

        stat_err[counter] = y[counter]*sqrt( pow(h1d_signal_stat->GetBinError(ibin)/h1d_signal_stat->GetBinContent(ibin), 2) + pow(hist_total->GetBinError(ibin)/hist_total->GetBinContent(ibin), 2) );
  
        sys_err[counter] = y[counter]*sqrt( pow(h1d_signal_sys->GetBinError(ibin)/h1d_signal_sys->GetBinContent(ibin), 2) );

        counter++;


    }
    

    TGraphErrors* graph_stat_ratio = new TGraphErrors(counter, x, y, 0, stat_err);
    graph_stat_ratio->SetMarkerStyle(20);
    graph_stat_ratio->SetMarkerSize(0.9);
    graph_stat_ratio->SetMarkerColor(kRed);
    graph_stat_ratio->SetLineColor(kRed);

    TGraphErrors* graph_sys_ratio = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    graph_sys_ratio->SetFillColorAlpha(kRed, 0.3);
    graph_sys_ratio->SetLineColor(kRed);

    pad2->cd();
    graph_stat_ratio->Draw("AP");
    graph_sys_ratio->Draw("e2same");




}