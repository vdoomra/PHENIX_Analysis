const double nevents = 4.216682807*pow(10,11)*(0.79);
double mass_lo = 0.87;
double mass_hi = 4.50;

const int npts = 10;
double x_lohi[npts] = {-999};
double x[npts] = {-999};
double y[npts] = {-999};
double stat_err[npts] = {-999};
double sys_err[npts] = {-999};

double combinedFunction(double *x, double *par) {

    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    double fn0_val = par[3] + par[4] * x[0] + par[5]*x[0]*x[0];
    double feoverp_val = par[6] + par[7]*x[0] + par[8]*x[0]*x[0];

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}

void azimuthal_correlations(){

    TF1* fcombined = new TF1("fcombined", combinedFunction, 0.0, 5.5, 9);
    fcombined->SetParameters(0.2272, -0.0867, 0.0215, 0.3217, -0.280, 0.0839, 0.3345, -0.1982, 0.05680);

    TFile* fin1 = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_FG_phi_redef_n0_2_with_correction.root","READ");
    TH1D* h1d_FG = (TH1D*)fin1->Get("h1d_ee_FG_unlike");
    h1d_FG->SetName("h1d_FG");
    h1d_FG->Scale(1./nevents);
    h1d_FG->Scale(1.0726*0.545*42.);

    TFile* fin2 = new TFile("plotting_BG_with_uncertainties.root","READ");
    TH1D* h1d_mlp_dataBG = (TH1D*)fin2->Get("h1d_mlp_dataBG");
    h1d_mlp_dataBG->SetName("h1d_mlp_dataBG");

    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
    h1d_FG->Draw();
    h1d_mlp_dataBG->Draw("same hist");

    double h1d_FG_error, h1d_mlp_dataBG_error;
    double S = h1d_FG->IntegralAndError(h1d_FG->FindBin(mass_lo+0.01), h1d_FG->FindBin(mass_hi-0.01), h1d_FG_error, "width");
    double B = h1d_mlp_dataBG->IntegralAndError(h1d_mlp_dataBG->FindBin(mass_lo+0.01), h1d_mlp_dataBG->FindBin(mass_hi-0.01), h1d_mlp_dataBG_error, "width");

    double bkg_fraction = B/S;
    double scale_error = bkg_fraction*sqrt(pow(h1d_FG_error/S, 2) + pow(h1d_mlp_dataBG_error/B, 2));

    TFile* f1 = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_FG_phi_redef_n0_2_with_correction.root","READ");
    TH2D* h2d_delta_phi_mass_FG = (TH2D*)f1->Get("h2d_delta_phi_pair_mass_unlike");
    h2d_delta_phi_mass_FG->SetName("h2d_delta_phi_mass_FG");

    h2d_delta_phi_mass_FG->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_FG = (TH1D*)h2d_delta_phi_mass_FG->ProjectionY();
    h1d_delta_phi_FG->SetName("h1d_delta_phi_FG");

    h1d_delta_phi_FG->Scale(1./nevents);
    h1d_delta_phi_FG->Scale(1.0726*0.545*42.);

    TFile* f2 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_jpsi_2k_gencut_with_correction.root","READ");
    TH2D* h2d_delta_phi_mass_jpsi = (TH2D*)f2->Get("h2d_delta_phi_pair_mass_unlike");
    h2d_delta_phi_mass_jpsi->SetName("h2d_delta_phi_mass_jpsi");

    h2d_delta_phi_mass_jpsi->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_jpsi = (TH1D*)h2d_delta_phi_mass_jpsi->ProjectionY();
    h1d_delta_phi_jpsi->SetName("h1d_delta_phi_jpsi");
    h1d_delta_phi_jpsi->Scale((0.05971*42.)/(25000.*2000.));

    TFile* f3 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_psip_2k_gencut_with_correction.root","READ");
    TH2D* h2d_delta_phi_mass_psip = (TH2D*)f3->Get("h2d_delta_phi_pair_mass_unlike");
    h2d_delta_phi_mass_psip->SetName("h2d_delta_phi_mass_psip");

    h2d_delta_phi_mass_psip->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_psip = (TH1D*)h2d_delta_phi_mass_psip->ProjectionY();
    h1d_delta_phi_psip->SetName("h1d_delta_phi_psip");
    h1d_delta_phi_psip->Scale((7.93e-3*42.)/(25000.*2000.));

    TFile* f4 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_phi_ee_2k_gencut_with_correction.root","READ");
    TH2D* h2d_delta_phi_mass_phi = (TH2D*)f4->Get("h2d_delta_phi_pair_mass_unlike");
    h2d_delta_phi_mass_phi->SetName("h2d_delta_phi_mass_phi");

    h2d_delta_phi_mass_phi->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_phi = (TH1D*)h2d_delta_phi_mass_phi->ProjectionY();
    h1d_delta_phi_phi->SetName("h1d_delta_phi_phi");
    h1d_delta_phi_phi->Scale((0.000297*42.)/(25000.*2000.));

    TFile* f5 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_jet_pairs_reweighted_with_correction.root","READ");
    TH2D* h2d_delta_phi_mass_dataBG = (TH2D*)f5->Get("h2d_delta_phi_pair_mass_unlike");
    h2d_delta_phi_mass_dataBG->SetName("h2d_delta_phi_mass_dataBG");

    h2d_delta_phi_mass_dataBG->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_dataBG = (TH1D*)h2d_delta_phi_mass_dataBG->ProjectionY();
    h1d_delta_phi_dataBG->SetName("h1d_delta_phi_dataBG");

    h1d_delta_phi_dataBG->Scale(h1d_delta_phi_FG->Integral()/h1d_delta_phi_dataBG->Integral());
    h1d_delta_phi_dataBG->Scale(bkg_fraction);

    TH1D* hsum = (TH1D*)h1d_delta_phi_dataBG->Clone();
    hsum->SetName("hsum");
    hsum->Add(h1d_delta_phi_jpsi);
    hsum->Add(h1d_delta_phi_psip);
    hsum->Add(h1d_delta_phi_phi);

    TCanvas* c2 = new TCanvas();
    h1d_delta_phi_FG->Draw();
    h1d_delta_phi_jpsi->Draw("same");
    h1d_delta_phi_psip->Draw("same");
    h1d_delta_phi_phi->Draw("same");
    h1d_delta_phi_dataBG->Draw("same");
    hsum->Draw("same");

    TH1D* h1d_signal = (TH1D*)h1d_delta_phi_FG->Clone();
    h1d_signal->SetName("h1d_signal");
    h1d_signal->Reset("ICESM");
    int nbins = h1d_signal->GetNbinsX();

    for(int ibin = 1; ibin < nbins + 1; ibin++){

        h1d_signal->SetBinContent(ibin, h1d_delta_phi_FG->GetBinContent(ibin) - hsum->GetBinContent(ibin));
        h1d_signal->SetBinError(ibin, sqrt( pow(h1d_delta_phi_FG->GetBinError(ibin), 2) ) );

    }

    int counter = 0;

    for(int ibin=2; ibin < nbins+1; ibin++){

        x[counter] = h1d_signal->GetBinCenter(ibin);
        y[counter] = h1d_signal->GetBinContent(ibin);
        x_lohi[counter] = 0.05;

        stat_err[counter] = h1d_signal->GetBinError(ibin);
        sys_err[counter] = y[counter]*fcombined->Eval(x[counter]);

        counter++;


    }

    TGraphErrors* graph_stat = new TGraphErrors(counter, x, y, 0, stat_err);
    graph_stat->SetMarkerStyle(20);
    graph_stat->SetMarkerSize(0.9);
    graph_stat->SetMarkerColor(kBlack);
    graph_stat->SetLineColor(kBlack);
    graph_stat->SetName("graph_stat");

    TGraphErrors* graph_sys = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    graph_sys->SetFillColorAlpha(kBlack, 0.3);
    graph_sys->SetLineColor(kBlack);
    graph_sys->SetName("graph_sys");

    TMultiGraph* mg = new TMultiGraph();
    mg->SetName("mg");
    
    mg->Add(graph_stat,"AP");
    mg->Add(graph_sys, "e2same");

    TCanvas* c3 = new TCanvas();
    mg->Draw("a");

    mg->GetXaxis()->SetTitle("#Delta#phi[rad]");
    mg->GetYaxis()->SetTitle("#frac{1}{N_{events}}#frac{#epsilon_{BBC}#sigma_{pp}}{#epsilon_{bias}}#frac{dN}{d#Delta#phi} [rad]^{-1} [In PHENIX Acceptance]");

    TLatex *   tex1 = new TLatex(2.0,5e-9,"Run 2015 p+p, #sqrt{s} = 200 GeV, |#eta| < 0.35");
    tex1->SetTextAlign(13);
    tex1->SetTextSize(0.052);
    tex1->SetLineWidth(1);
    tex1->SetTextFont(42);
    tex1->Draw();

    TLatex *   tex4 = new TLatex(3.5,1e-9, "p_{T}^{pair} < 5 GeV/c");
    tex4->SetTextAlign(13);
    tex4->SetTextSize(0.052);
    tex4->SetLineWidth(1);
    tex4->SetTextFont(42);
    tex4->Draw();

    TLatex *   tex5 = new TLatex(3.4,5e-9, "p_{T}^{trig, e} > 1.5 GeV/c");
    tex5->SetTextAlign(13);
    tex5->SetTextSize(0.052);
    tex5->SetLineWidth(1);
    tex5->SetTextFont(42);
    tex5->Draw();


    TFile* fstat = new TFile("scaling_hists.root","READ");

    //ccbar files

    TH2D* h2d_delta_phi_pair_mass_ccbar_central_stat = (TH2D*)fstat->Get("h2d_delta_phi_pair_mass_ccbar_central");
    h2d_delta_phi_pair_mass_ccbar_central_stat->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_ccbar_central_stat = (TH1D*)h2d_delta_phi_pair_mass_ccbar_central_stat->ProjectionY();
    h1d_delta_phi_ccbar_central_stat->SetName("h1d_delta_phi_ccbar_central_stat");

    TH2D* h2d_delta_phi_pair_mass_ccbar_lower_stat = (TH2D*)fstat->Get("h2d_delta_phi_pair_mass_ccbar_lower");
    h2d_delta_phi_pair_mass_ccbar_lower_stat->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_ccbar_lower_stat = (TH1D*)h2d_delta_phi_pair_mass_ccbar_lower_stat->ProjectionY();
    h1d_delta_phi_ccbar_lower_stat->SetName("h1d_delta_phi_ccbar_lower_stat");

    TH2D* h2d_delta_phi_pair_mass_ccbar_upper_stat = (TH2D*)fstat->Get("h2d_delta_phi_pair_mass_ccbar_upper");
    h2d_delta_phi_pair_mass_ccbar_upper_stat->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_ccbar_upper_stat = (TH1D*)h2d_delta_phi_pair_mass_ccbar_upper_stat->ProjectionY();
    h1d_delta_phi_ccbar_upper_stat->SetName("h1d_delta_phi_ccbar_upper_stat");

    // bbbar files

    TH2D* h2d_delta_phi_pair_mass_bbbar_central_stat = (TH2D*)fstat->Get("h2d_delta_phi_pair_mass_bbbar_central");
    h2d_delta_phi_pair_mass_bbbar_central_stat->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_bbbar_central_stat = (TH1D*)h2d_delta_phi_pair_mass_bbbar_central_stat->ProjectionY();
    h1d_delta_phi_bbbar_central_stat->SetName("h1d_delta_phi_bbbar_central_stat");

    TH2D* h2d_delta_phi_pair_mass_bbbar_lower_stat = (TH2D*)fstat->Get("h2d_delta_phi_pair_mass_bbbar_lower");
    h2d_delta_phi_pair_mass_bbbar_lower_stat->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_bbbar_lower_stat = (TH1D*)h2d_delta_phi_pair_mass_bbbar_lower_stat->ProjectionY();
    h1d_delta_phi_bbbar_lower_stat->SetName("h1d_delta_phi_bbbar_lower_stat");

    TH2D* h2d_delta_phi_pair_mass_bbbar_upper_stat = (TH2D*)fstat->Get("h2d_delta_phi_pair_mass_bbbar_upper");
    h2d_delta_phi_pair_mass_bbbar_upper_stat->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_bbbar_upper_stat = (TH1D*)h2d_delta_phi_pair_mass_bbbar_upper_stat->ProjectionY();
    h1d_delta_phi_bbbar_upper_stat->SetName("h1d_delta_phi_bbbar_upper_stat");

    // ccbar + bbbar

    TH1D* h1d_delta_phi_total_central_stat = (TH1D*)h1d_delta_phi_ccbar_central_stat->Clone();
    h1d_delta_phi_total_central_stat->SetName("h1d_delta_phi_total_central_stat");
    h1d_delta_phi_total_central_stat->Add(h1d_delta_phi_bbbar_central_stat);

    TH1D* h1d_delta_phi_total_lower_stat = (TH1D*)h1d_delta_phi_ccbar_lower_stat->Clone();
    h1d_delta_phi_total_lower_stat->SetName("h1d_delta_phi_total_lower_stat");
    h1d_delta_phi_total_lower_stat->Add(h1d_delta_phi_bbbar_lower_stat);

    TH1D* h1d_delta_phi_total_upper_stat = (TH1D*)h1d_delta_phi_ccbar_upper_stat->Clone();
    h1d_delta_phi_total_upper_stat->SetName("h1d_delta_phi_total_upper_stat");
    h1d_delta_phi_total_upper_stat->Add(h1d_delta_phi_bbbar_upper_stat);


    TGraphAsymmErrors* band_stat = new TGraphAsymmErrors((int)h1d_delta_phi_total_central_stat->GetNbinsX());
    band_stat->SetName("band_stat");

    for (int i = 1; i <= (int)h1d_delta_phi_total_central_stat->GetNbinsX(); ++i) {

        double x = h1d_delta_phi_total_central_stat->GetBinCenter(i);
        double y = h1d_delta_phi_total_central_stat->GetBinContent(i);

        double y_up   = h1d_delta_phi_total_upper_stat->GetBinContent(i);
        double y_down = h1d_delta_phi_total_lower_stat->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band_stat->SetPoint(i - 1, x, y);
        band_stat->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }

    TGraphAsymmErrors* band_stat_ccbar = new TGraphAsymmErrors((int)h1d_delta_phi_ccbar_central_stat->GetNbinsX());
    band_stat_ccbar->SetName("band_stat_ccbar");

    for (int i = 1; i <= (int)h1d_delta_phi_ccbar_central_stat->GetNbinsX(); ++i) {

        double x = h1d_delta_phi_ccbar_central_stat->GetBinCenter(i);
        double y = h1d_delta_phi_ccbar_central_stat->GetBinContent(i);

        double y_up   = h1d_delta_phi_ccbar_upper_stat->GetBinContent(i);
        double y_down = h1d_delta_phi_ccbar_lower_stat->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band_stat_ccbar->SetPoint(i - 1, x, y);
        band_stat_ccbar->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }

    TGraphAsymmErrors* band_stat_bbbar = new TGraphAsymmErrors((int)h1d_delta_phi_bbbar_central_stat->GetNbinsX());
    band_stat_bbbar->SetName("band_stat_bbbar");

    for (int i = 1; i <= (int)h1d_delta_phi_bbbar_central_stat->GetNbinsX(); ++i) {

        double x = h1d_delta_phi_bbbar_central_stat->GetBinCenter(i);
        double y = h1d_delta_phi_bbbar_central_stat->GetBinContent(i);

        double y_up   = h1d_delta_phi_bbbar_upper_stat->GetBinContent(i);
        double y_down = h1d_delta_phi_bbbar_lower_stat->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band_stat_bbbar->SetPoint(i - 1, x, y);
        band_stat_bbbar->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }

   TFile* fsys = new TFile("scaling_hists.root","READ");

    //ccbar files

    TH2D* h2d_delta_phi_pair_mass_ccbar_central_sys = (TH2D*)fsys->Get("h2d_delta_phi_pair_mass_ccbar_central");
    h2d_delta_phi_pair_mass_ccbar_central_sys->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_ccbar_central_sys = (TH1D*)h2d_delta_phi_pair_mass_ccbar_central_sys->ProjectionY();
    h1d_delta_phi_ccbar_central_sys->SetName("h1d_delta_phi_ccbar_central_sys");

    TH1D* h1d_delta_phi_ccbar_lower_sys = (TH1D*)h1d_delta_phi_ccbar_central_sys->Clone();
    h1d_delta_phi_ccbar_lower_sys->SetName("h1d_delta_phi_ccbar_lower_sys");
    h1d_delta_phi_ccbar_lower_sys->Scale(1.0-0.3337);

    TH1D* h1d_delta_phi_ccbar_upper_sys = (TH1D*)h1d_delta_phi_ccbar_central_sys->Clone();
    h1d_delta_phi_ccbar_upper_sys->SetName("h1d_delta_phi_ccbar_upper_sys");
    h1d_delta_phi_ccbar_upper_sys->Scale(1.0+0.1756);

    // bbbar files

    TH2D* h2d_delta_phi_pair_mass_bbbar_central_sys = (TH2D*)fsys->Get("h2d_delta_phi_pair_mass_bbbar_central");
    h2d_delta_phi_pair_mass_bbbar_central_sys->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_delta_phi_bbbar_central_sys = (TH1D*)h2d_delta_phi_pair_mass_bbbar_central_sys->ProjectionY();
    h1d_delta_phi_bbbar_central_sys->SetName("h1d_delta_phi_bbbar_central_sys");

    TH1D* h1d_delta_phi_bbbar_lower_sys = (TH1D*)h1d_delta_phi_bbbar_central_sys->Clone();
    h1d_delta_phi_bbbar_lower_sys->SetName("h1d_delta_phi_bbbar_lower_sys");
    h1d_delta_phi_bbbar_lower_sys->Scale(1.0-0.1272);

    TH1D* h1d_delta_phi_bbbar_upper_sys = (TH1D*)h1d_delta_phi_bbbar_central_sys->Clone();
    h1d_delta_phi_bbbar_upper_sys->SetName("h1d_delta_phi_bbbar_upper_sys");
    h1d_delta_phi_bbbar_upper_sys->Scale(1.0+0.2406);

    TGraphAsymmErrors* band_sys_ccbar = new TGraphAsymmErrors((int)h1d_delta_phi_ccbar_central_sys->GetNbinsX());
    band_sys_ccbar->SetName("band_sys_ccbar");

    for (int i = 1; i <= (int)h1d_delta_phi_ccbar_central_sys->GetNbinsX(); ++i) {

        double x = h1d_delta_phi_ccbar_central_sys->GetBinCenter(i);
        double y = h1d_delta_phi_ccbar_central_sys->GetBinContent(i);

        double y_up   = h1d_delta_phi_ccbar_upper_sys->GetBinContent(i);
        double y_down = h1d_delta_phi_ccbar_lower_sys->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band_sys_ccbar->SetPoint(i - 1, x, y);
        band_sys_ccbar->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }

    TGraphAsymmErrors* band_sys_bbbar = new TGraphAsymmErrors((int)h1d_delta_phi_bbbar_central_sys->GetNbinsX());
    band_sys_bbbar->SetName("band_sys_bbbar");

    for (int i = 1; i <= (int)h1d_delta_phi_bbbar_central_sys->GetNbinsX(); ++i) {

        double x = h1d_delta_phi_bbbar_central_sys->GetBinCenter(i);
        double y = h1d_delta_phi_bbbar_central_sys->GetBinContent(i);

        double y_up   = h1d_delta_phi_bbbar_upper_sys->GetBinContent(i);
        double y_down = h1d_delta_phi_bbbar_lower_sys->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band_sys_bbbar->SetPoint(i - 1, x, y);
        band_sys_bbbar->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }


    band_stat->SetFillColorAlpha(kRed, 0.3);
    band_stat->SetLineColor(kRed);
    band_stat->Draw("3same");

    band_stat_ccbar->SetFillColorAlpha(kBlack, 0.3);
    band_stat_ccbar->SetLineColor(kBlack);
    band_stat_ccbar->Draw("3same");

    band_sys_ccbar->SetFillColorAlpha(kBlack, 0.3);
    band_sys_ccbar->SetFillStyle(3002);
    band_sys_ccbar->SetLineColor(kBlack);
    band_sys_ccbar->Draw("3same");

    band_stat_bbbar->SetFillColorAlpha(kBlue-4, 0.3);
    band_stat_bbbar->SetLineColor(kBlue-4);
    band_stat_bbbar->Draw("3same");

    band_sys_bbbar->SetFillColorAlpha(kBlue-4, 0.3);
    band_sys_bbbar->SetFillStyle(3002);
    band_sys_bbbar->SetLineColor(kBlue-4);
    band_sys_bbbar->Draw("3same");

    TLegend* leg = new TLegend(0.2923706,0.6873418,0.7717984,0.8873418);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05164557);
    leg->SetNColumns(2);
    leg->AddEntry(graph_stat, "Data");
    leg->AddEntry(band_stat, "c#bar{c} + b#bar{b} #rightarrow e^{+}e^{-} (PYTHIA8) Stat.");
    leg->AddEntry(band_stat_ccbar, "c#bar{c} #rightarrow e^{+}e^{-} (PYTHIA8) Stat.");
    leg->AddEntry(band_sys_ccbar, " c#bar{c} #rightarrow e^{+}e^{-} (PYTHIA8) Sys.");
    leg->AddEntry(band_stat_bbbar, "b#bar{b} #rightarrow e^{+}e^{-} (PYTHIA8) Stat.");
    leg->AddEntry(band_sys_bbbar, " b#bar{b} #rightarrow e^{+}e^{-} (PYTHIA8) Sys.");
    leg->Draw("same");

    
    
}