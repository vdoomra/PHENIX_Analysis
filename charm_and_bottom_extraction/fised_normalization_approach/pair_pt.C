
const int nbins = 6;
double pt_bins[nbins+1] = {0.0, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0};

const double nevents = 4.216682807*pow(10,11)*(0.79);

double mass_lo = 0.87;
double mass_hi = 4.50;

const int npts = 7;
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

void pair_pt(){

    TF1* fcombined = new TF1("fcombined", combinedFunction, 0.0, 5.5, 9);
    fcombined->SetParameters(0.05814, 0.2883, -0.06865, 0.0229, 0.26705, -0.06064, 0.000573, 0.1994, -0.04178);

    TFile* f1 = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_FG_phi_redef_n0_2_with_correction.root","READ");
    TH2D* h2d_FG = (TH2D*)f1->Get("h2d_ee_FG_unlike");
    h2d_FG->SetName("h2d_FG");

    h2d_FG->Scale(1./nevents);
    h2d_FG->Scale(1.0726*0.545*42);

    h2d_FG->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_pair_pt_FG_temp = (TH1D*)h2d_FG->ProjectionY();
    h1d_pair_pt_FG_temp->SetName("h1d_pair_pt_FG_temp");
    h1d_pair_pt_FG_temp->Scale(10.);

    TH1D* h1d_pair_pt_FG = new TH1D("h1d_pair_pt_FG", "Pair pT Distribution", nbins, pt_bins);
    h1d_pair_pt_FG->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_pair_pt_FG_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_pair_pt_FG_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_pair_pt_FG_temp->GetBinContent(i)*h1d_pair_pt_FG_temp->GetBinWidth(i);
            bin_error += pow(h1d_pair_pt_FG_temp->GetBinError(i)*h1d_pair_pt_FG_temp->GetBinWidth(i), 2);

        }

        h1d_pair_pt_FG->SetBinContent(ibin+1, bin_content/h1d_pair_pt_FG->GetBinWidth(ibin+1));
        h1d_pair_pt_FG->SetBinError(ibin+1, sqrt(bin_error)/h1d_pair_pt_FG->GetBinWidth(ibin+1));


    }

    TFile* f2 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_jpsi_2k_gencut_with_correction.root","READ");
    TH2D* h2d_jpsi = (TH2D*)f2->Get("h2d_ee_FG_unlike");
    h2d_jpsi->SetName("h2d_jpsi");

    h2d_jpsi->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_pair_pt_jpsi_temp = (TH1D*)h2d_jpsi->ProjectionY();
    h1d_pair_pt_jpsi_temp->SetName("h1d_pair_pt_jpsi_temp");
    h1d_pair_pt_jpsi_temp->Scale(10.);

    TH1D* h1d_pair_pt_jpsi = new TH1D("h1d_pair_pt_jpsi", "Pair pT Distribution", nbins, pt_bins);
    h1d_pair_pt_jpsi->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_pair_pt_jpsi_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_pair_pt_jpsi_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_pair_pt_jpsi_temp->GetBinContent(i)*h1d_pair_pt_jpsi_temp->GetBinWidth(i);
            bin_error += pow(h1d_pair_pt_jpsi_temp->GetBinError(i)*h1d_pair_pt_jpsi_temp->GetBinWidth(i), 2);

        }

        h1d_pair_pt_jpsi->SetBinContent(ibin+1, bin_content/h1d_pair_pt_jpsi->GetBinWidth(ibin+1));
        h1d_pair_pt_jpsi->SetBinError(ibin+1, 0.0);


    }

    h1d_pair_pt_jpsi->Scale((42.*0.05971)/(2000.*25000.));

    TFile* f3 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_psip_2k_gencut_with_correction.root","READ");
    TH2D* h2d_psip = (TH2D*)f3->Get("h2d_ee_FG_unlike");
    h2d_psip->SetName("h2d_psip");

    h2d_psip->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_pair_pt_psip_temp = (TH1D*)h2d_psip->ProjectionY();
    h1d_pair_pt_psip_temp->SetName("h1d_pair_pt_psip_temp");
    h1d_pair_pt_psip_temp->Scale(10.);

    TH1D* h1d_pair_pt_psip = new TH1D("h1d_pair_pt_psip", "Pair pT Distribution", nbins, pt_bins);
    h1d_pair_pt_psip->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_pair_pt_psip_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_pair_pt_psip_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_pair_pt_psip_temp->GetBinContent(i)*h1d_pair_pt_psip_temp->GetBinWidth(i);
            bin_error += pow(h1d_pair_pt_psip_temp->GetBinError(i)*h1d_pair_pt_psip_temp->GetBinWidth(i), 2);

        }

        h1d_pair_pt_psip->SetBinContent(ibin+1, bin_content/h1d_pair_pt_psip->GetBinWidth(ibin+1));
        h1d_pair_pt_psip->SetBinError(ibin+1, 0.0);


    }

    h1d_pair_pt_psip->Scale((42.*7.93e-3)/(2000.*25000.));

    TFile* f4 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_phi_ee_2k_gencut_with_correction.root","READ");
    TH2D* h2d_phi = (TH2D*)f4->Get("h2d_ee_FG_unlike");
    h2d_phi->SetName("h2d_phi");

    h2d_phi->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_pair_pt_phi_temp = (TH1D*)h2d_phi->ProjectionY();
    h1d_pair_pt_phi_temp->SetName("h1d_pair_pt_phi_temp");
    h1d_pair_pt_phi_temp->Scale(10.);

    TH1D* h1d_pair_pt_phi = new TH1D("h1d_pair_pt_phi", "Pair pT Distribution", nbins, pt_bins);
    h1d_pair_pt_phi->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_pair_pt_phi_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_pair_pt_phi_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_pair_pt_phi_temp->GetBinContent(i)*h1d_pair_pt_phi_temp->GetBinWidth(i);
            bin_error += pow(h1d_pair_pt_phi_temp->GetBinError(i)*h1d_pair_pt_phi_temp->GetBinWidth(i), 2);

        }

        h1d_pair_pt_phi->SetBinContent(ibin+1, bin_content/h1d_pair_pt_phi->GetBinWidth(ibin+1));
        h1d_pair_pt_phi->SetBinError(ibin+1, sqrt(bin_error)/h1d_pair_pt_phi->GetBinWidth(ibin+1));


    }

    h1d_pair_pt_phi->Scale((42.*0.000297)/(2000.*25000.));

    TFile* fin1 = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_FG_phi_redef_n0_2_with_correction.root","READ");
    TH1D* h1d_FG = (TH1D*)fin1->Get("h1d_ee_FG_unlike");
    h1d_FG->SetName("h1d_FG");
    h1d_FG->Scale(1./nevents);
    h1d_FG->Scale(1.0726*0.545*42);

    TFile* fin2 = new TFile("plotting_BG_with_uncertainties.root","READ");
    TH1D* h1d_mlp_dataBG = (TH1D*)fin2->Get("h1d_mlp_dataBG");
    h1d_mlp_dataBG->SetName("h1d_mlp_dataBG");

    double bkg_fraction = h1d_mlp_dataBG->Integral(h1d_mlp_dataBG->FindBin(mass_lo+0.01), h1d_mlp_dataBG->FindBin(mass_hi-0.01), "width")/h1d_FG->Integral(h1d_FG->FindBin(mass_lo+0.01), h1d_FG->FindBin(mass_hi-0.01), "width");

    TH1D* h1d_pair_pt_BG_temp = new TH1D("h1d_pair_pt_BG_temp", "Pair pT Distribution", 100, 0.0, 5.0);
    double x1, y1, err_y1;

    std::ifstream inFile_data("/Users/vassu/Desktop/new_files/data_files/mlp_JetPairs_pair_pt.txt");

    while(!inFile_data.eof()){

        inFile_data >> x1 >> y1 >> err_y1;

        int binX = h1d_pair_pt_BG_temp->FindBin(x1);
        h1d_pair_pt_BG_temp->SetBinContent(binX, y1);
        h1d_pair_pt_BG_temp->SetBinError(binX, err_y1);

    }

    TH1D* h1d_pair_pt_BG = new TH1D("h1d_pair_pt_BG", "Pair pT Distribution", nbins, pt_bins);
    h1d_pair_pt_BG->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_pair_pt_BG_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_pair_pt_BG_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_pair_pt_BG_temp->GetBinContent(i)*h1d_pair_pt_BG_temp->GetBinWidth(i);
            bin_error += pow(h1d_pair_pt_BG_temp->GetBinError(i)*h1d_pair_pt_BG_temp->GetBinWidth(i), 2);

        }

        h1d_pair_pt_BG->SetBinContent(ibin+1, bin_content/h1d_pair_pt_BG->GetBinWidth(ibin+1));
        h1d_pair_pt_BG->SetBinError(ibin+1, 0.0);


    }

    h1d_pair_pt_BG->Scale(h1d_pair_pt_FG->Integral("width")/h1d_pair_pt_BG->Integral("width"));
    h1d_pair_pt_BG->Scale(bkg_fraction);


    TH1D* hsum = (TH1D*)h1d_pair_pt_jpsi->Clone();
    hsum->SetName("hsum");
    hsum->Add(h1d_pair_pt_BG);
    hsum->Add(h1d_pair_pt_psip);
    hsum->Add(h1d_pair_pt_phi);

    TH1D* h1d_signal = (TH1D*)h1d_pair_pt_FG->Clone();
    h1d_signal->SetName("h1d_signal");
    h1d_signal->Reset("ICESM");

    for(int ibin = 1; ibin < (int)h1d_signal->GetNbinsX()+1; ibin++){

        double bin_content = h1d_pair_pt_FG->GetBinContent(ibin) - hsum->GetBinContent(ibin);
        double bin_error = h1d_pair_pt_FG->GetBinError(ibin);
        h1d_signal->SetBinContent(ibin, bin_content);
        h1d_signal->SetBinError(ibin, bin_error);

    }

    TCanvas* c2 = new TCanvas();
    c2->SetLogy();
    h1d_pair_pt_FG->Draw();
    h1d_pair_pt_BG->Draw("same");
    h1d_pair_pt_jpsi->Draw("same");
    h1d_pair_pt_phi->Draw("same");
    h1d_pair_pt_psip->Draw("same");
    hsum->Draw("same");

    int counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

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

    mg->GetXaxis()->SetTitle("p_{T}^{pair}[GeV/c]");
    mg->GetYaxis()->SetTitle("#frac{1}{N_{events}}#frac{#epsilon_{BBC}#sigma_{pp}}{#epsilon_{bias}}#frac{dN}{dp_{T}^{pair}} [GeV/c]^{-1} [In PHENIX Acceptance]");

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


    TH1D* h1d_pair_pt_ccbar = new TH1D("h1d_pair_pt_ccbar", "Pair pT Distribution", 200, 0.0, 5.0);

    std::ifstream inFile_ccbar("mlp_pair_pt_ccbar.txt");

    while(!inFile_ccbar.eof()){

        inFile_ccbar >> x1 >> y1 >> err_y1;

        int binX = h1d_pair_pt_ccbar->FindBin(x1);
        h1d_pair_pt_ccbar->SetBinContent(binX, y1);
        h1d_pair_pt_ccbar->SetBinError(binX, err_y1);

    }

    TH1D* h1d_pair_pt_bbbar = new TH1D("h1d_pair_pt_bbbar", "Pair pT Distribution", 200, 0.0, 5.0);

    std::ifstream inFile_bbbar("mlp_pair_pt_bbbar.txt");

    while(!inFile_bbbar.eof()){

        inFile_bbbar >> x1 >> y1 >> err_y1;

        int binX = h1d_pair_pt_bbbar->FindBin(x1);
        h1d_pair_pt_bbbar->SetBinContent(binX, y1);
        h1d_pair_pt_bbbar->SetBinError(binX, err_y1);

    }

    h1d_pair_pt_ccbar->Scale(10.);
    h1d_pair_pt_bbbar->Scale(10.);

    TH1D* h1d_pair_pt_drellyan = new TH1D("h1d_pair_pt_drellyan", "Pair pT Distribution", 200, 0.0, 5.0);

    std::ifstream inFile_drellyan("mlp_pair_pt_drellyan.txt");

    while(!inFile_drellyan.eof()){

        inFile_drellyan >> x1 >> y1 >> err_y1;

        int binX = h1d_pair_pt_drellyan->FindBin(x1);
        h1d_pair_pt_drellyan->SetBinContent(binX, y1);
        h1d_pair_pt_drellyan->SetBinError(binX, err_y1);

    }

    h1d_pair_pt_drellyan->Scale(10.);
    h1d_pair_pt_drellyan->Scale(0.000246/10000.); // Multiplied by Cross section
    h1d_pair_pt_drellyan->Scale(4.8572);


    TH1D* h1d_pair_pt_total_central_stat = (TH1D*)h1d_pair_pt_ccbar->Clone();
    h1d_pair_pt_total_central_stat->SetName("h1d_pair_pt_total_central_stat");
    h1d_pair_pt_total_central_stat->Add(h1d_pair_pt_bbbar);
    h1d_pair_pt_total_central_stat->Add(h1d_pair_pt_drellyan);

    TH1D* h1d_pair_pt_total_lower_stat = (TH1D*)h1d_pair_pt_total_central_stat->Clone();
    h1d_pair_pt_total_lower_stat->SetName("h1d_pair_pt_total_lower_stat");
    h1d_pair_pt_total_lower_stat->Scale(1-0.1053);

    TH1D* h1d_pair_pt_total_upper_stat = (TH1D*)h1d_pair_pt_total_central_stat->Clone();
    h1d_pair_pt_total_upper_stat->SetName("h1d_pair_pt_total_upper_stat");
    h1d_pair_pt_total_upper_stat->Scale(1+0.1053);

    TGraphAsymmErrors* band_stat = new TGraphAsymmErrors((int)h1d_pair_pt_total_central_stat->GetNbinsX());
    band_stat->SetName("band_stat");

    for (int i = 1; i <= (int)h1d_pair_pt_total_central_stat->GetNbinsX(); ++i) {

        double x = h1d_pair_pt_total_central_stat->GetBinCenter(i);
        double y = h1d_pair_pt_total_central_stat->GetBinContent(i);

        double y_up   = h1d_pair_pt_total_upper_stat->GetBinContent(i);
        double y_down = h1d_pair_pt_total_lower_stat->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band_stat->SetPoint(i - 1, x, y);
        band_stat->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }

    TH1D* h1d_pair_pt_total_central_sys = (TH1D*)h1d_pair_pt_ccbar->Clone();
    h1d_pair_pt_total_central_sys->SetName("h1d_pair_pt_total_central_sys");
    h1d_pair_pt_total_central_sys->Add(h1d_pair_pt_bbbar);
    h1d_pair_pt_total_central_sys->Add(h1d_pair_pt_drellyan);

    TH1D* h1d_pair_pt_total_lower_sys = (TH1D*)h1d_pair_pt_total_central_sys->Clone();
    h1d_pair_pt_total_lower_sys->SetName("h1d_pair_pt_total_lower_sys");
    h1d_pair_pt_total_lower_sys->Scale(1-0.1749);

    TH1D* h1d_pair_pt_total_upper_sys = (TH1D*)h1d_pair_pt_total_central_sys->Clone();
    h1d_pair_pt_total_upper_sys->SetName("h1d_pair_pt_total_upper_sys");
    h1d_pair_pt_total_upper_sys->Scale(1+0.1749);

    TGraphAsymmErrors* band_sys = new TGraphAsymmErrors((int)h1d_pair_pt_total_central_sys->GetNbinsX());
    band_sys->SetName("band_sys");

    for (int i = 1; i <= (int)h1d_pair_pt_total_central_sys->GetNbinsX(); ++i) {

        double x = h1d_pair_pt_total_central_sys->GetBinCenter(i);
        double y = h1d_pair_pt_total_central_sys->GetBinContent(i);

        double y_up   = h1d_pair_pt_total_upper_sys->GetBinContent(i);
        double y_down = h1d_pair_pt_total_lower_sys->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band_sys->SetPoint(i - 1, x, y);
        band_sys->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }
    
    band_stat->SetFillColorAlpha(kBlack, 0.3);
    band_stat->SetLineColor(kBlack);
    band_stat->Draw("3same");

    band_sys->SetFillColorAlpha(kBlack, 0.3);
    band_sys->SetFillStyle(3002);
    band_sys->SetLineColor(kBlack);
    band_sys->Draw("3same");

    h1d_pair_pt_ccbar->SetLineColor(kRed);
    h1d_pair_pt_ccbar->SetLineStyle(2);
    h1d_pair_pt_ccbar->Draw("same hist");

    h1d_pair_pt_bbbar->SetLineColor(kGreen-7);
    h1d_pair_pt_bbbar->SetLineStyle(2);
    h1d_pair_pt_bbbar->Draw("same hist");

    h1d_pair_pt_drellyan->SetLineColor(kBlue-7);
    h1d_pair_pt_drellyan->SetLineStyle(2);
    h1d_pair_pt_drellyan->Draw("same hist");

    TLegend* leg = new TLegend(0.2923706,0.6873418,0.7717984,0.8873418);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05164557);
    leg->SetNColumns(2);
    leg->AddEntry(graph_stat, "Data");
    leg->AddEntry(band_stat, "c#bar{c} + b#bar{b} + DY #rightarrow e^{+}e^{-} (PYTHIA8) Stat.");
    leg->AddEntry(band_sys, "c#bar{c} + b#bar{b} + DY #rightarrow e^{+}e^{-} (PYTHIA8) Sys.");
    leg->AddEntry(h1d_pair_pt_ccbar, "c#bar{c} #rightarrow e^{+}e^{-} (PYTHIA8)");
    leg->AddEntry(h1d_pair_pt_bbbar, " b#bar{b} #rightarrow e^{+}e^{-} (PYTHIA8)");
    leg->AddEntry(h1d_pair_pt_drellyan, " DY #rightarrow e^{+}e^{-} (PYTHIA8)");
    leg->Draw("same");


    
}