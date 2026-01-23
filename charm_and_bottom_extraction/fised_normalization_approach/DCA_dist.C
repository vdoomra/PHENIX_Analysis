
const int nbins = 18;
double pt_bins[nbins+1] = {0.0, 40, 80, 120, 160, 200, 240, 280, 320, 400, 520, 640, 760, 880, 1000, 1160, 1280, 1400, 1520};

const double nevents = 4.216682807*pow(10,11)*(0.79);

double mass_lo = 0.87;
double mass_hi = 4.50;

const int npts = 18;
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

double combinedFunction_jpsi(double *x, double *par) {

    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    double fn0_val = par[3] + par[4] * x[0] + par[5]*x[0]*x[0];
    double feoverp_val = par[6] + par[7]*x[0] + par[8]*x[0]*x[0];

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}


void DCA_dist(){

     TCanvas *c1 = new TCanvas("c", "Unified Two Panel Plot", 800, 800);

    TPad *pad1 = new TPad("pad1", "Top pad", 0.0, 0.35, 1.0, 1.0);
    TPad *pad2 = new TPad("pad2", "Bottom pad", 0.0, 0.0, 1.0, 0.35);


    pad1->SetBottomMargin(0.0); 
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.35);
    pad1->SetTickx(0);

    pad1->Draw();
    pad2->Draw();

    TF1* fcombined = new TF1("fcombined", combinedFunction, 0, 1520, 9);
    fcombined->SetParameters(0.0732, -7.383e-5, 9.587e-8, 0.07187, -0.000144, 1.462e-7, 0.0455, 0.000142, -6.848e-8);

    TF1* fcombined_jpsi = new TF1("fcombined_jpsi", combinedFunction_jpsi, 0.0, 5.5, 9);
    fcombined_jpsi->SetParameters(0.0513, -6.75e-5, 9.95e-8, 0.0412, -1.09e-5, 8.959e-8, 0.0299, -3.89e-5, 1.55e-7);

    TFile* f1 = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_FG_phi_redef_n0_2_with_correction.root","READ");
    TH2D* h2d_FG = (TH2D*)f1->Get("hist_unlike");
    h2d_FG->SetName("h2d_FG");

    h2d_FG->Scale(1./nevents);
    h2d_FG->Scale(1.0726*0.545*42);

    h2d_FG->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_FG_temp = (TH1D*)h2d_FG->ProjectionY();
    h1d_dca_FG_temp->SetName("h1d_dca_FG_temp");
    h1d_dca_FG_temp->Scale(40.);

    TH1D* h1d_dca_FG = new TH1D("h1d_dca_FG", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_FG->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_FG_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_FG_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_FG_temp->GetBinContent(i)*h1d_dca_FG_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_FG_temp->GetBinError(i)*h1d_dca_FG_temp->GetBinWidth(i), 2);

        }

        h1d_dca_FG->SetBinContent(ibin+1, bin_content/h1d_dca_FG->GetBinWidth(ibin+1));
        h1d_dca_FG->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_FG->GetBinWidth(ibin+1));


    }

    TFile* f2 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_jpsi_2k_gencut_with_correction.root","READ");
    TH2D* h2d_jpsi = (TH2D*)f2->Get("hist_unlike");
    h2d_jpsi->SetName("h2d_jpsi");

    h2d_jpsi->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_jpsi_temp = (TH1D*)h2d_jpsi->ProjectionY();
    h1d_dca_jpsi_temp->SetName("h1d_dca_jpsi_temp");
    h1d_dca_jpsi_temp->Scale(40.);

    TH1D* h1d_dca_jpsi = new TH1D("h1d_dca_jpsi", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_jpsi->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_jpsi_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_jpsi_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_jpsi_temp->GetBinContent(i)*h1d_dca_jpsi_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_jpsi_temp->GetBinError(i)*h1d_dca_jpsi_temp->GetBinWidth(i), 2);

        }

        h1d_dca_jpsi->SetBinContent(ibin+1, bin_content/h1d_dca_jpsi->GetBinWidth(ibin+1));
        h1d_dca_jpsi->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_jpsi->GetBinWidth(ibin+1));


    }

    h1d_dca_jpsi->Scale((42.*0.05971)/(2000.*25000.));

    TFile* f3 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_psip_2k_gencut_with_correction.root","READ");
    TH2D* h2d_psip = (TH2D*)f3->Get("hist_unlike");
    h2d_psip->SetName("h2d_psip");

    h2d_psip->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_psip_temp = (TH1D*)h2d_psip->ProjectionY();
    h1d_dca_psip_temp->SetName("h1d_dca_psip_temp");
    h1d_dca_psip_temp->Scale(40.);

    TH1D* h1d_dca_psip = new TH1D("h1d_dca_psip", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_psip->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_psip_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_psip_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_psip_temp->GetBinContent(i)*h1d_dca_psip_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_psip_temp->GetBinError(i)*h1d_dca_psip_temp->GetBinWidth(i), 2);

        }

        h1d_dca_psip->SetBinContent(ibin+1, bin_content/h1d_dca_psip->GetBinWidth(ibin+1));
        h1d_dca_psip->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_psip->GetBinWidth(ibin+1));


    }

    h1d_dca_psip->Scale((42.*7.93e-3)/(2000.*25000.));

    TFile* f4 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_phi_ee_2k_gencut_with_correction.root","READ");
    TH2D* h2d_phi = (TH2D*)f4->Get("hist_unlike");
    h2d_phi->SetName("h2d_phi");

    h2d_phi->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_phi_temp = (TH1D*)h2d_phi->ProjectionY();
    h1d_dca_phi_temp->SetName("h1d_dca_phi_temp");
    h1d_dca_phi_temp->Scale(40.);

    TH1D* h1d_dca_phi = new TH1D("h1d_dca_phi", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_phi->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_phi_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_phi_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_phi_temp->GetBinContent(i)*h1d_dca_phi_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_phi_temp->GetBinError(i)*h1d_dca_phi_temp->GetBinWidth(i), 2);

        }

        h1d_dca_phi->SetBinContent(ibin+1, bin_content/h1d_dca_phi->GetBinWidth(ibin+1));
        h1d_dca_phi->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_phi->GetBinWidth(ibin+1));


    }

    h1d_dca_phi->Scale((42.*0.000297)/(2000.*25000.));

    TFile* fin1 = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_FG_phi_redef_n0_2_with_correction.root","READ");
    TH1D* h1d_FG = (TH1D*)fin1->Get("h1d_ee_FG_unlike");
    h1d_FG->SetName("h1d_FG");
    h1d_FG->Scale(1./nevents);
    h1d_FG->Scale(1.0726*0.545*42);

    TFile* fin2 = new TFile("plotting_BG_with_uncertainties.root","READ");
    TH1D* h1d_mlp_dataBG = (TH1D*)fin2->Get("h1d_mlp_dataBG");
    h1d_mlp_dataBG->SetName("h1d_mlp_dataBG");

    double bkg_fraction = h1d_mlp_dataBG->Integral(h1d_mlp_dataBG->FindBin(mass_lo+0.01), h1d_mlp_dataBG->FindBin(mass_hi-0.01), "width")/h1d_FG->Integral(h1d_FG->FindBin(mass_lo+0.01), h1d_FG->FindBin(mass_hi-0.01), "width");

    TFile* f5 = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_BG_phi_redef_n0_2_with_correction.root","READ");
    TH2D* h2d_BG = (TH2D*)f5->Get("hist_unlike");
    h2d_BG->SetName("h2d_BG");

    h2d_BG->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_BG_temp = (TH1D*)h2d_BG->ProjectionY();
    h1d_dca_BG_temp->SetName("h1d_dca_BG_temp");
    h1d_dca_BG_temp->Scale(40.);

    TH1D* h1d_dca_BG = new TH1D("h1d_dca_BG", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_BG->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_BG_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_BG_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_BG_temp->GetBinContent(i)*h1d_dca_BG_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_BG_temp->GetBinError(i)*h1d_dca_BG_temp->GetBinWidth(i), 2);

        }

        h1d_dca_BG->SetBinContent(ibin+1, bin_content/h1d_dca_BG->GetBinWidth(ibin+1));
        h1d_dca_BG->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_BG->GetBinWidth(ibin+1));


    }

    h1d_dca_BG->Scale(h1d_dca_FG->Integral("width")/h1d_dca_BG->Integral("width"));
    h1d_dca_BG->Scale(bkg_fraction);

    TFile* f6 = new TFile("heavy_flavor_continuum_stat.root","READ");
    TH2D* h2d_ccbar_central = (TH2D*)f6->Get("hist_dca_ccbar_central");
    h2d_ccbar_central->SetName("h2d_ccbar_central");

    h2d_ccbar_central->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_ccbar_central_temp = (TH1D*)h2d_ccbar_central->ProjectionY();
    h1d_dca_ccbar_central_temp->SetName("h1d_dca_ccbar_central_temp");
    h1d_dca_ccbar_central_temp->Scale(40.);

    TH1D* h1d_dca_ccbar_central = new TH1D("h1d_dca_ccbar_central", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_ccbar_central->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_ccbar_central_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_ccbar_central_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_ccbar_central_temp->GetBinContent(i)*h1d_dca_ccbar_central_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_ccbar_central_temp->GetBinError(i)*h1d_dca_ccbar_central_temp->GetBinWidth(i), 2);

        }

        h1d_dca_ccbar_central->SetBinContent(ibin+1, bin_content/h1d_dca_ccbar_central->GetBinWidth(ibin+1));
        h1d_dca_ccbar_central->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_ccbar_central->GetBinWidth(ibin+1));


    }


    TH2D* h2d_bbbar_central = (TH2D*)f6->Get("hist_dca_bbbar_central");
    h2d_bbbar_central->SetName("h2d_bbbar_central");

    h2d_bbbar_central->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_bbbar_central_temp = (TH1D*)h2d_bbbar_central->ProjectionY();
    h1d_dca_bbbar_central_temp->SetName("h1d_dca_bbbar_central_temp");
    h1d_dca_bbbar_central_temp->Scale(40.);

    TH1D* h1d_dca_bbbar_central = new TH1D("h1d_dca_bbbar_central", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_bbbar_central->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_bbbar_central_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_bbbar_central_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_bbbar_central_temp->GetBinContent(i)*h1d_dca_bbbar_central_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_bbbar_central_temp->GetBinError(i)*h1d_dca_bbbar_central_temp->GetBinWidth(i), 2);

        }

        h1d_dca_bbbar_central->SetBinContent(ibin+1, bin_content/h1d_dca_bbbar_central->GetBinWidth(ibin+1));
        h1d_dca_bbbar_central->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_bbbar_central->GetBinWidth(ibin+1));


    }

    TH2D* h2d_hf_central = (TH2D*)f6->Get("hist_dca_total_central");
    h2d_hf_central->SetName("h2d_hf_central");

    h2d_hf_central->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_hf_central_temp = (TH1D*)h2d_hf_central->ProjectionY();
    h1d_dca_hf_central_temp->SetName("h1d_dca_hf_central_temp");
    h1d_dca_hf_central_temp->Scale(40.);

    TH1D* h1d_dca_hf_central = new TH1D("h1d_dca_hf_central", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_hf_central->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_hf_central_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_hf_central_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_hf_central_temp->GetBinContent(i)*h1d_dca_hf_central_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_hf_central_temp->GetBinError(i)*h1d_dca_hf_central_temp->GetBinWidth(i), 2);

        }

        h1d_dca_hf_central->SetBinContent(ibin+1, bin_content/h1d_dca_hf_central->GetBinWidth(ibin+1));
        h1d_dca_hf_central->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_hf_central->GetBinWidth(ibin+1));


    }

    TH2D* h2d_hf_lower = (TH2D*)f6->Get("hist_dca_total_lower");
    h2d_hf_lower->SetName("h2d_hf_lower");

    h2d_hf_lower->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_hf_lower_temp = (TH1D*)h2d_hf_lower->ProjectionY();
    h1d_dca_hf_lower_temp->SetName("h1d_dca_hf_lower_temp");
    h1d_dca_hf_lower_temp->Scale(40.);

    TH1D* h1d_dca_hf_lower = new TH1D("h1d_dca_hf_lower", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_hf_lower->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_hf_lower_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_hf_lower_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_hf_lower_temp->GetBinContent(i)*h1d_dca_hf_lower_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_hf_lower_temp->GetBinError(i)*h1d_dca_hf_lower_temp->GetBinWidth(i), 2);

        }

        h1d_dca_hf_lower->SetBinContent(ibin+1, bin_content/h1d_dca_hf_lower->GetBinWidth(ibin+1));
        h1d_dca_hf_lower->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_hf_lower->GetBinWidth(ibin+1));


    }

    TH2D* h2d_hf_upper = (TH2D*)f6->Get("hist_dca_total_upper");
    h2d_hf_upper->SetName("h2d_hf_upper");

    h2d_hf_upper->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_hf_upper_temp = (TH1D*)h2d_hf_upper->ProjectionY();
    h1d_dca_hf_upper_temp->SetName("h1d_dca_hf_upper_temp");
    h1d_dca_hf_upper_temp->Scale(40.);

    TH1D* h1d_dca_hf_upper = new TH1D("h1d_dca_hf_upper", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_hf_upper->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_hf_upper_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_hf_upper_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_hf_upper_temp->GetBinContent(i)*h1d_dca_hf_upper_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_hf_upper_temp->GetBinError(i)*h1d_dca_hf_upper_temp->GetBinWidth(i), 2);

        }

        h1d_dca_hf_upper->SetBinContent(ibin+1, bin_content/h1d_dca_hf_upper->GetBinWidth(ibin+1));
        h1d_dca_hf_upper->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_hf_upper->GetBinWidth(ibin+1));


    }

    TH1D* h1d_dca_hf = new TH1D("h1d_dca_hf","Mass Distribution", nbins, pt_bins);
    h1d_dca_hf->Sumw2();

    for(int ibin = 0; ibin < nbins; ibin++){
        
        double bin_center = (pt_bins[ibin] + pt_bins[ibin+1])/2.;
        h1d_dca_hf->SetBinContent(ibin+1, h1d_dca_hf_central->GetBinContent(h1d_dca_hf_central->FindBin(bin_center)));
        h1d_dca_hf->SetBinError(ibin+1, h1d_dca_hf_central->GetBinContent(h1d_dca_hf_central->FindBin(bin_center)) - h1d_dca_hf_lower->GetBinContent(h1d_dca_hf_lower->FindBin(bin_center)));

    }


    TFile* f7 = new TFile("heavy_flavor_continuum_sys.root","READ");

    TH2D* h2d_hf_central_sys = (TH2D*)f7->Get("hist_dca_total_central");
    h2d_hf_central_sys->SetName("h2d_hf_central_sys");

    h2d_hf_central_sys->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_hf_central_sys_temp = (TH1D*)h2d_hf_central_sys->ProjectionY();
    h1d_dca_hf_central_sys_temp->SetName("h1d_dca_hf_central_sys_temp");
    h1d_dca_hf_central_sys_temp->Scale(40.);

    TH1D* h1d_dca_hf_central_sys = new TH1D("h1d_dca_hf_central_sys", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_hf_central_sys->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_hf_central_sys_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_hf_central_sys_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_hf_central_sys_temp->GetBinContent(i)*h1d_dca_hf_central_sys_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_hf_central_sys_temp->GetBinError(i)*h1d_dca_hf_central_sys_temp->GetBinWidth(i), 2);

        }

        h1d_dca_hf_central_sys->SetBinContent(ibin+1, bin_content/h1d_dca_hf_central_sys->GetBinWidth(ibin+1));
        h1d_dca_hf_central_sys->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_hf_central_sys->GetBinWidth(ibin+1));


    }

    TH2D* h2d_hf_lower_sys = (TH2D*)f7->Get("hist_dca_total_lower");
    h2d_hf_lower_sys->SetName("h2d_hf_lower_sys");

    h2d_hf_lower_sys->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_hf_lower_sys_temp = (TH1D*)h2d_hf_lower_sys->ProjectionY();
    h1d_dca_hf_lower_sys_temp->SetName("h1d_dca_hf_lower_sys_temp");
    h1d_dca_hf_lower_sys_temp->Scale(40.);

    TH1D* h1d_dca_hf_lower_sys = new TH1D("h1d_dca_hf_lower_sys", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_hf_lower_sys->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_hf_lower_sys_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_hf_lower_sys_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_hf_lower_sys_temp->GetBinContent(i)*h1d_dca_hf_lower_sys_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_hf_lower_sys_temp->GetBinError(i)*h1d_dca_hf_lower_sys_temp->GetBinWidth(i), 2);

        }

        h1d_dca_hf_lower_sys->SetBinContent(ibin+1, bin_content/h1d_dca_hf_lower_sys->GetBinWidth(ibin+1));
        h1d_dca_hf_lower_sys->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_hf_lower_sys->GetBinWidth(ibin+1));


    }

    TH2D* h2d_hf_upper_sys = (TH2D*)f7->Get("hist_dca_total_upper");
    h2d_hf_upper_sys->SetName("h2d_hf_upper_sys");

    h2d_hf_upper_sys->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_dca_hf_upper_sys_temp = (TH1D*)h2d_hf_upper_sys->ProjectionY();
    h1d_dca_hf_upper_sys_temp->SetName("h1d_dca_hf_upper_sys_temp");
    h1d_dca_hf_upper_sys_temp->Scale(40.);

    TH1D* h1d_dca_hf_upper_sys = new TH1D("h1d_dca_hf_upper_sys", "Pair pT Distribution", nbins, pt_bins);
    h1d_dca_hf_upper_sys->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_dca_hf_upper_sys_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_dca_hf_upper_sys_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dca_hf_upper_sys_temp->GetBinContent(i)*h1d_dca_hf_upper_sys_temp->GetBinWidth(i);
            bin_error += pow(h1d_dca_hf_upper_sys_temp->GetBinError(i)*h1d_dca_hf_upper_sys_temp->GetBinWidth(i), 2);

        }

        h1d_dca_hf_upper_sys->SetBinContent(ibin+1, bin_content/h1d_dca_hf_upper_sys->GetBinWidth(ibin+1));
        h1d_dca_hf_upper_sys->SetBinError(ibin+1, sqrt(bin_error)/h1d_dca_hf_upper_sys->GetBinWidth(ibin+1));


    }

    TH1D* h1d_dca_hf_sys = new TH1D("h1d_dca_hf_sys","Mass Distribution", nbins, pt_bins);
    h1d_dca_hf_sys->Sumw2();

    for(int ibin = 0; ibin < nbins; ibin++){
        
        double bin_center = (pt_bins[ibin] + pt_bins[ibin+1])/2.;
        h1d_dca_hf_sys->SetBinContent(ibin+1, h1d_dca_hf_central_sys->GetBinContent(h1d_dca_hf_central_sys->FindBin(bin_center)));
        h1d_dca_hf_sys->SetBinError(ibin+1, h1d_dca_hf_central_sys->GetBinContent(h1d_dca_hf_central_sys->FindBin(bin_center)) - h1d_dca_hf_lower_sys->GetBinContent(h1d_dca_hf_lower_sys->FindBin(bin_center)));

    }

    

    TH1D* hsum = (TH1D*)h1d_dca_jpsi->Clone();
    hsum->SetName("hsum");
    hsum->Add(h1d_dca_BG);
    hsum->Add(h1d_dca_psip);
    hsum->Add(h1d_dca_phi);
    hsum->Add(h1d_dca_hf);

    int counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = h1d_dca_FG->GetBinCenter(ibin);
        y[counter] = h1d_dca_FG->GetBinContent(ibin);
        x_lohi[counter] = h1d_dca_FG->GetBinWidth(ibin)/4;

        stat_err[counter] = h1d_dca_FG->GetBinError(ibin);
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

    pad1->cd();
    mg->Draw("a");
    mg->GetXaxis()->SetRangeUser(10, 1500);

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

    TLatex *   tex6 = new TLatex(3.4,5e-9, "0.87 < m_{ee} << 4.50 GeV/c^{2}");
    tex6->SetTextAlign(13);
    tex6->SetTextSize(0.052);
    tex6->SetLineWidth(1);
    tex6->SetTextFont(42);
    tex6->Draw();

    counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = hsum->GetBinCenter(ibin);
        y[counter] = hsum->GetBinContent(ibin);
        x_lohi[counter] = hsum->GetBinWidth(ibin)/4;

        stat_err[counter] = hsum->GetBinError(ibin);

        double sys_jpsi = fcombined_jpsi->Eval(x[counter])*h1d_dca_jpsi->GetBinContent(ibin);
        double sys_bg = h1d_mlp_dataBG->GetBinContent(ibin)*0.05212;
        double sys_hf = h1d_dca_hf_sys->GetBinError(ibin);

        sys_err[counter] = sqrt( pow(sys_jpsi, 2) + pow(sys_bg, 2) + pow(sys_hf, 2) );

        counter++;


    }

    TGraphErrors* graph_stat_cocktail = new TGraphErrors(counter, x, y, 0, stat_err);
    graph_stat_cocktail->SetMarkerStyle(20);
    graph_stat_cocktail->SetMarkerSize(0.9);
    graph_stat_cocktail->SetMarkerColor(kViolet-4);
    graph_stat_cocktail->SetLineColor(kViolet-4);

    TGraphErrors* graph_sys_cocktail = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    graph_sys_cocktail->SetFillColorAlpha(kViolet-4, 0.3);
    graph_sys_cocktail->SetLineColor(kViolet-4);

    graph_stat_cocktail->Draw("psame");
    graph_sys_cocktail->Draw("e2same");

    h1d_dca_bbbar_central->Smooth();


    h1d_dca_jpsi->Draw("same hist");
    h1d_dca_phi->Draw("same hist");
    h1d_dca_psip->Draw("same hist");
    h1d_dca_BG->Draw("same hist");
    h1d_dca_ccbar_central->Draw("same hist");
    h1d_dca_bbbar_central->Draw("same hist");

    TLegend* leg = new TLegend(0.2923706,0.6873418,0.7717984,0.8873418);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05164557);
    leg->SetNColumns(2);
    leg->AddEntry(graph_stat, "FG_{e^{+}e^{-}}");
    leg->AddEntry(h1d_dca_BG, "BG_{e^{+}e^{-}}");
    leg->AddEntry(h1d_dca_jpsi, "J/#psi #rightarrow e^{+}e^{-}");
    leg->AddEntry(h1d_dca_phi, "#phi #rightarrow e^{+}e^{-}");
    leg->AddEntry(h1d_dca_psip, "#psi' #rightarrow e^{+}e^{-}");
    leg->AddEntry(h1d_dca_ccbar_central, "c#bar{c} #rightarrow e^{+}e^{-}");
    leg->AddEntry(h1d_dca_bbbar_central, "b#bar{b} #rightarrow e^{+}e^{-}");
    leg->AddEntry(graph_stat_cocktail, "Cocktail");
    leg->Draw("same");


    counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = hsum->GetBinCenter(ibin);
        y[counter] = h1d_dca_FG->GetBinContent(h1d_dca_FG->FindBin(x[counter]))/hsum->GetBinContent(ibin);

        x_lohi[counter] = hsum->GetBinWidth(ibin)/2;

        stat_err[counter] = y[counter]*sqrt( pow(h1d_dca_FG->GetBinError(h1d_dca_FG->FindBin(x[counter]))/h1d_dca_FG->GetBinContent(h1d_dca_FG->FindBin(x[counter])), 2) + pow(hsum->GetBinError(ibin)/hsum->GetBinContent(ibin), 2) );
  
        double sys_fg = fcombined->Eval(x[counter])*h1d_dca_FG->GetBinContent(ibin);
        double sys_jpsi = fcombined_jpsi->Eval(x[counter])*h1d_dca_jpsi->GetBinContent(ibin);
        double sys_bg = h1d_mlp_dataBG->GetBinContent(ibin)*0.05212;
        double sys_hf = h1d_dca_hf_sys->GetBinError(ibin);

        double sys_err_cocktail = sqrt( pow(sys_jpsi, 2) +  pow(sys_bg, 2) + pow(sys_hf, 2) );

        sys_err[counter] = y[counter]*sqrt( pow(sys_fg/h1d_dca_FG->GetBinContent(ibin), 2) + pow(sys_err_cocktail/hsum->GetBinContent(ibin), 2));

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