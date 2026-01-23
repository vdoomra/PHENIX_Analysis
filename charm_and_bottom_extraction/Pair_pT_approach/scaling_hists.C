const int nbins = 6;
double pt_bins[nbins+1] = {0.0, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0};

const double nevents = 4.216682807*pow(10,11)*(0.79);

double mass_lo = 0.87;
double mass_hi = 4.50;

void scaling_hists(){


    TFile* f1 = new TFile("pair_pt.root","READ");
    TH1D* h1d_pair_pt_FG = (TH1D*)f1->Get("h1d_signal");

    TFile* f2 = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_ccbar_bbbar_SoftQCD_inelastic_75925_gencut_with_correction.root","READ");

    TH2D* h2d_FG_ccbar = (TH2D*)f2->Get("h2d_ee_FG_unlike_ccbar");
    h2d_FG_ccbar->SetName("h2d_FG_ccbar");

    h2d_FG_ccbar->Scale(42./75925.);
    h2d_FG_ccbar->Scale(4.8572);

    h2d_FG_ccbar->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_pair_pt_ccbar_temp = (TH1D*)h2d_FG_ccbar->ProjectionY();
    h1d_pair_pt_ccbar_temp->SetName("h1d_pair_pt_ccbar_temp");
    h1d_pair_pt_ccbar_temp->Scale(10.);

    TH1D* h1d_pair_pt_ccbar = new TH1D("h1d_pair_pt_ccbar", "Pair pT Distribution", nbins, pt_bins);
    h1d_pair_pt_ccbar->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_pair_pt_ccbar_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_pair_pt_ccbar_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_pair_pt_ccbar_temp->GetBinContent(i)*h1d_pair_pt_ccbar_temp->GetBinWidth(i);
            bin_error += pow(h1d_pair_pt_ccbar_temp->GetBinError(i)*h1d_pair_pt_ccbar_temp->GetBinWidth(i), 2);

        }

        h1d_pair_pt_ccbar->SetBinContent(ibin+1, bin_content/h1d_pair_pt_ccbar->GetBinWidth(ibin+1));
        h1d_pair_pt_ccbar->SetBinError(ibin+1, sqrt(bin_error)/h1d_pair_pt_ccbar->GetBinWidth(ibin+1));


    }

    TH2D* h2d_FG_bbbar = (TH2D*)f2->Get("h2d_ee_FG_unlike_bbbar");
    h2d_FG_bbbar->SetName("h2d_FG_bbbar");
    h2d_FG_bbbar->Scale(42./75925.);
    h2d_FG_bbbar->Scale(4.8572);

    h2d_FG_bbbar->GetXaxis()->SetRangeUser(mass_lo, mass_hi);
    TH1D* h1d_pair_pt_bbbar_temp = (TH1D*)h2d_FG_bbbar->ProjectionY();
    h1d_pair_pt_bbbar_temp->SetName("h1d_pair_pt_bbbar_temp");
    h1d_pair_pt_bbbar_temp->Scale(10.);

    TH1D* h1d_pair_pt_bbbar = new TH1D("h1d_pair_pt_bbbar", "Pair pT Distribution", nbins, pt_bins);
    h1d_pair_pt_bbbar->Sumw2();

    for(int ibin=0; ibin < nbins; ibin++){

        int bin_lo = h1d_pair_pt_bbbar_temp->FindBin(pt_bins[ibin] + 0.01);
        int bin_hi = h1d_pair_pt_bbbar_temp->FindBin(pt_bins[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_pair_pt_bbbar_temp->GetBinContent(i)*h1d_pair_pt_bbbar_temp->GetBinWidth(i);
            bin_error += pow(h1d_pair_pt_bbbar_temp->GetBinError(i)*h1d_pair_pt_bbbar_temp->GetBinWidth(i), 2);

        }

        h1d_pair_pt_bbbar->SetBinContent(ibin+1, bin_content/h1d_pair_pt_bbbar->GetBinWidth(ibin+1));
        h1d_pair_pt_bbbar->SetBinError(ibin+1, sqrt(bin_error)/h1d_pair_pt_bbbar->GetBinWidth(ibin+1));


    }

    TCanvas* c1 = new TCanvas();
    h1d_pair_pt_ccbar->Draw();
    h1d_pair_pt_bbbar->Draw("same");

    TCanvas* c2 = new TCanvas();

    TFile* f3 = new TFile("template_fits_pair_pt.root","READ");
    TH1D* h1d_para_ccbar = (TH1D*)f3->Get("h1d_para_ccbar");
    TH1D* h1d_para_bbbar = (TH1D*)f3->Get("h1d_para_bbbar");

    TF1* func_ccbar = new TF1("func_ccbar", "gaus", h1d_para_ccbar->GetMean()-1.5*h1d_para_ccbar->GetStdDev(), h1d_para_ccbar->GetMean()+1.5*h1d_para_ccbar->GetStdDev());
    TF1* func_bbbar = new TF1("func_bbbar", "gaus", h1d_para_bbbar->GetMean()-1.5*h1d_para_bbbar->GetStdDev(), h1d_para_bbbar->GetMean()+1.5*h1d_para_bbbar->GetStdDev());

    h1d_para_ccbar->Fit(func_ccbar, "R+");
    h1d_para_bbbar->Fit(func_bbbar, "R+");

    double scale_ccbar = func_ccbar->GetParameter(1)*(h1d_pair_pt_FG->Integral("width")/h1d_pair_pt_ccbar->Integral("width"));
    double scale_error_ccbar = func_ccbar->GetParameter(2)*(h1d_pair_pt_FG->Integral("width")/h1d_pair_pt_ccbar->Integral("width"));

    double scale_bbbar = func_bbbar->GetParameter(1)*(h1d_pair_pt_FG->Integral("width")/h1d_pair_pt_bbbar->Integral("width"));
    double scale_error_bbbar = func_bbbar->GetParameter(2)*(h1d_pair_pt_FG->Integral("width")/h1d_pair_pt_bbbar->Integral("width"));

    cout << scale_ccbar << "\t" << scale_error_ccbar << "\t" << h1d_pair_pt_ccbar->Integral("width")/h1d_pair_pt_FG->Integral("width") << endl;
    cout << scale_bbbar << "\t" << scale_error_bbbar << "\t" << h1d_pair_pt_bbbar->Integral("width")/h1d_pair_pt_FG->Integral("width") << endl;


    TH2D* h2d_delta_phi_pair_mass_ccbar = (TH2D*)f2->Get("h2d_delta_phi_pair_mass_unlike_ccbar");
    h2d_delta_phi_pair_mass_ccbar->Scale(42./75925.);
    h2d_delta_phi_pair_mass_ccbar->Scale(4.8572);

    TH2D* h2d_delta_phi_pair_mass_ccbar_central = (TH2D*)h2d_delta_phi_pair_mass_ccbar->Clone();
    h2d_delta_phi_pair_mass_ccbar_central->SetName("h2d_delta_phi_pair_mass_ccbar_central");
    h2d_delta_phi_pair_mass_ccbar_central->Scale(scale_ccbar);

    TH2D* h2d_delta_phi_pair_mass_ccbar_lower = (TH2D*)h2d_delta_phi_pair_mass_ccbar->Clone();
    h2d_delta_phi_pair_mass_ccbar_lower->SetName("h2d_delta_phi_pair_mass_ccbar_lower");
    h2d_delta_phi_pair_mass_ccbar_lower->Scale(scale_ccbar - scale_error_ccbar);

    TH2D* h2d_delta_phi_pair_mass_ccbar_upper = (TH2D*)h2d_delta_phi_pair_mass_ccbar->Clone();
    h2d_delta_phi_pair_mass_ccbar_upper->SetName("h2d_delta_phi_pair_mass_ccbar_upper");
    h2d_delta_phi_pair_mass_ccbar_upper->Scale(scale_ccbar + scale_error_ccbar);

    TH2D* h2d_delta_phi_pair_mass_bbbar = (TH2D*)f2->Get("h2d_delta_phi_pair_mass_unlike_bbbar");
    h2d_delta_phi_pair_mass_bbbar->Scale(42./75925.);
    h2d_delta_phi_pair_mass_bbbar->Scale(4.8572);

    TH2D* h2d_delta_phi_pair_mass_bbbar_central = (TH2D*)h2d_delta_phi_pair_mass_bbbar->Clone();
    h2d_delta_phi_pair_mass_bbbar_central->SetName("h2d_delta_phi_pair_mass_bbbar_central");
    h2d_delta_phi_pair_mass_bbbar_central->Scale(scale_bbbar);

    TH2D* h2d_delta_phi_pair_mass_bbbar_lower = (TH2D*)h2d_delta_phi_pair_mass_bbbar->Clone();
    h2d_delta_phi_pair_mass_bbbar_lower->SetName("h2d_delta_phi_pair_mass_bbbar_lower");
    h2d_delta_phi_pair_mass_bbbar_lower->Scale(scale_bbbar - scale_error_bbbar);

    TH2D* h2d_delta_phi_pair_mass_bbbar_upper = (TH2D*)h2d_delta_phi_pair_mass_bbbar->Clone();
    h2d_delta_phi_pair_mass_bbbar_upper->SetName("h2d_delta_phi_pair_mass_bbbar_upper");
    h2d_delta_phi_pair_mass_bbbar_upper->Scale(scale_bbbar + scale_error_bbbar);

    TH2D* hist_ccbar = (TH2D*)f2->Get("hist_unlike_ccbar");
    hist_ccbar->Scale(42./75925.);
    hist_ccbar->Scale(4.8572);

    TH2D* hist_ccbar_central = (TH2D*)hist_ccbar->Clone();
    hist_ccbar_central->SetName("hist_ccbar_central");
    hist_ccbar_central->Scale(scale_ccbar);

    TH2D* hist_ccbar_lower = (TH2D*)hist_ccbar->Clone();
    hist_ccbar_lower->SetName("hist_ccbar_lower");
    hist_ccbar_lower->Scale(scale_ccbar - scale_error_ccbar);

    TH2D* hist_ccbar_upper = (TH2D*)hist_ccbar->Clone();
    hist_ccbar_upper->SetName("hist_ccbar_upper");
    hist_ccbar_upper->Scale(scale_ccbar + scale_error_ccbar);

    TH2D* hist_bbbar = (TH2D*)f2->Get("hist_unlike_bbbar");
    hist_bbbar->Scale(42./75925.);
    hist_bbbar->Scale(4.8572);

    TH2D* hist_bbbar_central = (TH2D*)hist_bbbar->Clone();
    hist_bbbar_central->SetName("hist_bbbar_central");
    hist_bbbar_central->Scale(scale_bbbar);

    TH2D* hist_bbbar_lower = (TH2D*)hist_bbbar->Clone();
    hist_bbbar_lower->SetName("hist_bbbar_lower");
    hist_bbbar_lower->Scale(scale_bbbar - scale_error_bbbar);

    TH2D* hist_bbbar_upper = (TH2D*)hist_bbbar->Clone();
    hist_bbbar_upper->SetName("hist_bbbar_upper");
    hist_bbbar_upper->Scale(scale_bbbar + scale_error_bbbar);

    TFile* f4 = new TFile("/Users/vassu/Desktop/new_files/data_files/producing_templates_without_20percent_reduction.root","READ");

    TH1D* h1d_ee_FG_ccbar_central = (TH1D*)f4->Get("h1d_ee_FG_ccbar");
    h1d_ee_FG_ccbar_central->SetName("h1d_ee_FG_ccbar_central");
    h1d_ee_FG_ccbar_central->Scale(scale_ccbar);

    TH1D* h1d_ee_FG_ccbar_lower = (TH1D*)f4->Get("h1d_ee_FG_ccbar");
    h1d_ee_FG_ccbar_lower->SetName("h1d_ee_FG_ccbar_lower");
    h1d_ee_FG_ccbar_lower->Scale(scale_ccbar - scale_error_ccbar);

    TH1D* h1d_ee_FG_ccbar_upper = (TH1D*)f4->Get("h1d_ee_FG_ccbar");
    h1d_ee_FG_ccbar_upper->SetName("h1d_ee_FG_ccbar_upper");
    h1d_ee_FG_ccbar_upper->Scale(scale_ccbar + scale_error_ccbar);

    TH1D* h1d_ee_FG_bbbar_central = (TH1D*)f4->Get("h1d_ee_FG_bbbar");
    h1d_ee_FG_bbbar_central->SetName("h1d_ee_FG_bbbar_central");
    h1d_ee_FG_bbbar_central->Scale(scale_bbbar);

    TH1D* h1d_ee_FG_bbbar_lower = (TH1D*)f4->Get("h1d_ee_FG_bbbar");
    h1d_ee_FG_bbbar_lower->SetName("h1d_ee_FG_bbbar_lower");
    h1d_ee_FG_bbbar_lower->Scale(scale_bbbar - scale_error_bbbar);

    TH1D* h1d_ee_FG_bbbar_upper = (TH1D*)f4->Get("h1d_ee_FG_bbbar");
    h1d_ee_FG_bbbar_upper->SetName("h1d_ee_FG_bbbar_upper");
    h1d_ee_FG_bbbar_upper->Scale(scale_bbbar + scale_error_bbbar);

    TFile* fout = new TFile("scaling_hists.root","RECREATE");
    fout->cd();
   
    h2d_delta_phi_pair_mass_ccbar_central->Write();
    h2d_delta_phi_pair_mass_ccbar_lower->Write();
    h2d_delta_phi_pair_mass_ccbar_upper->Write();

    h2d_delta_phi_pair_mass_bbbar_central->Write();
    h2d_delta_phi_pair_mass_bbbar_lower->Write();
    h2d_delta_phi_pair_mass_bbbar_upper->Write();

    h1d_ee_FG_ccbar_central->Write();
    h1d_ee_FG_ccbar_lower->Write();
    h1d_ee_FG_ccbar_upper->Write();

    h1d_ee_FG_bbbar_central->Write();
    h1d_ee_FG_bbbar_lower->Write();
    h1d_ee_FG_bbbar_upper->Write();

    hist_ccbar_central->Write();
    hist_ccbar_lower->Write();
    hist_ccbar_upper->Write();

    hist_bbbar_central->Write();
    hist_bbbar_lower->Write();
    hist_bbbar_upper->Write();

    fout->Close();
    
}