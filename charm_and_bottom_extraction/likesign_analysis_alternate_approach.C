const int nbins_mass = 11;
const double mass_range[nbins_mass+1] = { 0.80, 0.95, 1.15, 1.45, 1.75, 2.05, 2.40, 2.80, 3.20, 3.45, 4.00, 4.5};

const int nbins_mass_bbbar = 3;
const double mass_range_bbbar[nbins_mass_bbbar+1] = {2.80, 3.45, 4.00, 4.5};

const double nevents = 4.216682807*pow(10,11)*(0.79);

double combinedFunction(double *x, double *par) {

    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    if(x[0] > 3.50) fERT_val = 0.0;
    double fn0_val = par[3];
    double feoverp_val = par[4] + par[5]*x[0] + par[6]*x[0]*x[0];

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}

void likesign_analysis_alternate_approach(const int stat){

    TFile* f1 = new TFile("plotting_BG_with_uncertainties.root","READ");
    TH1D* h1d_mlp_dataBG_unlike_temp = (TH1D*)f1->Get("h1d_mlp_dataBG");
    h1d_mlp_dataBG_unlike_temp->SetName("h1d_mlp_dataBG_unlike_temp");

    TH1D* h1d_mlp_dataBG_unlike = new TH1D("h1d_mlp_dataBG_unlike","Mass Distribution", nbins_mass, mass_range);
    h1d_mlp_dataBG_unlike->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_mlp_dataBG_unlike_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_mlp_dataBG_unlike_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_mlp_dataBG_unlike_temp->GetBinContent(i)*h1d_mlp_dataBG_unlike_temp->GetBinWidth(i);
            bin_error += pow(h1d_mlp_dataBG_unlike_temp->GetBinError(i)*h1d_mlp_dataBG_unlike_temp->GetBinWidth(i), 2);

        }

        h1d_mlp_dataBG_unlike->SetBinContent(ibin+1, bin_content/h1d_mlp_dataBG_unlike->GetBinWidth(ibin+1));
        h1d_mlp_dataBG_unlike->SetBinError(ibin+1, sqrt(bin_error)/h1d_mlp_dataBG_unlike->GetBinWidth(ibin+1));


    }

    TFile* f2 = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_BG_phi_redef_n0_2_with_correction.root","READ");
    TH1D* h1d_BG_unlike = (TH1D*)f2->Get("h1d_ee_FG_unlike");
    h1d_BG_unlike->SetName("h1d_BG_unlike");

    TH1D* h1d_BG_like = (TH1D*)f2->Get("h1d_ee_FG_like");
    h1d_BG_like->SetName("h1d_BG_like");

    TH1D* h1d_correction_temp = (TH1D*)h1d_BG_unlike->Clone();
    h1d_correction_temp->SetName("h1d_correction_temp");
    h1d_correction_temp->Divide(h1d_BG_like);

    TCanvas* c_correction = new TCanvas();
    c_correction->cd();
    h1d_correction_temp->Draw();

    TH1D* h1d_correction = new TH1D("h1d_correction","DCA Distribution", nbins_mass, mass_range);
    h1d_correction->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        double bin_center = (mass_range[ibin] + mass_range[ibin+1])/2.;

        if(bin_center < 2.60){

            h1d_correction->SetBinContent(ibin+1, h1d_correction_temp->GetBinContent(h1d_correction_temp->FindBin(bin_center)));
            h1d_correction->SetBinError(ibin+1, h1d_correction_temp->GetBinError(h1d_correction_temp->FindBin(bin_center)));

        }

        else {

            h1d_correction->SetBinContent(ibin+1, 0.959819);
            h1d_correction->SetBinError(ibin+1, 0.00359);

        }

        
    }

    TH1D* h1d_mlp_dataBG_like = (TH1D*)h1d_mlp_dataBG_unlike->Clone();
    h1d_mlp_dataBG_like->SetName("h1d_mlp_dataBG_like");
    h1d_mlp_dataBG_like->Divide(h1d_correction);

    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
    h1d_mlp_dataBG_unlike->Draw();
    h1d_mlp_dataBG_like->Draw("same");

    TFile* f3 = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_FG_phi_redef_n0_2_with_correction.root","READ");
    TH1D* h1d_like_FG_temp = (TH1D*)f3->Get("h1d_ee_FG_like");
    h1d_like_FG_temp->SetName("h1d_like_FG_temp");

    TH1D* h1d_like_FG = new TH1D("h1d_like_FG","Mass Distribution", nbins_mass, mass_range);
    h1d_like_FG->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_like_FG_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_like_FG_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_like_FG_temp->GetBinContent(i)*h1d_like_FG_temp->GetBinWidth(i);
            bin_error += pow(h1d_like_FG_temp->GetBinError(i)*h1d_like_FG_temp->GetBinWidth(i), 2);

        }

        h1d_like_FG->SetBinContent(ibin+1, bin_content/h1d_like_FG->GetBinWidth(ibin+1));
        h1d_like_FG->SetBinError(ibin+1, sqrt(bin_error)/h1d_like_FG->GetBinWidth(ibin+1));


    }

    h1d_like_FG->Scale(1./nevents);
    h1d_like_FG->Scale(1.0726*0.545*42);

    TF1* fcombined = new TF1("fcombined", combinedFunction, 0.0, 5.5, 7);
    fcombined->SetParameters(0.08895, 0.0334, -0.0137, 0.06154, 0.0837, 0.01724, -0.00625);

    TH1D* h1d_signal_temp = (TH1D*)h1d_like_FG->Clone();
    h1d_signal_temp->SetName("h1d_signal_temp");
    h1d_signal_temp->Reset("ICESM");

    for(int ibin = 1; ibin < nbins_mass+1; ibin++){

        double bin_content = h1d_like_FG->GetBinContent(ibin) - h1d_mlp_dataBG_like->GetBinContent(ibin);
        h1d_signal_temp->SetBinContent(ibin, bin_content);
        h1d_signal_temp->SetBinError(ibin, sqrt( pow(h1d_like_FG->GetBinError(ibin), 2) + pow(h1d_mlp_dataBG_like->GetBinError(ibin), 2) ));


    }

    TCanvas* c2 = new TCanvas();
    c2->SetLogy();
    h1d_like_FG->Draw();
    h1d_mlp_dataBG_like->Draw("same");

    TCanvas* c3 = new TCanvas();
    c3->SetLogy();
    h1d_signal_temp->Draw();

    TH1D* h1d_signal_stat = new TH1D("h1d_signal_stat","Mass Distribution", nbins_mass_bbbar, mass_range_bbbar);
    h1d_signal_stat->Sumw2();

    TH1D* h1d_signal_sys = new TH1D("h1d_signal_sys","Mass Distribution", nbins_mass_bbbar, mass_range_bbbar);
    h1d_signal_sys->Sumw2();

    for(int ibin=0; ibin < nbins_mass_bbbar; ibin++){

        int bin_lo = h1d_signal_temp->FindBin(mass_range_bbbar[ibin] + 0.001);
        int bin_hi = h1d_signal_temp->FindBin(mass_range_bbbar[ibin+1] - 0.001);

        double bin_center = (mass_range_bbbar[ibin] + mass_range_bbbar[ibin+1])/2;

        double bin_content = 0; double bin_error = 0;
        double bin_content_fg = 0; double bin_content_bg = 0;
        double bin_error_stat = 0; double bin_error_sys = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_signal_temp->GetBinContent(i)*h1d_signal_temp->GetBinWidth(i);
            bin_content_fg += h1d_like_FG->GetBinContent(i)*h1d_like_FG->GetBinWidth(i);
            bin_content_bg += h1d_mlp_dataBG_like->GetBinContent(i)*h1d_mlp_dataBG_like->GetBinWidth(i);
            bin_error_stat += pow(h1d_signal_temp->GetBinError(i)*h1d_signal_temp->GetBinWidth(i), 2);

        }

        double fg_error = bin_content_fg*fcombined->Eval(bin_center);
        double bg_error = bin_content_bg*sqrt( 0.02*0.02 + 0.01*0.01 + 0.05212*0.05212 );
            
        bin_error_sys = sqrt( pow(fg_error, 2) + pow(bg_error, 2));

        h1d_signal_stat->SetBinContent(ibin+1, bin_content/h1d_signal_stat->GetBinWidth(ibin+1));
        h1d_signal_sys->SetBinContent(ibin+1, bin_content/h1d_signal_sys->GetBinWidth(ibin+1));
        h1d_signal_stat->SetBinError(ibin+1, sqrt(bin_error_stat)/h1d_signal_stat->GetBinWidth(ibin+1));
        h1d_signal_sys->SetBinError(ibin+1, bin_error_sys/h1d_signal_sys->GetBinWidth(ibin+1));


    }

    TFile* f5 = new TFile("/Users/vassu/Desktop/new_files/data_files/producing_templates_without_20percent_reduction.root","READ");
    TH1D* h1d_bbbar_pythia_like_temp = (TH1D*)f5->Get("h1d_ee_FG_bbbar_like");
    h1d_bbbar_pythia_like_temp->SetName("h1d_bbbar_pythia_like_temp");

    TH1D* h1d_bbbar_pythia_like = new TH1D("h1d_bbbar_pythia_like","Mass Distribution", nbins_mass_bbbar, mass_range_bbbar);
    h1d_bbbar_pythia_like->Sumw2();

    for(int ibin=0; ibin < nbins_mass_bbbar; ibin++){

        int bin_lo = h1d_bbbar_pythia_like_temp->FindBin(mass_range_bbbar[ibin] + 0.001);
        int bin_hi = h1d_bbbar_pythia_like_temp->FindBin(mass_range_bbbar[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_bbbar_pythia_like_temp->GetBinContent(i)*h1d_bbbar_pythia_like_temp->GetBinWidth(i);
            bin_error += pow(h1d_bbbar_pythia_like_temp->GetBinError(i)*h1d_bbbar_pythia_like_temp->GetBinWidth(i), 2);

        }

        h1d_bbbar_pythia_like->SetBinContent(ibin+1, bin_content/h1d_bbbar_pythia_like->GetBinWidth(ibin+1));
        h1d_bbbar_pythia_like->SetBinError(ibin+1, 0.0);


    }

    h1d_bbbar_pythia_like->Draw();

    double scale, scale_error;

    if(stat){

        double h1d_signal_error, h1d_bbbar_pythia_like_error;
        double S = h1d_signal_stat->IntegralAndError(1, h1d_signal_stat->GetNbinsX(), h1d_signal_error, "width");
        double B = h1d_bbbar_pythia_like->IntegralAndError(1, h1d_bbbar_pythia_like->GetNbinsX(), h1d_bbbar_pythia_like_error, "width");

        scale = S/B;
        scale_error = scale*sqrt(pow(h1d_signal_error/S, 2) + pow(h1d_bbbar_pythia_like_error/B, 2));

    }

    else{

        TH1D* hist = new TH1D("hist", "Distribution", 200, -10, 20);

        for(int itr = 0; itr < 2e5; itr++){

            double r = gRandom->Gaus(0, 1);

            TH1D* h1d_signal_iterative = new TH1D(Form("h1d_signal_iterative_%d", itr),"Mass Distribution", nbins_mass_bbbar, mass_range_bbbar);
            h1d_signal_iterative->Sumw2();
        
            for (int i = 1; i < h1d_signal_iterative->GetNbinsX()+1; ++i) {
    
                double new_y = (h1d_signal_sys->GetBinContent(i) + r * h1d_signal_sys->GetBinError(i));
                h1d_signal_iterative->SetBinContent(i, new_y);
                h1d_signal_iterative->SetBinError(i, h1d_signal_stat->GetBinError(i));
    
            }
    
            double scale = h1d_signal_iterative->Integral("width")/h1d_bbbar_pythia_like->Integral("width");
    
            hist->Fill(scale);
    
        }

        TF1* func = new TF1("func","gaus", hist->GetMean()-2*hist->GetStdDev(), hist->GetMean()+2*hist->GetStdDev());
        hist->Fit(func,"R+");

        scale = func->GetParameter(1);
        scale_error = func->GetParameter(2);

    }

    cout << scale << "\t" << scale_error << endl;

    TH1D* h1d_bbbar_pythia_like_scaled_central = (TH1D*)h1d_bbbar_pythia_like->Clone();
    h1d_bbbar_pythia_like_scaled_central->SetName("h1d_bbbar_pythia_like_scaled_central");
    h1d_bbbar_pythia_like_scaled_central->Scale(scale);

    TH1D* h1d_bbbar_pythia_like_scaled_lower = (TH1D*)h1d_bbbar_pythia_like->Clone();
    h1d_bbbar_pythia_like_scaled_lower->SetName("h1d_bbbar_pythia_like_scaled_lower");
    h1d_bbbar_pythia_like_scaled_lower->Scale(scale-scale_error);

    TH1D* h1d_bbbar_pythia_like_scaled_upper = (TH1D*)h1d_bbbar_pythia_like->Clone();
    h1d_bbbar_pythia_like_scaled_upper->SetName("h1d_bbbar_pythia_like_scaled_upper");
    h1d_bbbar_pythia_like_scaled_upper->Scale(scale+scale_error);


    TH1D* h1d_bbbar_pythia_unlike_dataBinning_central = (TH1D*)f5->Get("h1d_ee_FG_bbbar");
    h1d_bbbar_pythia_unlike_dataBinning_central->SetName("h1d_bbbar_pythia_unlike_dataBinning_central");
    h1d_bbbar_pythia_unlike_dataBinning_central->Scale(scale);

    TH1D* h1d_bbbar_pythia_unlike_dataBinning_lower = (TH1D*)f5->Get("h1d_ee_FG_bbbar");
    h1d_bbbar_pythia_unlike_dataBinning_lower->SetName("h1d_bbbar_pythia_unlike_dataBinning_lower");
    h1d_bbbar_pythia_unlike_dataBinning_lower->Scale(scale-scale_error);

    TH1D* h1d_bbbar_pythia_unlike_dataBinning_upper = (TH1D*)f5->Get("h1d_ee_FG_bbbar");
    h1d_bbbar_pythia_unlike_dataBinning_upper->SetName("h1d_bbbar_pythia_unlike_dataBinning_upper");
    h1d_bbbar_pythia_unlike_dataBinning_upper->Scale(scale+scale_error);

    TCanvas* c5 = new TCanvas();
    c5->SetLogy();
    h1d_signal_stat->Draw();
    h1d_bbbar_pythia_like_scaled_central->Draw("same hist");
    h1d_bbbar_pythia_like_scaled_lower->Draw("same hist");
    h1d_bbbar_pythia_like_scaled_upper->Draw("same hist");



    TGraphAsymmErrors* band = new TGraphAsymmErrors((int)h1d_bbbar_pythia_unlike_dataBinning_central->GetNbinsX());
    band->SetName("band");

    for (int i = 1; i <= (int)h1d_bbbar_pythia_unlike_dataBinning_central->GetNbinsX(); ++i) {

        double x = h1d_bbbar_pythia_unlike_dataBinning_central->GetBinCenter(i);
        double y = h1d_bbbar_pythia_unlike_dataBinning_central->GetBinContent(i);

        double y_up   = h1d_bbbar_pythia_unlike_dataBinning_upper->GetBinContent(i);
        double y_down = h1d_bbbar_pythia_unlike_dataBinning_lower->GetBinContent(i);

        double ey_high = y_up   - y;
        double ey_low  = y      - y_down;

        band->SetPoint(i - 1, x, y);
        band->SetPointError(i - 1, 0, 0, ey_low, ey_high);

    }



    TCanvas* c6 = new TCanvas();
    c6->SetLogy();
    band->SetFillColorAlpha(kBlue, 0.3);
    band->SetLineColor(kBlue);
    band->Draw("A3");

    TFile* fout = new TFile("dummy.root","RECREATE");
    fout->cd();
    h1d_signal_stat->Write();
    h1d_signal_sys->Write();
    band->Write();
    h1d_bbbar_pythia_like->Write();

    h1d_bbbar_pythia_like_scaled_central->Write();
    h1d_bbbar_pythia_like_scaled_lower->Write();
    h1d_bbbar_pythia_like_scaled_upper->Write();

    h1d_bbbar_pythia_unlike_dataBinning_central->Write();
    h1d_bbbar_pythia_unlike_dataBinning_lower->Write();
    h1d_bbbar_pythia_unlike_dataBinning_upper->Write();
    fout->Close();


}