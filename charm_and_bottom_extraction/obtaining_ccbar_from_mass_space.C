
const int nbins_mass = 26;
const double mass_range[nbins_mass+1] = {0.87, 0.95, 1.05, 1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.05, 3.10, 3.15, 3.20, 3.30, 3.45, 3.60, 3.70, 3.80, 4.00, 4.25, 4.5};

const double nevents = 4.216682807*pow(10,11)*(0.79);

const int nbins_mass_ccbar = 4;
const double mass_range_ccbar[nbins_mass_ccbar+1] =  { 1.90, 2.20, 2.80, 3.10, 3.45};

const int npts = 40;
double x_lohi[npts] = {-999};
double x[npts] = {-999};
double y[npts] = {-999};
double stat_err[npts] = {-999};
double sys_err[npts] = {-999};

double combinedFunction(double *x, double *par) {

    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    if(x[0] > 3.50) fERT_val = 0.0;
    double fn0_val = par[3];
    double feoverp_val = par[4] + par[5]*x[0] + par[6]*x[0]*x[0];
    if(x[0] > 3.20) feoverp_val = 0.03262;

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}

void obtaining_ccbar_from_mass_space(const int stat){

    TFile* f1 = new TFile("/Users/vassu/Desktop/new_files/data_files/DCA_data_FG_phi_redef_n0_2_with_correction.root","READ");
    TH1D* h1d_dataFG_temp = (TH1D*)f1->Get("h1d_ee_FG_unlike");
    h1d_dataFG_temp->Scale(1./nevents);
    h1d_dataFG_temp->Scale(1.0726*0.545*42);

    TH1D* h1d_dataFG = new TH1D("h1d_dataFG","Mass Distribution", nbins_mass, mass_range);
    h1d_dataFG->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_dataFG_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_dataFG_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dataFG_temp->GetBinContent(i)*h1d_dataFG_temp->GetBinWidth(i);
            bin_error += pow(h1d_dataFG_temp->GetBinError(i)*h1d_dataFG_temp->GetBinWidth(i), 2);

        }

        h1d_dataFG->SetBinContent(ibin+1, bin_content/h1d_dataFG->GetBinWidth(ibin+1));
        h1d_dataFG->SetBinError(ibin+1, sqrt(bin_error)/h1d_dataFG->GetBinWidth(ibin+1));


    }

    TFile* f2 = new TFile("plotting_BG_with_uncertainties.root","READ");
    TH1D* h1d_dataBG_temp = (TH1D*)f2->Get("h1d_mlp_dataBG");

    TH1D* h1d_dataBG = new TH1D("h1d_dataBG","Mass Distribution", nbins_mass, mass_range);
    h1d_dataBG->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_dataBG_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_dataBG_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_dataBG_temp->GetBinContent(i)*h1d_dataBG_temp->GetBinWidth(i);
            bin_error += pow(h1d_dataBG_temp->GetBinError(i)*h1d_dataBG_temp->GetBinWidth(i), 2);

        }

        h1d_dataBG->SetBinContent(ibin+1, bin_content/h1d_dataBG->GetBinWidth(ibin+1));
        h1d_dataBG->SetBinError(ibin+1, sqrt(bin_error)/h1d_dataBG->GetBinWidth(ibin+1));


    }

    TFile* fjpsi = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_jpsi_2k_gencut_with_correction.root","READ");
    TH1D* h1d_ee_FG_jpsi_temp = (TH1D*)fjpsi->Get("h1d_ee_FG_unlike");
    h1d_ee_FG_jpsi_temp->SetName("h1d_ee_FG_jpsi_temp");
    h1d_ee_FG_jpsi_temp->Scale(42.*0.05971);
    h1d_ee_FG_jpsi_temp->Scale(1./(2000.*25000.));

    TH1D* h1d_ee_FG_jpsi = new TH1D("h1d_ee_FG_jpsi","Mass Distribution", nbins_mass, mass_range);
    h1d_ee_FG_jpsi->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_ee_FG_jpsi_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_ee_FG_jpsi_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_ee_FG_jpsi_temp->GetBinContent(i)*h1d_ee_FG_jpsi_temp->GetBinWidth(i);
            bin_error += pow(h1d_ee_FG_jpsi_temp->GetBinError(i)*h1d_ee_FG_jpsi_temp->GetBinWidth(i), 2);

        }

        h1d_ee_FG_jpsi->SetBinContent(ibin+1, bin_content/h1d_ee_FG_jpsi->GetBinWidth(ibin+1));
        h1d_ee_FG_jpsi->SetBinError(ibin+1, sqrt(bin_error)/h1d_ee_FG_jpsi->GetBinWidth(ibin+1));


    }

    TFile* fpsip = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_psip_2k_gencut_with_correction.root","READ");
    TH1D* h1d_ee_FG_psip_temp = (TH1D*)fpsip->Get("h1d_ee_FG_unlike");
    h1d_ee_FG_psip_temp->SetName("h1d_ee_FG_psip_temp");
    h1d_ee_FG_psip_temp->Scale(42.*(7.93e-3));
    h1d_ee_FG_psip_temp->Scale(1./(2000.*25000.));

    TH1D* h1d_ee_FG_psip = new TH1D("h1d_ee_FG_psip","Mass Distribution", nbins_mass, mass_range);
    h1d_ee_FG_psip->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_ee_FG_psip_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_ee_FG_psip_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_ee_FG_psip_temp->GetBinContent(i)*h1d_ee_FG_psip_temp->GetBinWidth(i);
            bin_error += pow(h1d_ee_FG_psip_temp->GetBinError(i)*h1d_ee_FG_psip_temp->GetBinWidth(i), 2);

        }

        h1d_ee_FG_psip->SetBinContent(ibin+1, bin_content/h1d_ee_FG_psip->GetBinWidth(ibin+1));
        h1d_ee_FG_psip->SetBinError(ibin+1, sqrt(bin_error)/h1d_ee_FG_psip->GetBinWidth(ibin+1));


    }

    TFile* fbbbar = new TFile("normalizing_two_TGraphs_bbbar_stat.root","READ");
    TH1D* h1d_bbbar_central = (TH1D*)fbbbar->Get("h1d_bbbar_central");
    TH1D* h1d_bbbar_lower = (TH1D*)fbbbar->Get("h1d_bbbar_lower");
    TH1D* h1d_bbbar_upper = (TH1D*)fbbbar->Get("h1d_bbbar_upper");

    TH1D* h1d_ee_FG_bbbar = new TH1D("h1d_ee_FG_bbbar","Mass Distribution", nbins_mass, mass_range);
    h1d_ee_FG_bbbar->Sumw2();

    for(int ibin = 0; ibin < nbins_mass; ibin++){
        
        double bin_center = (mass_range[ibin] + mass_range[ibin+1])/2.;
        h1d_ee_FG_bbbar->SetBinContent(ibin+1, h1d_bbbar_central->GetBinContent(h1d_bbbar_central->FindBin(bin_center)));
        h1d_ee_FG_bbbar->SetBinError(ibin+1, h1d_bbbar_central->GetBinContent(h1d_bbbar_central->FindBin(bin_center)) - h1d_bbbar_lower->GetBinContent(h1d_bbbar_lower->FindBin(bin_center)));

    }

    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
    h1d_dataFG->Draw();
    h1d_dataBG->Draw("same");
    h1d_ee_FG_jpsi->Draw("same hist");
    h1d_ee_FG_psip->Draw("same hist");
    h1d_ee_FG_bbbar->Draw("same");

    TH1D* hsum = (TH1D*)h1d_dataBG->Clone();
    hsum->SetName("hsum");
    hsum->Add(h1d_ee_FG_jpsi);
    hsum->Add(h1d_ee_FG_psip);
    hsum->Add(h1d_ee_FG_bbbar);

    TCanvas* c2 = new TCanvas();
    c2->SetLogy();
    h1d_dataFG->Draw();
    hsum->Draw("same");

    TH1D* h1d_signal_temp = (TH1D*)h1d_dataFG->Clone();
    h1d_signal_temp->SetName("h1d_signal_temp");
    h1d_signal_temp->Reset("ICESM");

    for(int ibin = 1; ibin < nbins_mass+1; ibin++){

        double bin_content = h1d_dataFG->GetBinContent(ibin) - hsum->GetBinContent(ibin);
        h1d_signal_temp->SetBinContent(ibin, bin_content);
        h1d_signal_temp->SetBinError(ibin, sqrt( pow(h1d_dataFG->GetBinError(ibin), 2) + pow(h1d_dataBG->GetBinError(ibin), 2) + pow(h1d_ee_FG_bbbar->GetBinError(ibin), 2)));


    }


    TCanvas* c3 = new TCanvas();
    c3->SetLogy();
    h1d_signal_temp->Draw();

    TF1* func = new TF1("func","pol1", 1.0, 4.5);
    func->SetParameters(0.4891, -0.1212); // bbbar systematic parameters

    TF1* fcombined = new TF1("fcombined", combinedFunction, 0.0, 5.5, 7);
    fcombined->SetParameters(0.102, -0.0022, -0.00492, 0.06424, 0.06622, 0.08956, -0.03127);

    TH1D* h1d_signal = new TH1D("h1d_signal","Mass Distribution", nbins_mass_ccbar, mass_range_ccbar);
    h1d_signal->Sumw2();

    TH1D* h1d_signal_sys = new TH1D("h1d_signal_sys","Mass Distribution", nbins_mass_ccbar, mass_range_ccbar);
    h1d_signal_sys->Sumw2();

    for(int ibin=0; ibin < nbins_mass_ccbar; ibin++){

        int bin_lo = h1d_signal_temp->FindBin(mass_range_ccbar[ibin] + 0.001);
        int bin_hi = h1d_signal_temp->FindBin(mass_range_ccbar[ibin+1] - 0.001);

        double bin_center = (mass_range_ccbar[ibin] + mass_range_ccbar[ibin+1])/2.;

        double bin_content = 0; double bin_error = 0;
        double bin_content_fg = 0;
        double bin_content_bg = 0;
        double bin_content_bbbar = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_signal_temp->GetBinContent(i)*h1d_signal_temp->GetBinWidth(i);
            bin_error += pow(h1d_signal_temp->GetBinError(i)*h1d_signal_temp->GetBinWidth(i), 2);
            bin_content_fg += h1d_dataFG->GetBinContent(i)*h1d_dataFG->GetBinWidth(i);
            bin_content_bg += h1d_dataBG->GetBinContent(i)*h1d_dataBG->GetBinWidth(i);
            bin_content_bbbar += h1d_ee_FG_bbbar->GetBinContent(i)*h1d_ee_FG_bbbar->GetBinWidth(i);

        }

        double sys_err_fg = (fcombined->Eval(bin_center)*bin_content_fg)/h1d_signal->GetBinWidth(ibin+1);
        double sys_err_bg = (0.05212*bin_content_bg)/h1d_signal->GetBinWidth(ibin+1);
        double sys_err_bbbar = (func->Eval(bin_center)*bin_content_bbbar)/h1d_signal->GetBinWidth(ibin+1);

        double sys_err = sqrt( pow(sys_err_fg, 2) + pow(sys_err_bg, 2)); //I am not adding the bbbar systematics

        h1d_signal->SetBinContent(ibin+1, bin_content/h1d_signal->GetBinWidth(ibin+1));
        h1d_signal->SetBinError(ibin+1, sqrt(bin_error)/h1d_signal->GetBinWidth(ibin+1));

        h1d_signal_sys->SetBinContent(ibin+1, bin_content/h1d_signal->GetBinWidth(ibin+1));
        h1d_signal_sys->SetBinError(ibin+1, sys_err);


    }

    TCanvas* c4 = new TCanvas();
    c4->SetLogy();
    h1d_signal->Draw();

    TF1* fcorrection = new TF1("fcorrection","pol3", 1.0, 4.5);
    fcorrection->SetParameters(0.0591339, -0.0391315, 0.0196653, -0.00210255);

    int nbins = h1d_signal->GetNbinsX();
    int counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = h1d_signal->GetBinCenter(ibin);
        y[counter] = h1d_signal->GetBinContent(ibin)/fcorrection->Eval(x[counter]);
        x_lohi[counter] = 0.05;

        stat_err[counter] = h1d_signal->GetBinError(ibin)/fcorrection->Eval(x[counter]);
        sys_err[counter] = h1d_signal_sys->GetBinError(ibin)/fcorrection->Eval(x[counter]);

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

    TCanvas* c12 = new TCanvas();
    mg->Draw("a");

    mg->GetXaxis()->SetTitle("m_{ee}[GeV/c^{2}]");
    mg->GetYaxis()->SetTitle("#frac{1}{N_{events}}#frac{#epsilon_{BBC}#sigma_{pp}}{#epsilon_{bias}}#frac{dN}{dm_{ee}} [GeV/c^{2}]^{-1} [In PHENIX Acceptance]");

    

    TFile* fout = new TFile("obtaining_ccbar_from_mass_space_corrected.root","RECREATE");
    fout->cd();
    h1d_signal->Write();
    h1d_signal_sys->Write();
    graph_stat->Write();
    graph_sys->Write();
    fout->Close();

    TFile* f3 = new TFile("/Users/vassu/Desktop/new_files/data_files/producing_templates_without_20percent_reduction.root","READ");
    TH1D* h1d_ccbar_pythia_temp = (TH1D*)f3->Get("h1d_ee_FG_ccbar");
    h1d_ccbar_pythia_temp->SetName("h1d_ccbar_pythia_temp");

    TH1D* h1d_ccbar_pythia = new TH1D("h1d_ccbar_pythia","Mass Distribution", nbins_mass_ccbar, mass_range_ccbar);
    h1d_ccbar_pythia->Sumw2();

    for(int ibin=0; ibin < nbins_mass_ccbar; ibin++){

        int bin_lo = h1d_ccbar_pythia_temp->FindBin(mass_range_ccbar[ibin] + 0.001);
        int bin_hi = h1d_ccbar_pythia_temp->FindBin(mass_range_ccbar[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_ccbar_pythia_temp->GetBinContent(i)*h1d_ccbar_pythia_temp->GetBinWidth(i);
            bin_error += pow(h1d_ccbar_pythia_temp->GetBinError(i)*h1d_ccbar_pythia_temp->GetBinWidth(i), 2);

        }

        h1d_ccbar_pythia->SetBinContent(ibin+1, bin_content/h1d_ccbar_pythia->GetBinWidth(ibin+1));
        h1d_ccbar_pythia->SetBinError(ibin+1, 0.0);


    }

    double scale, scale_error;

    if(stat){

        double h1d_signal_error, h1d_ccbar_pythia_error;
        double S = h1d_signal->IntegralAndError(1, h1d_signal->GetNbinsX(), h1d_signal_error, "width");
        double B = h1d_ccbar_pythia->IntegralAndError(1, h1d_ccbar_pythia->GetNbinsX(), h1d_ccbar_pythia_error, "width");

        scale = S/B;
        scale_error = scale*sqrt(pow(h1d_signal_error/S, 2) + pow(h1d_ccbar_pythia_error/B, 2));

    }

    else{

        TH1D* hist = new TH1D("hist", "Distribution", 200, -10, 20);

        for(int itr = 0; itr < 2e5; itr++){

            double r = gRandom->Gaus(0, 1);

            TH1D* h1d_signal_iterative = new TH1D(Form("h1d_signal_iterative_%d", itr),"Mass Distribution", nbins_mass_ccbar, mass_range_ccbar);
            h1d_signal_iterative->Sumw2();
        
            for (int i = 1; i < h1d_signal_iterative->GetNbinsX()+1; ++i) {
    
                double new_y = (h1d_signal_sys->GetBinContent(i) + r * h1d_signal_sys->GetBinError(i));
                h1d_signal_iterative->SetBinContent(i, new_y);
                h1d_signal_iterative->SetBinError(i, h1d_signal->GetBinError(i));
    
            }
    
            double scale = h1d_signal_iterative->Integral("width")/h1d_ccbar_pythia->Integral("width");
    
            hist->Fill(scale);
    
        }

        TF1* func = new TF1("func","gaus", hist->GetMean()-2*hist->GetStdDev(), hist->GetMean()+2*hist->GetStdDev());
        hist->Fit(func,"R+");

        scale = func->GetParameter(1);
        scale_error = func->GetParameter(2);

    }

    cout << scale << "\t" << scale_error << endl;

    TH1D* h1d_ccbar_pythia_scaled_central = (TH1D*)h1d_ccbar_pythia->Clone();
    h1d_ccbar_pythia_scaled_central->SetName("h1d_ccbar_pythia_scaled_central");
    h1d_ccbar_pythia_scaled_central->Scale(scale);

    TH1D* h1d_ccbar_pythia_scaled_lower = (TH1D*)h1d_ccbar_pythia->Clone();
    h1d_ccbar_pythia_scaled_lower->SetName("h1d_ccbar_pythia_scaled_lower");
    h1d_ccbar_pythia_scaled_lower->Scale(scale-scale_error);

    TH1D* h1d_ccbar_pythia_scaled_upper = (TH1D*)h1d_ccbar_pythia->Clone();
    h1d_ccbar_pythia_scaled_upper->SetName("h1d_ccbar_pythia_scaled_upper");
    h1d_ccbar_pythia_scaled_upper->Scale(scale+scale_error);

    TCanvas* c5 = new TCanvas();
    c5->SetLogy();
    h1d_signal->Draw();
    h1d_ccbar_pythia_scaled_central->Draw("same hist");
    h1d_ccbar_pythia_scaled_lower->Draw("same hist");
    h1d_ccbar_pythia_scaled_upper->Draw("same hist");

}