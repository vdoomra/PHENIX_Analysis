
const int numTemplates = 2;

double initial_fraction[numTemplates] = { 0.5, 0.5};
double min_fraction[numTemplates] = { 0.0, 0.0};
double max_fraction[numTemplates] = { 1.0, 1.0};

const double step_size = 0.001;

TObjArray* templates = NULL;
TH1D* h1d_dataFG = NULL;
TH1D* h1d_ccbar = NULL;
TH1D* h1d_bbbar = NULL;

std::vector<double> hData;
std::vector<double> parameters_bbbar;
std::vector<double> parameters_ccbar;

const int nbins_mass = 6;
const double mass_range[nbins_mass+1] = {0.0, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0 };

const double nevents = 4.216682807*pow(10,11)*(0.79);

double NegativeLogLikelihood(const double* params) {

    double logLikelihood = 0.0;
    
    for (int i = 1; i <= nbins_mass; ++i) {

        double observed = hData[i-1];

        double expected = 0;
                
        for (int j = 0; j < numTemplates; ++j) {

            TH1D* templateHist = (TH1D*)templates->At(j);
            expected += params[j] * templateHist->GetBinContent(i);
        }

        if (expected > 0) logLikelihood -= observed * log(expected) - expected;
    }

    return logLikelihood;
}

void template_fits_pair_pt(){

    TFile* f2 = new TFile("pair_pt.root","READ");
    TH1D* h1d_pair_pt_data_original = (TH1D*)f2->Get("h1d_signal");
    h1d_pair_pt_data_original->Scale(1./h1d_pair_pt_data_original->Integral("width"));

    TH1D* h1d_pair_pt_data_temp = new TH1D("h1d_pair_pt_data_temp", "Pair pT Distribution", 200, 0.0, 5.0);

    std::ifstream inFile_data("mlp_pair_pt_data.txt");
    double x1, y1, err_y1;

    while(!inFile_data.eof()){

        inFile_data >> x1 >> y1 >> err_y1;

        int binX = h1d_pair_pt_data_temp->FindBin(x1);
        h1d_pair_pt_data_temp->SetBinContent(binX, y1);
        h1d_pair_pt_data_temp->SetBinError(binX, err_y1);

    }

    TH1D* h1d_pair_pt_data = new TH1D("h1d_pair_pt_data","h1d_pair_pt_data", nbins_mass, mass_range);
    h1d_pair_pt_data->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        double bin_center = (mass_range[ibin] + mass_range[ibin+1])/2.;

        h1d_pair_pt_data->SetBinContent(ibin+1, h1d_pair_pt_data_temp->GetBinContent(h1d_pair_pt_data_temp->FindBin(bin_center)));
        h1d_pair_pt_data->SetBinError(ibin+1, h1d_pair_pt_data_temp->GetBinError(h1d_pair_pt_data_temp->FindBin(bin_center)));

    }

    h1d_pair_pt_data->Scale(1./h1d_pair_pt_data->Integral("width"));

    TH1D* h1d_pair_pt_ccbar_temp = new TH1D("h1d_pair_pt_ccbar_temp", "Pair pT Distribution", 200, 0.0, 5.0);

    std::ifstream inFile_ccbar("mlp_pair_pt_ccbar.txt");

    while(!inFile_ccbar.eof()){

        inFile_ccbar >> x1 >> y1 >> err_y1;

        int binX = h1d_pair_pt_ccbar_temp->FindBin(x1);
        h1d_pair_pt_ccbar_temp->SetBinContent(binX, y1);
        h1d_pair_pt_ccbar_temp->SetBinError(binX, err_y1);

    }

    h1d_pair_pt_ccbar_temp->Scale(10.);

    TH1D* h1d_pair_pt_ccbar = new TH1D("h1d_pair_pt_ccbar", "Pair pT Distribution", nbins_mass, mass_range);
    h1d_pair_pt_ccbar->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_pair_pt_ccbar_temp->FindBin(mass_range[ibin] + 0.01);
        int bin_hi = h1d_pair_pt_ccbar_temp->FindBin(mass_range[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_pair_pt_ccbar_temp->GetBinContent(i)*h1d_pair_pt_ccbar_temp->GetBinWidth(i);
            bin_error += pow(h1d_pair_pt_ccbar_temp->GetBinError(i)*h1d_pair_pt_ccbar_temp->GetBinWidth(i), 2);

        }

        h1d_pair_pt_ccbar->SetBinContent(ibin+1, bin_content/h1d_pair_pt_ccbar->GetBinWidth(ibin+1));
        h1d_pair_pt_ccbar->SetBinError(ibin+1, 0.0);


    }

    h1d_pair_pt_ccbar->Scale(1./h1d_pair_pt_ccbar->Integral("width"));

    TH1D* h1d_pair_pt_bbbar_temp = new TH1D("h1d_pair_pt_bbbar_temp", "Pair pT Distribution", 200, 0.0, 5.0);

    std::ifstream inFile_bbbar("mlp_pair_pt_bbbar.txt");

    while(!inFile_bbbar.eof()){

        inFile_bbbar >> x1 >> y1 >> err_y1;

        int binX = h1d_pair_pt_bbbar_temp->FindBin(x1);
        h1d_pair_pt_bbbar_temp->SetBinContent(binX, y1);
        h1d_pair_pt_bbbar_temp->SetBinError(binX, err_y1);

    }

    h1d_pair_pt_bbbar_temp->Scale(10.);

    TH1D* h1d_pair_pt_bbbar = new TH1D("h1d_pair_pt_bbbar", "Pair pT Distribution", nbins_mass, mass_range);
    h1d_pair_pt_bbbar->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_pair_pt_bbbar_temp->FindBin(mass_range[ibin] + 0.01);
        int bin_hi = h1d_pair_pt_bbbar_temp->FindBin(mass_range[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_pair_pt_bbbar_temp->GetBinContent(i)*h1d_pair_pt_bbbar_temp->GetBinWidth(i);
            bin_error += pow(h1d_pair_pt_bbbar_temp->GetBinError(i)*h1d_pair_pt_bbbar_temp->GetBinWidth(i), 2);

        }

        h1d_pair_pt_bbbar->SetBinContent(ibin+1, bin_content/h1d_pair_pt_bbbar->GetBinWidth(ibin+1));
        h1d_pair_pt_bbbar->SetBinError(ibin+1, 0.0);


    }

    h1d_pair_pt_bbbar->Scale(1./h1d_pair_pt_bbbar->Integral("width"));


    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
    c1->cd();

    h1d_pair_pt_data->Draw();
    h1d_pair_pt_ccbar->Draw("same");
    h1d_pair_pt_bbbar->Draw("same");

    templates = new TObjArray();
    templates->Add(h1d_pair_pt_ccbar);
    templates->Add(h1d_pair_pt_bbbar);

    parameters_bbbar.clear();
    parameters_ccbar.clear();
        
    for(int itr = 0; itr < 1e5; itr++){

        if(itr%50000==0) cout << "Iterations " << itr << " Completed" << endl;

        hData.clear();

        int nbins = h1d_pair_pt_data->GetNbinsX();

        for (int ii = 1; ii <= nbins; ++ii) { hData.push_back(gRandom->Gaus(h1d_pair_pt_data->GetBinContent(ii), h1d_pair_pt_data->GetBinError(ii))); }

        ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        min->SetTolerance(1e-6);

        ROOT::Math::Functor functor;

        functor = ROOT::Math::Functor(NegativeLogLikelihood, numTemplates);

        min->SetFunction(functor);

        for (int ii = 0; ii < numTemplates; ++ii) {

            min->SetLimitedVariable(ii, Form("f%d", ii + 1), initial_fraction[ii], step_size, min_fraction[ii], max_fraction[ii]); // Initial guess and step size
        }

        min->Minimize();
        const double* fractions = min->X();

        bool fit_is_bad = false;
            for (int j = 0; j < numTemplates; ++j) {
                if (fractions[j] < 0.05) {
                    fit_is_bad = true;
                    break;
                }
            }

        if (fit_is_bad) continue;

        parameters_ccbar.push_back(fractions[0]);
        parameters_bbbar.push_back(fractions[1]);
        

    } 


    TH1D* h1d_para_bbbar = new TH1D("h1d_para_bbbar","Parameters Distribution",100,0,1);
    for(int ii = 0; ii < parameters_bbbar.size(); ii++) h1d_para_bbbar->Fill(parameters_bbbar[ii]);

    TH1D* h1d_para_ccbar = new TH1D("h1d_para_ccbar","Parameters Distribution",100,0,1);
    for(int ii = 0; ii < parameters_ccbar.size(); ii++) h1d_para_ccbar->Fill(parameters_ccbar[ii]);

    TCanvas* c2 = new TCanvas();
    c2->cd();
    h1d_para_ccbar->Draw();

    TCanvas* c3 = new TCanvas();
    c3->cd();
    h1d_para_bbbar->Draw();

    TCanvas* c5 = new TCanvas();
    c5->SetLogy();
    h1d_pair_pt_data->Draw();

    h1d_pair_pt_ccbar->Scale(h1d_para_ccbar->GetMean());
    h1d_pair_pt_bbbar->Scale(h1d_para_bbbar->GetMean());

    h1d_pair_pt_ccbar->Draw("same");
    h1d_pair_pt_bbbar->Draw("same");

    TH1D* hsum = (TH1D*)h1d_pair_pt_ccbar->Clone();
    hsum->SetName("hsum");
    hsum->Add(h1d_pair_pt_bbbar);
    hsum->Draw("same");

    TFile* fout = new TFile("template_fits_pair_pt_MLP_data.root","RECREATE");
    fout->cd();
    h1d_pair_pt_data_original->Write();
    h1d_pair_pt_data->Write();
    h1d_pair_pt_ccbar->Write();
    h1d_pair_pt_bbbar->Write();
    h1d_para_ccbar->Write();
    h1d_para_bbbar->Write();
    hsum->Write();
    fout->Close();

}