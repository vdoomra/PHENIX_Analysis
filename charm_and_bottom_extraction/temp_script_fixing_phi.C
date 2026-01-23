const int numTemplates = 4;
int nbins = 9;

double initial_fraction[numTemplates] = {0.5, 0.5, 0.5, 0.5}; 
double min_fraction[numTemplates] = {0.0, 0.0, 0.0, 0.0};
double max_fraction[numTemplates] = {1.0, 1.0, 1.0, 1.0};

const int mass_regions = 1;
double mass_lo[mass_regions] = {0.87};
double mass_hi[mass_regions] = {1.30};

const char* mass_lo_name[mass_regions] = {"087"}; 
const char* mass_hi_name[mass_regions] = {"130"};

const char* region_name[mass_regions] = {"Region1"};
int templates_used[mass_regions] = {4};
const double step_size = 0.0001;

TH1D* hsum = NULL;

TObjArray* templates = NULL;
TH1D* h1d_prompt = NULL;

TFile* fin = NULL;
TFile* fout = NULL;

std::vector<double> hData;
std::vector<double> parameters_BG;
std::vector<double> parameters_prompt;
std::vector<double> parameters_bbbar;
std::vector<double> parameters_ccbar;

double NegativeLogLikelihood(const double* params) {

    double logLikelihood = 0.0;
    
    for (int i = 1; i <= nbins; ++i) {

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


void temp_script_fixing_phi(){

    fout = new TFile("template_fits_fixing_phi.root","RECREATE");

    fin = new TFile("data_files/producing_templates.root","READ");

    // These histograms will be used to calculate the initial fractions
    
    TH1D* h1d_ee_FG_unlike = (TH1D*)fin->Get("h1d_ee_FG_unlike_original");
    h1d_ee_FG_unlike->SetName("h1d_ee_FG_unlike");

    TH1D* h1d_ee_FG_phi = (TH1D*)fin->Get("h1d_ee_FG_phi");
    h1d_ee_FG_phi->SetName("h1d_ee_FG_phi");

    TFile* fbbbar = new TFile("/Users/vassu/Desktop/new_files/likesign_analysis_alternate_approach.root","READ");
    TH1D* h1d_bbbar_dataBinning_central = (TH1D*)fbbbar->Get("h1d_bbbar_pythia_unlike_dataBinning_central");
    TH1D* h1d_bbbar_dataBinning_lower = (TH1D*)fbbbar->Get("h1d_bbbar_pythia_unlike_dataBinning_lower");
    TH1D* h1d_bbbar_dataBinning_upper = (TH1D*)fbbbar->Get("h1d_bbbar_pythia_unlike_dataBinning_upper");


    for(int i = 0; i < mass_regions; i++){

        cout << "Starting process in Mass Range: " << mass_lo[i] << " - " << mass_hi[i] << endl;

        TString hname_dataFG = Form("h1d_FG_%s_%s", mass_lo_name[i], mass_hi_name[i]);
        TH1D* h1d_dataFG = (TH1D*)fin->Get(hname_dataFG);

        TString hname_dataFG_original = Form("h1d_FG_%s_%s_original", mass_lo_name[i], mass_hi_name[i]);
        TH1D* h1d_dataFG_original = (TH1D*)fin->Get(hname_dataFG_original);

        TString hname_dataBG = Form("h1d_BG_%s_%s", mass_lo_name[i], mass_hi_name[i]);
        TH1D* h1d_dataBG = (TH1D*)fin->Get(hname_dataBG);

        TString hname_bbbar = Form("h1d_bbbar_%s_%s", mass_lo_name[i], mass_hi_name[i]);
        TH1D* h1d_bbbar = (TH1D*)fin->Get(hname_bbbar);

        TString hname_ccbar = Form("h1d_ccbar_%s_%s", mass_lo_name[i], mass_hi_name[i]);
        TH1D* h1d_ccbar = (TH1D*)fin->Get(hname_ccbar);

        TString hname_phi = Form("h1d_phi_%s_%s", mass_lo_name[i], mass_hi_name[i]);
        TH1D* h1d_phi = (TH1D*)fin->Get(hname_phi);


        for(int ii = 0; ii < numTemplates; ii++){

            if(ii == 3){

                initial_fraction[ii] = h1d_bbbar_dataBinning_central->Integral(h1d_bbbar_dataBinning_central->FindBin(mass_lo[i]+0.01), h1d_bbbar_dataBinning_central->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width");
                min_fraction[ii] = h1d_bbbar_dataBinning_lower->Integral(h1d_bbbar_dataBinning_lower->FindBin(mass_lo[i]+0.01), h1d_bbbar_dataBinning_lower->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width");
                max_fraction[ii] = h1d_bbbar_dataBinning_upper->Integral(h1d_bbbar_dataBinning_upper->FindBin(mass_lo[i]+0.01), h1d_bbbar_dataBinning_upper->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width");

            }

            else{

                initial_fraction[ii] = 0.5;
                min_fraction[ii] = initial_fraction[ii] - 0.5;
                max_fraction[ii] = initial_fraction[ii] + 0.5;

            }

        }

        templates = new TObjArray();

        h1d_prompt = (TH1D*)h1d_phi->Clone();

        for(int cc=0; cc < numTemplates; cc++){  cout << initial_fraction[cc] << "\t" << min_fraction[cc] << "\t" << max_fraction[cc]  << endl;  }

        templates->Add(h1d_dataBG);
        templates->Add(h1d_prompt);
        templates->Add(h1d_ccbar);
        templates->Add(h1d_bbbar);
    
        parameters_BG.clear();
        parameters_prompt.clear();
        parameters_bbbar.clear();
        parameters_ccbar.clear();
        
        for(int itr = 0; itr < 1e5; itr++){

            if(itr%50000==0) cout << "Iterations " << itr << " Completed" << endl;

            hData.clear();

            int num = h1d_dataFG->GetNbinsX();
    
            for (int ii = 1; ii <= num; ++ii) { 
                hData.push_back(gRandom->Gaus(h1d_dataFG->GetBinContent(ii), h1d_dataFG->GetBinError(ii))); 
            }

            ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            min->SetTolerance(1e-5);

            ROOT::Math::Functor functor;

            functor = ROOT::Math::Functor(NegativeLogLikelihood, templates_used[i]);

            min->SetFunction(functor);

            for (int ii = 0; ii < templates_used[i]; ++ii) {

                min->SetLimitedVariable(ii, Form("f%d", ii + 1), initial_fraction[ii], step_size, min_fraction[ii], max_fraction[ii]); 
            }
    
            min->Minimize();
            const double* fractions = min->X();

            bool fit_is_bad = false;
            for (int j = 0; j < templates_used[i]; ++j) {
                if (fractions[j] < 0.005) {
                    fit_is_bad = true;
                    break;
                }
            }

            if (fit_is_bad) continue;

            parameters_BG.push_back(fractions[0]);
            parameters_prompt.push_back(fractions[1]);
            parameters_ccbar.push_back(fractions[2]);
            parameters_bbbar.push_back(fractions[3]);
            
           
        } 

        TH1D* h1d_para_BG = new TH1D(Form("h1d_para_BG_%s", region_name[i]),"Parameters Distribution",100,0.0,1);
        for(int ii = 0; ii < parameters_BG.size(); ii++) h1d_para_BG->Fill(parameters_BG[ii]);

        TH1D* h1d_para_prompt = new TH1D(Form("h1d_para_prompt_%s", region_name[i]),"Parameters Distribution",100,0.0,1);
        for(int ii = 0; ii < parameters_prompt.size(); ii++) h1d_para_prompt->Fill(parameters_prompt[ii]);

        TH1D* h1d_para_ccbar = new TH1D(Form("h1d_para_ccbar_%s", region_name[i]),"Parameters Distribution",100,0.0,1.0);
        for(int ii = 0; ii < parameters_ccbar.size(); ii++) h1d_para_ccbar->Fill(parameters_ccbar[ii]);

        TH1D* h1d_para_bbbar = new TH1D(Form("h1d_para_bbbar_%s", region_name[i]),"Parameters Distribution",100,0.0,1.0);
        for(int ii = 0; ii < parameters_bbbar.size(); ii++) h1d_para_bbbar->Fill(parameters_bbbar[ii]);

        fout->cd();
        h1d_dataFG->Write(Form("DCA_dataFG_%s",region_name[i]));
        h1d_dataFG_original->Write(Form("DCA_dataFG_%s_original",region_name[i]));
        h1d_dataBG->Scale(h1d_para_BG->GetMean());
        h1d_dataBG->Write(Form("DCA_dataBG_%s",region_name[i]));
        h1d_prompt->Scale(h1d_para_prompt->GetMean());
        h1d_prompt->Write(Form("DCA_prompt_%s",region_name[i]));
        h1d_ccbar->Scale(h1d_para_ccbar->GetMean());
        h1d_ccbar->Write(Form("DCA_ccbar_%s",region_name[i]));
        h1d_bbbar->Scale(h1d_para_bbbar->GetMean());
        h1d_bbbar->Write(Form("DCA_bbbar_%s",region_name[i]));

        TH1D* hDCA_sum = (TH1D*)h1d_dataBG->Clone();
        hDCA_sum->SetName(Form("hDCA_sum_%s",region_name[i]));
        hDCA_sum->Add(h1d_prompt);
        hDCA_sum->Add(h1d_ccbar);
        hDCA_sum->Add(h1d_bbbar);
        hDCA_sum->Write();

        h1d_para_BG->Write(Form("h1d_para_BG_%s",region_name[i]));
        h1d_para_prompt->Write(Form("h1d_para_prompt_%s",region_name[i]));
        h1d_para_ccbar->Write(Form("h1d_para_ccbar_%s",region_name[i]));
        h1d_para_bbbar->Write(Form("h1d_para_bbbar_%s",region_name[i]));

        TF1* func = new TF1("func","gaus", h1d_para_prompt->GetMean()-1.5*h1d_para_prompt->GetStdDev(), h1d_para_prompt->GetMean()+1.5*h1d_para_prompt->GetStdDev());
        h1d_para_prompt->Fit(func,"R+");

        double initial_value = h1d_ee_FG_phi->Integral(h1d_ee_FG_phi->FindBin(mass_lo[i]+0.01), h1d_ee_FG_phi->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width");

        double scale_factor_central = func->GetParameter(1)/initial_value;

        TH1D* h1d_ee_FG_phi_central = (TH1D*)h1d_ee_FG_phi->Clone();
        h1d_ee_FG_phi_central->SetName("h1d_ee_FG_phi_central");
        h1d_ee_FG_phi_central->Scale(scale_factor_central);

        double scale_factor_lower = (func->GetParameter(1)-func->GetParameter(2))/initial_value;

        TH1D* h1d_ee_FG_phi_lower = (TH1D*)h1d_ee_FG_phi->Clone();
        h1d_ee_FG_phi_lower->SetName("h1d_ee_FG_phi_lower");
        h1d_ee_FG_phi_lower->Scale(scale_factor_lower);

        double scale_factor_upper = (func->GetParameter(1)+func->GetParameter(2))/initial_value;

        TH1D* h1d_ee_FG_phi_upper = (TH1D*)h1d_ee_FG_phi->Clone();
        h1d_ee_FG_phi_upper->SetName("h1d_ee_FG_phi_upper");
        h1d_ee_FG_phi_upper->Scale(scale_factor_upper);

        h1d_ee_FG_phi->Write();
        h1d_ee_FG_phi_central->Write();
        h1d_ee_FG_phi_lower->Write();
        h1d_ee_FG_phi_upper->Write();
      

    }

    fout->Close();


}