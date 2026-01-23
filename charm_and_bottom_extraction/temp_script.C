const int numTemplates = 4;
int nbins = 9;

double initial_fraction[numTemplates] = {0.5, 0.5, 0.5, 0.5}; 
double min_fraction[numTemplates] = {0.0, 0.0, 0.0, 0.0};
double max_fraction[numTemplates] = {1.0, 1.0, 1.0, 1.0};

const int mass_regions = 13;
double mass_lo[mass_regions] = { 0.70, 0.87, 1.30, 1.75, 2.20, 2.40, 1.60, 2.40, 3.30, 2.60, 3.10, 3.30, 1.15};
double mass_hi[mass_regions] = { 1.15, 1.30, 2.05, 2.40, 2.80, 2.95, 2.80, 3.45, 4.50, 4.00, 4.25, 5.50, 1.75};

const char* mass_lo_name[mass_regions] = { "070", "087", "130", "175", "220", "240", "160", "240", "330", "260", "310", "330", "115"}; 
const char* mass_hi_name[mass_regions] = { "115", "130", "205", "240", "280", "295", "280", "345", "450", "400", "425", "550", "175"};

const char* region_name[mass_regions] = {"Region1", "Region2", "Region3", "Region4", "Region5", "Region6", "Region7", "Region8", "Region9", "Region10", "Region11", "Region12", "Region13"};
int templates_used[mass_regions] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
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


void temp_script(){

    fout = new TFile("template_fits.root","RECREATE");

    fin = new TFile("/Users/vassu/Desktop/new_files/data_files/producing_templates.root","READ");

    // These histograms will be used to calculate the initial fractions
    
    TH1D* h1d_ee_FG_unlike = (TH1D*)fin->Get("h1d_ee_FG_unlike_original");
    h1d_ee_FG_unlike->SetName("h1d_ee_FG_unlike");

    TH1D* h1d_ee_FG_jpsi = (TH1D*)fin->Get("h1d_ee_FG_jpsi");
    h1d_ee_FG_jpsi->SetName("h1d_ee_FG_jpsi");

    TH1D* h1d_ee_FG_phi = (TH1D*)fin->Get("h1d_ee_FG_phi");
    h1d_ee_FG_phi->SetName("h1d_ee_FG_phi");

    TH1D* h1d_ee_FG_omega = (TH1D*)fin->Get("h1d_ee_FG_omega");
    h1d_ee_FG_omega->SetName("h1d_ee_FG_omega");

    TH1D* h1d_ee_FG_psip = (TH1D*)fin->Get("h1d_ee_FG_psip");
    h1d_ee_FG_psip->SetName("h1d_ee_FG_psip");


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

        TH1D* h1d_phi = NULL;
        TH1D* h1d_jpsi = NULL;
        TH1D* h1d_psip = NULL;

        if(i < 2){

            TString hname_phi = Form("h1d_phi_%s_%s", mass_lo_name[i], mass_hi_name[i]);
            h1d_phi = (TH1D*)fin->Get(hname_phi);

        }

        else if( (i > 1 && i < 8) || i == mass_regions-4 || i == mass_regions-1){

            TString hname_jpsi = Form("h1d_jpsi_%s_%s", mass_lo_name[i], mass_hi_name[i]);
            h1d_jpsi = (TH1D*)fin->Get(hname_jpsi);

        }

        else {

            TString hname_psip = Form("h1d_psip_%s_%s", mass_lo_name[i], mass_hi_name[i]);
            h1d_psip = (TH1D*)fin->Get(hname_psip);

        }

        for(int ii = 0; ii < numTemplates; ii++){

            if(ii == 1) continue; //Prompt component

            else{

                initial_fraction[ii] = 0.5;
                min_fraction[ii] = initial_fraction[ii] - 0.5;
                max_fraction[ii] = initial_fraction[ii] + 0.5;

            }

        }

        if(i==0){

            double total = (h1d_ee_FG_phi->Integral(h1d_ee_FG_phi->FindBin(mass_lo[i]+0.01), h1d_ee_FG_phi->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width")) + (h1d_ee_FG_omega->Integral(h1d_ee_FG_omega->FindBin(mass_lo[i]+0.01), h1d_ee_FG_omega->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width"));

            initial_fraction[1] = total;
            min_fraction[1] = initial_fraction[1] - 0.05;
            if(min_fraction[1] < 0) min_fraction[1] = 0.0;
            max_fraction[1] = initial_fraction[1] + 0.05;


        }

        else if(i==1){

            initial_fraction[1] = h1d_ee_FG_phi->Integral(h1d_ee_FG_phi->FindBin(mass_lo[i]+0.01), h1d_ee_FG_phi->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width");
            min_fraction[1] = initial_fraction[1] - 0.05;
            if(min_fraction[1] < 0) min_fraction[1] = 0.0;
            max_fraction[1] = initial_fraction[1] + 0.05;

        }

        else if(i > 1 && i <8 || i == mass_regions-1){
            
            initial_fraction[1] = h1d_ee_FG_jpsi->Integral(h1d_ee_FG_jpsi->FindBin(mass_lo[i]+0.01), h1d_ee_FG_jpsi->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width");
            min_fraction[1] = initial_fraction[1] - 0.05;
            if(min_fraction[1] < 0) min_fraction[1] = 0.0;
            max_fraction[1] = initial_fraction[1] + 0.05;

        }

        else if(i==mass_regions-3 || i == mass_regions-4){

            double total = h1d_ee_FG_jpsi->Integral(h1d_ee_FG_jpsi->FindBin(mass_lo[i]+0.01), h1d_ee_FG_jpsi->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width") + h1d_ee_FG_psip->Integral(h1d_ee_FG_psip->FindBin(mass_lo[i]+0.01), h1d_ee_FG_psip->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width");
            
            initial_fraction[1] = total;
            min_fraction[1] = initial_fraction[1] - 0.05;
            if(min_fraction[1] < 0) min_fraction[1] = 0.0;
            max_fraction[1] = initial_fraction[1] + 0.05;

        }

        else{
            
            initial_fraction[1] = h1d_ee_FG_psip->Integral(h1d_ee_FG_psip->FindBin(mass_lo[i]+0.01), h1d_ee_FG_psip->FindBin(mass_hi[i]-0.01),"width")/h1d_ee_FG_unlike->Integral(h1d_ee_FG_unlike->FindBin(mass_lo[i]+0.01), h1d_ee_FG_unlike->FindBin(mass_hi[i]-0.01),"width");
            min_fraction[1] = initial_fraction[1] - 0.05;
            if(min_fraction[1] < 0) min_fraction[1] = 0.0;
            max_fraction[1] = initial_fraction[1] + 0.05;

        }

        templates = new TObjArray();

        if(i < 2) h1d_prompt = (TH1D*)h1d_phi->Clone();
        else if((i > 1 && i < 8) || i == mass_regions-4 || i == mass_regions-1) h1d_prompt = (TH1D*)h1d_jpsi->Clone();
        else h1d_prompt = (TH1D*)h1d_psip->Clone();

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

        TH1D* h1d_para_ccbar = new TH1D(Form("h1d_para_ccbar_%s", region_name[i]),"Parameters Distribution",100,0.0,1);
        for(int ii = 0; ii < parameters_ccbar.size(); ii++) h1d_para_ccbar->Fill(parameters_ccbar[ii]);

        TH1D* h1d_para_bbbar = new TH1D(Form("h1d_para_bbbar_%s", region_name[i]),"Parameters Distribution",100,0.0,1);
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
      

    }

    fout->Close();


}