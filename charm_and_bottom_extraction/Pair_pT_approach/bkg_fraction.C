const double mass_lo = 0.87;
const double mass_hi = 4.50;
const double nevents = 4.216682807*pow(10,11)*(0.79);  

const int nfiles = 3;
const char* fileNames[nfiles] = { "phi_ee_2k", "jpsi_2k", "psip_2k"};
const double num[nfiles] = { 2000., 2000., 2000.,};
const double BR[nfiles] = { 0.000297, 0.05971, 7.93e-3};
const int color[nfiles] = { kRed, kGreen-7, kViolet,};
const int line_style[nfiles] = {2, 2, 2};

void bkg_fraction(){
    
    TFile* f1 = new TFile("/Users/vassu/Desktop/new_files/data_files/h1d_ee_FG_unlike.root","READ");
    TH1D* h1d_data = (TH1D*)f1->Get("h1d_ee_FG_cent0");
    h1d_data->SetName("h1d_data");

    for(int ifile = 0; ifile < nfiles; ifile++){

        fin = new TFile(Form("/Users/vassu/Desktop/new_files/data_files/dca_output_%s_gencut_with_correction.root",fileNames[ifile]), "READ");
        hist_temp[ifile] = (TH1D*)fin->Get("h1d_ee_FG_unlike");
        hist_temp[ifile]->SetName(Form("h1d_ee_%s",fileNames[ifile]));

        hist_temp[ifile]->Scale((42.*BR[ifile])/(25000.*num[ifile]));
        hist_temp[ifile]->SetLineColor(color[ifile]);
        hist_temp[ifile]->SetLineStyle(line_style[ifile]);

        hist[ifile] = new TH1D(Form("hist_%s", fileNames[ifile]) ,"hist", nbins_mass, mass_range);
        hist[ifile]->Sumw2();

        for(int ibin=0; ibin < nbins_mass; ibin++){

            int bin_lo = hist_temp[ifile]->FindBin(mass_range[ibin] + 0.001);
            int bin_hi = hist_temp[ifile]->FindBin(mass_range[ibin+1] - 0.001);

            double bin_content = 0; double bin_error = 0;

            for(int i = bin_lo; i <= bin_hi; i++){

                bin_content += hist_temp[ifile]->GetBinContent(i)*hist_temp[ifile]->GetBinWidth(i);
                bin_error += pow(hist_temp[ifile]->GetBinError(i)*hist_temp[ifile]->GetBinWidth(i), 2);

            }

            hist[ifile]->SetBinContent(ibin+1, bin_content/hist[ifile]->GetBinWidth(ibin+1));
            hist[ifile]->SetBinError(ibin+1, 0.0);


        }

    }


    TFile* fccbar = new TFile("/Users/vassu/Desktopnormalizing_two_TGraphs_ccbar_stat.root","READ");
    TH1D* h1d_ccbar_central = (TH1D*)fccbar->Get("h1d_ccbar_central");
    TH1D* h1d_ccbar_lower = (TH1D*)fccbar->Get("h1d_ccbar_lower");
    TH1D* h1d_ccbar_upper = (TH1D*)fccbar->Get("h1d_ccbar_upper");

    TH1D* h1d_ccbar = new TH1D("h1d_ccbar","Mass Distribution", nbins_mass, mass_range);
    h1d_ccbar->Sumw2();

    for(int ibin = 0; ibin < nbins_mass; ibin++){
        
        double bin_center = (mass_range[ibin] + mass_range[ibin+1])/2.;
        h1d_ccbar->SetBinContent(ibin+1, h1d_ccbar_central->GetBinContent(h1d_ccbar_central->FindBin(bin_center)));
        h1d_ccbar->SetBinError(ibin+1, h1d_ccbar_central->GetBinContent(h1d_ccbar_central->FindBin(bin_center)) - h1d_ccbar_lower->GetBinContent(h1d_ccbar_lower->FindBin(bin_center)));

    }


    TFile* fbbbar = new TFile("normalizing_two_TGraphs_bbbar_stat.root","READ");
    TH1D* h1d_bbbar_central = (TH1D*)fbbbar->Get("h1d_bbbar_central");
    TH1D* h1d_bbbar_lower = (TH1D*)fbbbar->Get("h1d_bbbar_lower");
    TH1D* h1d_bbbar_upper = (TH1D*)fbbbar->Get("h1d_bbbar_upper");

    TH1D* h1d_bbbar = new TH1D("h1d_bbbar","Mass Distribution", nbins_mass, mass_range);
    h1d_bbbar->Sumw2();

    for(int ibin = 0; ibin < nbins_mass; ibin++){
        
        double bin_center = (mass_range[ibin] + mass_range[ibin+1])/2.;
        h1d_bbbar->SetBinContent(ibin+1, h1d_bbbar_central->GetBinContent(h1d_bbbar_central->FindBin(bin_center)));
        h1d_bbbar->SetBinError(ibin+1, h1d_bbbar_central->GetBinContent(h1d_bbbar_central->FindBin(bin_center)) - h1d_bbbar_lower->GetBinContent(h1d_bbbar_lower->FindBin(bin_center)));

    }
    

    TFile* fbg = new TFile("plotting_BG_with_uncertainties.root","READ");
    TH1D* h1d_mlp_dataBG = (TH1D*)fbg->Get("h1d_mlp_dataBG");

    h1d_mlp_dataBG->Draw("samep");

}

