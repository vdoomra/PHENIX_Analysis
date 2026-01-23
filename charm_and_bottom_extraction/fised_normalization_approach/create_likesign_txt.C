#include "TFile.h"
#include "TH1.h"

const int nbins_dca = 20;
double dca_range[nbins_dca+1] = {0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0};

void create_likesign_txt(){

    TFile* fdata_FG = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_drell_yan_10k_gencut_with_correction.root","READ");
    TH2D* h2d_ee_FG_bbbar_central = (TH2D*)fdata_FG->Get("h2d_ee_FG_unlike");
    h2d_ee_FG_bbbar_central->GetXaxis()->SetRangeUser(0.87, 4.50);

    TH1D* h1d_pair_pt_bbbar_central_temp = (TH1D*)h2d_ee_FG_bbbar_central->ProjectionY();
    h1d_pair_pt_bbbar_central_temp->SetName("h1d_pair_pt_bbbar_central_temp");

    TH1D* h1d_pair_pt_bbbar_central = new TH1D("h1d_pair_pt_bbbar_central","Mass Distribution", nbins_dca, dca_range);
    h1d_pair_pt_bbbar_central->Sumw2();

    for(int ibin=0; ibin < nbins_dca; ibin++){

        int bin_lo = h1d_pair_pt_bbbar_central_temp->FindBin(dca_range[ibin] + 0.01);
        int bin_hi = h1d_pair_pt_bbbar_central_temp->FindBin(dca_range[ibin+1] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_pair_pt_bbbar_central_temp->GetBinContent(i);
            bin_error += pow(h1d_pair_pt_bbbar_central_temp->GetBinError(i), 2);

        }

        h1d_pair_pt_bbbar_central->SetBinContent(ibin+1, bin_content);
        h1d_pair_pt_bbbar_central->SetBinError(ibin+1, sqrt(bin_error));


    }

    h1d_pair_pt_bbbar_central->Scale(h1d_pair_pt_bbbar_central_temp->Integral("width")/h1d_pair_pt_bbbar_central->Integral("width"));

    TCanvas* c1 = new TCanvas();
    c1->cd();
    c1->SetLogy();
    h1d_pair_pt_bbbar_central->Draw();
    h1d_pair_pt_bbbar_central_temp->Draw("same");

    std::ofstream datafile;
    datafile.open("pair_pt_drellyan.txt",std::ofstream::app);

    int nbinsX = h1d_pair_pt_bbbar_central->GetNbinsX();
    for(int j = 1; j < nbinsX+1; j++){

        if(h1d_pair_pt_bbbar_central->GetBinContent(j) == 0) continue;
        datafile << h1d_pair_pt_bbbar_central->GetBinCenter(j) << "\t" << h1d_pair_pt_bbbar_central->GetBinContent(j) << "\t" << h1d_pair_pt_bbbar_central->GetBinError(j) << endl;

    }


}