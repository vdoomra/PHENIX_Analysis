#include "TFile.h"
#include "TH1.h"

void create_likesign_txt(){

    TFile* f1 = new TFile("pair_pt.root","READ");
    TH1D* h1d_pair_pt_FG = (TH1D*)f1->Get("h1d_signal");

    std::ofstream datafile;
    datafile.open("pair_pt_data.txt",std::ofstream::app);

    int nbinsX = h1d_pair_pt_FG->GetNbinsX();
    for(int j = 1; j < nbinsX+1; j++){

        if(h1d_pair_pt_FG->GetBinContent(j) == 0) continue;
        datafile << h1d_pair_pt_FG->GetBinCenter(j) << "\t" << h1d_pair_pt_FG->GetBinContent(j) << "\t" << h1d_pair_pt_FG->GetBinError(j) << endl;

    }


}