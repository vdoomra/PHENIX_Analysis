const int npts = 26;
double x_lohi[npts] = {-999};
double x[npts] = {-999};
double y[npts] = {-999};
double stat_err[npts] = {-999};
double sys_err[npts] = {-999};

const int nbins_mass = 26;
const double mass_range[nbins_mass+1] = {0.87, 0.95, 1.05, 1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.05, 3.10, 3.15, 3.20, 3.30, 3.45, 3.60, 3.70, 3.80, 4.00, 4.25, 4.5};

void plotting_BG_with_uncertainties(){

    TFile* fbg = new TFile("getting_final_background_BG.root","READ");
    TH1D* h1d_mlp_dataBG_temp = (TH1D*)fbg->Get("h1d_mlp_dataBG");
    h1d_mlp_dataBG_temp->SetName("h1d_mlp_dataBG_temp");

    int counter = 0;

    for(int ibin=1; ibin < (int)h1d_mlp_dataBG_temp->GetNbinsX()+1; ibin++){

        x_lohi[counter] = h1d_mlp_dataBG_temp->GetBinWidth(ibin)/2;
        x[counter] = h1d_mlp_dataBG_temp->GetBinCenter(ibin);
        y[counter] = h1d_mlp_dataBG_temp->GetBinContent(ibin);

        sys_err[counter] = 0.05212*h1d_mlp_dataBG_temp->GetBinContent(ibin)*h1d_mlp_dataBG_temp->GetBinWidth(ibin);
        stat_err[counter] = sqrt( pow(h1d_mlp_dataBG_temp->GetBinError(ibin)*h1d_mlp_dataBG_temp->GetBinWidth(ibin), 2) - pow(sys_err[counter], 2) )/h1d_mlp_dataBG_temp->GetBinWidth(ibin);

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
    
    mg->Add(graph_stat,"AP");
    mg->Add(graph_sys, "e2same");
    mg->SetName("mg");

    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
    mg->Draw("a");

    mg->GetXaxis()->SetTitle("m_{ee}[GeV/c^{2}]");
    mg->GetYaxis()->SetTitle("#frac{1}{N_{events}}#frac{#epsilon_{BBC}#sigma_{pp}}{#epsilon_{bias}}#frac{dN}{dm_{ee}} [GeV/c^{2}]^{-1} [In PHENIX Acceptance]");

    TH1D* h1d_mlp_dataBG = new TH1D("h1d_mlp_dataBG","DCA Distribution", nbins_mass, mass_range);
    h1d_mlp_dataBG->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        h1d_mlp_dataBG->SetBinContent(ibin+1, y[ibin]);
        h1d_mlp_dataBG->SetBinError(ibin+1, stat_err[ibin]);

    }

    TH1D* h1d_mlp_dataBG_sys = new TH1D("h1d_mlp_dataBG_sys","DCA Distribution", nbins_mass, mass_range);
    h1d_mlp_dataBG_sys->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        h1d_mlp_dataBG_sys->SetBinContent(ibin+1, y[ibin]);
        h1d_mlp_dataBG_sys->SetBinError(ibin+1, sys_err[ibin]);

    }


    TFile* fout = new TFile("plotting_BG_with_uncertainties.root","RECREATE");
    fout->cd();
    graph_stat->Write();
    graph_sys->Write();
    h1d_mlp_dataBG->Write();
    h1d_mlp_dataBG_sys->Write();
    mg->Write();
    fout->Close();


    
}