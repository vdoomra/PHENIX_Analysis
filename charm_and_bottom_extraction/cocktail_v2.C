#include <TF1.h>
#include <TFile.h>
#include <iostream>
#include <TGraphErrors.h>

const int npts = 40;
double x_lohi[npts] = {-999};
double x[npts] = {-999};
double y[npts] = {-999};
double stat_err[npts] = {-999};
double sys_err[npts] = {-999};

const int nfiles = 2;
const char* fileNames[nfiles] = { "jpsi_2k", "psip_2k"};
const double num[nfiles] = { 2000., 2000.,};
const double BR[nfiles] = { 0.05971, 7.93e-3};
const int color[nfiles] = { kGreen-7, kViolet,};
const int line_style[nfiles] = {1,1};

TH1D* hist_temp[nfiles] = {NULL};
TH1D* hist[nfiles] = {NULL};
TFile* fin = NULL;

const double nevents = 4.216682807*pow(10,11)*(0.79);

const int nbins_mass = 15;
const double mass_range[nbins_mass+1] = {1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.05, 3.10, 3.15, 3.20, 3.30, 3.45};


//These are the systematics on the FG pairs

double combinedFunction(double *x, double *par) {
    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    if(x[0] > 3.50) fERT_val = 0.0;
    double fn0_val = par[3];
    double feoverp_val = par[4] + par[5]*x[0] + par[6]*x[0]*x[0];
    if(x[0] > 3.20) feoverp_val = 0.03262;

    // 2% uncertainty on the number of events

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}

double combinedFunction_jpsi(double *x, double *par) {

    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    double fn0_val = par[3] + par[4] * x[0] + par[5]*x[0]*x[0];
    double feoverp_val = par[6] + par[7]*x[0] + par[8]*x[0]*x[0];

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}

double combinedFunction_psip(double *x, double *par) {
    
    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    double fn0_val = par[3] + par[4] * x[0] + par[5]*x[0]*x[0];
    double feoverp_val = par[6] + par[7]*x[0] + par[8]*x[0]*x[0];

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val);
}

void cocktail_v2(){

    TCanvas *c1 = new TCanvas("c", "Unified Two Panel Plot", 800, 800);

    TPad *pad1 = new TPad("pad1", "Top pad", 0.0, 0.35, 1.0, 1.0);
    TPad *pad2 = new TPad("pad2", "Bottom pad", 0.0, 0.0, 1.0, 0.35);


    pad1->SetBottomMargin(0.0); 
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.35);

    pad1->SetTickx(0);

    pad1->Draw();
    pad2->Draw();


    TFile* f1 = new TFile("/Users/vassu/Desktop/new_files/data_files/h1d_ee_FG_unlike.root","READ");
    TH1D* h1d_data = (TH1D*)f1->Get("h1d_ee_FG_cent0");
    h1d_data->SetName("h1d_data");

    TF1* fcombined = new TF1("fcombined", combinedFunction, 0.0, 5.5, 7);
    fcombined->SetParameters(0.102, -0.0022, -0.00492, 0.06424, 0.06622, 0.08956, -0.03127);

    TF1* fcombined_jpsi = new TF1("fcombined_jpsi", combinedFunction_jpsi, 0.0, 5.5, 9);
    fcombined_jpsi->SetParameters(0.0973, 0.00915, -0.00925, 0.289, -0.1783, 0.0318, 0.412, -0.268, 0.0465);

    TF1* fcombined_psip = new TF1("fcombined_psip", combinedFunction_psip, 0.0, 5.5, 9);
    fcombined_psip->SetParameters(0.091, -0.0108, -0.0025, 0.306, -0.1561, 0.0234, 0.4744, -0.272, 0.0404);

    int nbins = h1d_data->GetNbinsX();
    int counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = h1d_data->GetBinCenter(ibin);
        y[counter] = h1d_data->GetBinContent(ibin);
        x_lohi[counter] = h1d_data->GetBinWidth(ibin)/2;

        stat_err[counter] = h1d_data->GetBinError(ibin);
        sys_err[counter] = fcombined->Eval(h1d_data->GetBinCenter(ibin))*h1d_data->GetBinContent(ibin);

        counter++;


    }

    TGraphErrors* graph_stat = new TGraphErrors(counter, x, y, 0, stat_err);
    graph_stat->SetMarkerStyle(20);
    graph_stat->SetMarkerSize(0.9);
    graph_stat->SetMarkerColor(kBlack);
    graph_stat->SetLineColor(kBlack);

    TGraphErrors* graph_sys = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    graph_sys->SetFillColorAlpha(kBlack, 0.3);
    graph_sys->SetLineColor(kBlack);

    TMultiGraph* mg = new TMultiGraph();
    
    mg->Add(graph_stat,"AP");
    mg->Add(graph_sys, "e2same");
    mg->SetName("mg");

    mg->GetXaxis()->SetRangeUser(1.45, 3.45);

    pad1->cd();
    mg->Draw("a");

    mg->GetXaxis()->SetTitle("m_{ee}[GeV/c^{2}]");
    mg->GetYaxis()->SetTitle("#frac{1}{N_{events}}#frac{#epsilon_{BBC}#sigma_{pp}}{#epsilon_{bias}}#frac{dN}{dm_{ee}} [GeV/c^{2}]^{-1} [In PHENIX Acceptance]");

    TLatex *   tex1 = new TLatex(1.2,1e-4,"Run 2015 p+p, #sqrt{s} = 200 GeV, |#eta| < 0.35");
    tex1->SetTextAlign(13);
    tex1->SetTextSize(0.052);
    tex1->SetLineWidth(1);
    tex1->SetTextFont(42);
    tex1->Draw();

    TLatex *   tex4 = new TLatex(3.5,1e-5, "p_{T}^{pair} < 5 GeV/c");
    tex4->SetTextAlign(13);
    tex4->SetTextSize(0.052);
    tex4->SetLineWidth(1);
    tex4->SetTextFont(42);
    tex4->Draw();

    TLatex *   tex5 = new TLatex(3.5,5e-5, "p_{T}^{trig, e} > 1.5 GeV/c");
    tex5->SetTextAlign(13);
    tex5->SetTextSize(0.052);
    tex5->SetLineWidth(1);
    tex5->SetTextFont(42);
    tex5->Draw();

    for(int ifile = 0; ifile < nfiles; ifile++){

        fin = new TFile(Form("/Users/vassu/Desktop/new_files/data_files/dca_output_%s.root",fileNames[ifile]), "READ");
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
            hist[ifile]->SetBinError(ibin+1, sqrt(bin_error)/hist[ifile]->GetBinWidth(ibin+1));


        }

        hist[ifile]->Draw("same hist");

    }

    TFile* fdrellyan = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_drell_yan_10k_gencut.root","READ");
    TH1D* h1d_drell_yan_temp = (TH1D*)fdrellyan->Get("h1d_ee_FG_unlike");
    h1d_drell_yan_temp->SetName("h1d_drell_yan_temp");
    h1d_drell_yan_temp->Scale(0.000246/10000.);
    h1d_drell_yan_temp->Scale(4.8572);

    TH1D* h1d_drell_yan = new TH1D("h1d_drell_yan","h1d_drell_yan", nbins_mass, mass_range);
    h1d_drell_yan->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_drell_yan_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_drell_yan_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_drell_yan_temp->GetBinContent(i)*h1d_drell_yan_temp->GetBinWidth(i);
            bin_error += pow(h1d_drell_yan_temp->GetBinError(i)*h1d_drell_yan_temp->GetBinWidth(i), 2);

        }

        h1d_drell_yan->SetBinContent(ibin+1, bin_content/h1d_drell_yan->GetBinWidth(ibin+1));
        h1d_drell_yan->SetBinError(ibin+1, sqrt(bin_error)/h1d_drell_yan->GetBinWidth(ibin+1));


    }


    h1d_drell_yan->SetLineColor(kBlue);
    h1d_drell_yan->SetLineStyle(3);
    h1d_drell_yan->Draw("same hist");


    TH1D* h1d_mlp_ccbar_temp = new TH1D("h1d_mlp_ccbar_temp", "DCA Distribution", 200, 1.50, 3.50);

    std::ifstream inFile_data("mlp_full_ccbar.txt");
    double x1, y1, err_y1;

    while(!inFile_data.eof()){

        inFile_data >> x1 >> y1 >> err_y1;

        int binX = h1d_mlp_ccbar_temp->FindBin(x1);
        h1d_mlp_ccbar_temp->SetBinContent(binX, y1);
        h1d_mlp_ccbar_temp->SetBinError(binX, err_y1);

    }

    TH1D* h1d_mlp_ccbar = new TH1D("h1d_mlp_ccbar","DCA Distribution", nbins_mass, mass_range);

    h1d_mlp_ccbar->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        double bin_center = (mass_range[ibin] + mass_range[ibin+1])/2.;

        h1d_mlp_ccbar->SetBinContent(ibin+1, h1d_mlp_ccbar_temp->GetBinContent(h1d_mlp_ccbar_temp->FindBin(bin_center)));
        h1d_mlp_ccbar->SetBinError(ibin+1, h1d_mlp_ccbar_temp->GetBinError(h1d_mlp_ccbar_temp->FindBin(bin_center)));

    }

    nbins = h1d_mlp_ccbar->GetNbinsX();
    counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = h1d_mlp_ccbar->GetBinCenter(ibin);
        y[counter] = h1d_mlp_ccbar->GetBinContent(ibin);
        x_lohi[counter] = h1d_mlp_ccbar->GetBinWidth(ibin)/2;

        stat_err[counter] = h1d_mlp_ccbar->GetBinError(ibin);

        sys_err[counter] = stat_err[counter];

        counter++;


    }

    TGraphErrors* band_ccbar_stat = new TGraphErrors(counter, x, y, 0, stat_err);
    band_ccbar_stat->SetMarkerStyle(20);
    band_ccbar_stat->SetMarkerSize(0.9);
    band_ccbar_stat->SetMarkerColor(kRed);
    band_ccbar_stat->SetLineColor(kRed);

    TGraphErrors* band_ccbar_sys = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    band_ccbar_sys->SetFillColorAlpha(kRed, 0.3);
    band_ccbar_sys->SetLineColor(kRed);

    band_ccbar_stat->Draw("e3same");
    band_ccbar_sys->Draw("e3same");


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

    TGraphErrors* band_bbbar_stat = (TGraphErrors*)fbbbar->Get("band");
    band_bbbar_stat->SetName("band_bbbar_stat");

    TFile* f5 = new TFile("normalizing_two_TGraphs_bbbar_sys.root","READ"); 
    TGraphErrors* band_bbbar_sys = (TGraphErrors*)f5->Get("band");
    band_bbbar_sys->SetName("band_bbbar_sys");

    band_bbbar_stat->SetFillColorAlpha(kBlue-4, 0.3);
    band_bbbar_stat->SetLineColor(kBlue-4);
    band_bbbar_stat->SetMarkerSize(0.0);
    band_bbbar_stat->SetLineWidth(0);

    band_bbbar_sys->SetFillColorAlpha(kBlue-4, 0.3);
    band_bbbar_sys->SetFillStyle(3022);
    band_bbbar_sys->SetMarkerSize(0.0);
    band_bbbar_sys->SetLineWidth(0);

    band_bbbar_stat->Draw("e3same");
    band_bbbar_sys->Draw("e3same");
    

    TH1D* h1d_mlp_dataBG_temp = new TH1D("h1d_mlp_dataBG_temp", "DCA Distribution", 200, 0.9, 4.50);

    std::ifstream inFile_data_BG("mlp_DCA_Fits_BG.txt");

    while(!inFile_data_BG.eof()){

        inFile_data_BG >> x1 >> y1 >> err_y1;

        int binX = h1d_mlp_dataBG_temp->FindBin(x1);
        h1d_mlp_dataBG_temp->SetBinContent(binX, y1);
        h1d_mlp_dataBG_temp->SetBinError(binX, err_y1);

    }

    TH1D* h1d_mlp_dataBG = new TH1D("h1d_mlp_dataBG","DCA Distribution", nbins_mass, mass_range);

    h1d_mlp_dataBG->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        double bin_center = (mass_range[ibin] + mass_range[ibin+1])/2.;

        h1d_mlp_dataBG->SetBinContent(ibin+1, h1d_mlp_dataBG_temp->GetBinContent(h1d_mlp_dataBG_temp->FindBin(bin_center)));
        h1d_mlp_dataBG->SetBinError(ibin+1, h1d_mlp_dataBG_temp->GetBinError(h1d_mlp_dataBG_temp->FindBin(bin_center)));

    }

    h1d_mlp_dataBG->Draw("same");

    TH1D* hist_total = (TH1D*)hist[0]->Clone();
    hist_total->SetName("hist_total");

    for(int ifile = 1; ifile < nfiles; ifile++){  hist_total->Add(hist[ifile]);}

    hist_total->Add(h1d_drell_yan);
    hist_total->Add(h1d_mlp_dataBG);
    hist_total->Add(h1d_mlp_ccbar);
    hist_total->Add(h1d_bbbar);

    nbins = hist_total->GetNbinsX();
    counter = 0;

    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = hist_total->GetBinCenter(ibin);
        y[counter] = hist_total->GetBinContent(ibin);
        x_lohi[counter] = hist_total->GetBinWidth(ibin)/2;

        stat_err[counter] = hist_total->GetBinError(ibin);

        double sys_jpsi = fcombined_jpsi->Eval(hist[0]->GetBinCenter(ibin))*hist[0]->GetBinContent(ibin);
        double sys_psip = fcombined_psip->Eval(hist[1]->GetBinCenter(ibin))*hist[1]->GetBinContent(ibin);
        double sys_bg = h1d_mlp_dataBG->GetBinContent(ibin)*0.05212;

        sys_err[counter] = sqrt( pow(sys_jpsi, 2) + pow(sys_psip, 2) + pow(sys_bg, 2) );

        counter++;


    }

    TGraphErrors* graph_stat_cocktail = new TGraphErrors(counter, x, y, 0, stat_err);
    graph_stat_cocktail->SetMarkerStyle(20);
    graph_stat_cocktail->SetMarkerSize(0.9);
    graph_stat_cocktail->SetMarkerColor(kRed);
    graph_stat_cocktail->SetLineColor(kRed);

    TGraphErrors* graph_sys_cocktail = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    graph_sys_cocktail->SetFillColorAlpha(kRed, 0.3);
    graph_sys_cocktail->SetLineColor(kRed);

    graph_stat_cocktail->Draw("psame");
    graph_sys_cocktail->Draw("e2same");

    TLegend* leg = new TLegend(0.2923706,0.6873418,0.7717984,0.8873418);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05164557);
    leg->SetNColumns(2);
    leg->AddEntry(mg, "FG_{e^{+}e^{-}}");
    leg->AddEntry(hist[0], "J/#psi #rightarrow e^{+}e^{-}");
    leg->AddEntry(hist[1], "#phi' #rightarrow e^{+}e^{-}");
    leg->AddEntry(band_ccbar_stat, "c#bar{c} #rightarrow e^{+}e^{-} (PYTHIA8) Stat.");
    leg->AddEntry(band_ccbar_sys, "c#bar{c} #rightarrow e^{+}e^{-} (PYTHIA8) Sys.");
    leg->AddEntry(band_bbbar_stat, "b#bar{b} #rightarrow e^{+}e^{-} (PYTHIA8) Stat.");
    leg->AddEntry(band_bbbar_sys, "b#bar{b} #rightarrow e^{+}e^{-} (PYTHIA8) Sys.");
    leg->AddEntry(h1d_drell_yan, "DY #rightarrow e^{+}e^{-} (PYTHIA8)");
    leg->AddEntry(graph_stat_cocktail, "Cocktail Sum");
    leg->Draw("same");

    nbins = hist_total->GetNbinsX();
    counter = 0;


    for(int ibin=1; ibin < nbins+1; ibin++){

        x[counter] = hist_total->GetBinCenter(ibin);
        y[counter] = h1d_data->GetBinContent(h1d_data->FindBin(x[counter]))/hist_total->GetBinContent(ibin);

        x_lohi[counter] = hist_total->GetBinWidth(ibin)/2;

        stat_err[counter] = y[counter]*sqrt( pow(h1d_data->GetBinError(h1d_data->FindBin(x[counter]))/h1d_data->GetBinContent(h1d_data->FindBin(x[counter])), 2) + pow(hist_total->GetBinError(ibin)/hist_total->GetBinContent(ibin), 2) );

        double sys_err_data = fcombined->Eval(h1d_data->GetBinCenter(h1d_data->FindBin(x[counter])))*h1d_data->GetBinContent(h1d_data->FindBin(x[counter]));

        double sys_jpsi = fcombined_jpsi->Eval(hist[0]->GetBinCenter(ibin))*hist[0]->GetBinContent(ibin);
        double sys_psip = fcombined_psip->Eval(hist[1]->GetBinCenter(ibin))*hist[1]->GetBinContent(ibin);
        double sys_bg = h1d_mlp_dataBG->GetBinContent(ibin)*0.05212;

        double sys_err_cocktail = sqrt( pow(sys_jpsi, 2) + pow(sys_psip, 2) + pow(sys_bg, 2) );
  
        sys_err[counter] = y[counter]*sqrt( pow(sys_err_data/h1d_data->GetBinContent(h1d_data->FindBin(x[counter])), 2) + pow(sys_err_cocktail/hist_total->GetBinContent(ibin), 2) );

        counter++;


    }
    

    TGraphErrors* graph_stat_ratio = new TGraphErrors(counter, x, y, 0, stat_err);
    graph_stat_ratio->SetMarkerStyle(20);
    graph_stat_ratio->SetMarkerSize(0.9);
    graph_stat_ratio->SetMarkerColor(kRed);
    graph_stat_ratio->SetLineColor(kRed);

    TGraphErrors* graph_sys_ratio = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    graph_sys_ratio->SetFillColorAlpha(kRed, 0.3);
    graph_sys_ratio->SetLineColor(kRed);

    pad2->cd();
    graph_stat_ratio->Draw("AP");
    graph_sys_ratio->Draw("e2same");

    
    TFile* fout = new TFile("dummy.root","RECREATE");
    fout->cd();
    h1d_data->Write();
    graph_stat->Write();
    graph_sys->Write();
    for(int ifile = 0; ifile < nfiles; ifile++){  hist[ifile]->Write(); }
    h1d_drell_yan->Write();
    band_ccbar_stat->Write();
    band_ccbar_sys->Write();
    band_bbbar_stat->Write();
    band_bbbar_sys->Write();
    graph_stat_cocktail->Write();
    graph_sys_cocktail->Write();
    mg->Write();
    c1->Write();
    fout->Close();

} 
