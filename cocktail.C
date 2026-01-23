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

const int nfiles = 10;
const char* fileNames[nfiles] = { "pi0_gee_2k", "eta_gee_1999" ,"etap_gee_2k", "omega_ee_1998", "omega_pi0ee_2k", "phi_ee_2k" ,"phi_etaee_1992", "jpsi_2k", "psip_2k", "rho_1972"};
const double num[nfiles] = {2000., 1999., 2000., 1998., 2000., 2000., 1992., 2000., 2000., 1972.};
const double BR[nfiles] = {1-0.98823, 0.0069, 0.000491, 0.0000739, 0.00077, 0.000297, 0.000108, 0.05971, 7.93e-3, 4.72e-5};
const int color[nfiles] = {kRed, kBlue, kGreen-7, kViolet, kBlack, kCyan-4, kPink, kAzure, kRed, kBlue};
const int line_style[nfiles] = {1,1,1,1,1,1,1,1,2,2};

TH1D* hist[nfiles] = {NULL};
TFile* fin = NULL;

const double nevents = 4.216682807*pow(10,11)*(0.79);

double combinedFunction(double *x, double *par) {
    double fERT_val = par[0] + par[1] * x[0] + par[2]*x[0]*x[0];
    if(x[0] > 3.50) fERT_val = 0.0;
    double fn0_val = par[3];
    double feoverp_val = par[4] + par[5]*x[0] + par[6]*x[0]*x[0];
    if(x[0] > 3.20) feoverp_val = 0.03262;

    // 2% uncertainty on the number of events

    return TMath::Sqrt(fERT_val * fERT_val + fn0_val * fn0_val + feoverp_val * feoverp_val + 0.02*0.02);
}

void cocktail(){


    TFile* f1 = new TFile("/Users/vassu/Desktop/new_files/data_files/h1d_ee_FG_unlike.root","READ");
    TH1D* h1d_data = (TH1D*)f1->Get("h1d_ee_FG_cent0");
    h1d_data->SetName("h1d_data");

    TF1* fcombined = new TF1("fcombined", combinedFunction, 0.0, 5.5, 7);
    fcombined->SetParameters(0.102, -0.0022, -0.00492, 0.06424, 0.06622, 0.08956, -0.03127);

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

    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
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
        hist[ifile] = (TH1D*)fin->Get("h1d_ee_FG_unlike");
        hist[ifile]->SetName(Form("h1d_ee_%s",fileNames[ifile]));

        hist[ifile]->Scale((42.*BR[ifile])/(25000.*num[ifile]));
        hist[ifile]->SetLineColor(color[ifile]);
        hist[ifile]->SetLineStyle(line_style[ifile]);

        hist[ifile]->Draw("same hist");

    }

    TFile* fdrellyan = new TFile("/Users/vassu/Desktop/new_files/data_files/dca_output_drell_yan_10k_gencut.root","READ");
    TH1D* h1d_drell_yan = (TH1D*)fdrellyan->Get("h1d_ee_FG_unlike");
    h1d_drell_yan->SetName("h1d_drell_yan");
    h1d_drell_yan->Scale(0.000246/10000.); // Multiplied by Cross section

    h1d_drell_yan->SetLineColor(kBlue);
    h1d_drell_yan->SetLineStyle(3);
    h1d_drell_yan->Draw("same hist");

    TH1D* hist_total = (TH1D*)hist[0]->Clone();
    hist_total->SetName("hist_total");

    for(int ifile = 1; ifile < nfiles; ifile++){  hist_total->Add(hist[ifile]);}

    hist_total->Add(h1d_drell_yan);

    hist_total->SetLineStyle(1);
    hist_total->SetLineWidth(3);
    hist_total->SetLineColor(kViolet);
    hist_total->Draw("same hist");

    TLegend* leg = new TLegend(0.2923706,0.6873418,0.7717984,0.8873418);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05164557);
    leg->SetNColumns(2);
    leg->AddEntry(mg, "FG_{e^{+}e^{-}}");
    leg->AddEntry(hist[0], "#pi^{0} #rightarrow #gamma e^{+}e^{-}");
    leg->AddEntry(hist[1], "#eta #rightarrow #gamma e^{+}e^{-}");
    leg->AddEntry(hist[2], "#eta' #rightarrow #gamma e^{+}e^{-}");
    leg->AddEntry(hist[3], "#omega #rightarrow e^{+}e^{-}");
    leg->AddEntry(hist[4], "#omega #rightarrow #pi^{0} e^{+}e^{-}");
    leg->AddEntry(hist[5], "#phi #rightarrow e^{+}e^{-}");
    leg->AddEntry(hist[6], "#phi #rightarrow #eta e^{+}e^{-}");
    leg->AddEntry(hist[7], "J/#psi #rightarrow e^{+}e^{-}");
    leg->AddEntry(hist[8], "#psi' #rightarrow e^{+}e^{-}");
    leg->AddEntry(hist[9], "#rho #rightarrow e^{+}e^{-}");
    leg->AddEntry(h1d_drell_yan, "DY #rightarrow e^{+}e^{-} (PYTHIA8)");
    leg->AddEntry(hist_total, "Cocktail Sum");
    leg->Draw("same");
    
    TFile* fout = new TFile("cocktail_without_heavy_flavor.root","RECREATE");
    fout->cd();
    h1d_data->Write();
    for(int ifile = 0; ifile < nfiles; ifile++){  hist[ifile]->Write(); }
    h1d_drell_yan->Write();
    hist_total->Write();
    graph_stat->Write();
    graph_sys->Write();
    mg->Write();
    c1->Write();
    fout->Close();
    

    
}