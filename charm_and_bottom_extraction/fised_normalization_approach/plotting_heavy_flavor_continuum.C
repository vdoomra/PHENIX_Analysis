
const int npts = 8;
double x_lohi[npts] = {-999};
double x[npts] = {-999};
double y[npts] = {-999};
double stat_err[npts] = {-999};
double sys_err[npts] = {-999};

void plotting_heavy_flavor_continuum(){

    TFile* f1 = new TFile("heavy_flavor_continuum_stat.root","READ");
    TGraphErrors* graph_stat = (TGraphErrors*)f1->Get("graph_stat");
    TGraphErrors* graph_sys = (TGraphErrors*)f1->Get("graph_sys");

    graph_stat->SetMarkerStyle(20);
    graph_stat->SetMarkerSize(0.9);
    graph_stat->SetMarkerColor(kBlack);
    graph_stat->SetLineColor(kBlack);

    graph_sys->SetFillColorAlpha(kBlack, 0.3);
    graph_sys->SetLineColor(kBlack);


    TMultiGraph* mg = new TMultiGraph();
    mg->SetName("mg");
    
    mg->Add(graph_stat,"AP");
    mg->Add(graph_sys, "e2same");

    TCanvas* c1 = new TCanvas();
    mg->Draw("a");

    mg->GetXaxis()->SetTitle("m_{ee}[GeV/c^{2}]");
    mg->GetYaxis()->SetTitle("#frac{1}{N_{events}}#frac{#epsilon_{BBC}#sigma_{pp}}{#epsilon_{bias}}#frac{dN}{dm_{ee}} [GeV/c^{2}]^{-1} [In PHENIX Acceptance]");

    TLatex *   tex1 = new TLatex(2.0,5e-7,"Run 2015 p+p, #sqrt{s} = 200 GeV, |#eta| < 0.35");
    tex1->SetTextAlign(13);
    tex1->SetTextSize(0.052);
    tex1->SetLineWidth(1);
    tex1->SetTextFont(42);
    tex1->Draw();

    TLatex *   tex4 = new TLatex(3.5,1e-8, "p_{T}^{pair} < 5 GeV/c");
    tex4->SetTextAlign(13);
    tex4->SetTextSize(0.052);
    tex4->SetLineWidth(1);
    tex4->SetTextFont(42);
    tex4->Draw();

    TLatex *   tex5 = new TLatex(3.4,5e-8, "p_{T}^{trig, e} > 1.5 GeV/c");
    tex5->SetTextAlign(13);
    tex5->SetTextSize(0.052);
    tex5->SetLineWidth(1);
    tex5->SetTextFont(42);
    tex5->Draw();

    TGraphErrors* band_stat = (TGraphErrors*)f1->Get("band");
    band_stat->SetName("band_stat");

    TFile* f2 = new TFile("heavy_flavor_continuum_sys.root","READ"); 
    TGraphErrors* band_sys = (TGraphErrors*)f2->Get("band");
    band_sys->SetName("band_sys");

    band_stat->SetFillColorAlpha(kBlack, 0.3);
    band_stat->SetMarkerSize(0);
    band_stat->SetLineWidth(0);

    band_sys->SetFillColorAlpha(kBlack, 0.3);
    band_sys->SetFillStyle(3022);
    band_stat->SetMarkerSize(0);
    band_sys->SetLineWidth(0);

    band_stat->Draw("e3same");
    band_sys->Draw("e3same");

    TH1D* hist_ccbar_dataBinning_central = (TH1D*)f1->Get("hist_ccbar_dataBinning_central");
    hist_ccbar_dataBinning_central->SetLineColor(kGreen-4);
    hist_ccbar_dataBinning_central->SetLineStyle(2);
    hist_ccbar_dataBinning_central->Draw("same hist");

    TH1D* hist_bbbar_dataBinning_central = (TH1D*)f1->Get("hist_bbbar_dataBinning_central");
    hist_bbbar_dataBinning_central->SetLineColor(kGreen-4);
    hist_bbbar_dataBinning_central->SetLineStyle(2);
    hist_bbbar_dataBinning_central->Draw("same hist");

    TH1D* hist_drell_yan_dataBinning_central = (TH1D*)f1->Get("hist_drell_yan_dataBinning_central");
    hist_drell_yan_dataBinning_central->SetLineColor(kGreen-4);
    hist_drell_yan_dataBinning_central->SetLineStyle(2);
    hist_drell_yan_dataBinning_central->Draw("same hist");



    TLegend* leg = new TLegend(0.2923706,0.6873418,0.7717984,0.8873418);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05164557);
    leg->SetNColumns(2);
    leg->AddEntry(graph_stat, "Data");
    leg->AddEntry(band_stat, "c#bar{c} + b#bar{b} + DY #rightarrow e^{+}e^{-} (PYTHIA8) Stat.");
    leg->AddEntry(band_sys, "c#bar{c} + b#bar{b} + DY #rightarrow e^{+}e^{-} (PYTHIA8) Sys.");
    leg->AddEntry(hist_ccbar_dataBinning_central, "c#bar{c} #rightarrow e^{+}e^{-} (PYTHIA8)");
    leg->AddEntry(hist_bbbar_dataBinning_central, " b#bar{b} #rightarrow e^{+}e^{-} (PYTHIA8)");
    leg->AddEntry(hist_drell_yan_dataBinning_central, "DY #rightarrow e^{+}e^{-} (PYTHIA8)");
    leg->Draw("same");

    
}