
const int npts = 4;
double x_lohi[npts] = {-999};
double x[npts] = {-999};
double y[npts] = {-999};
double stat_err[npts] = {-999};
double sys_err[npts] = {-999};

void plotting_full_ccbar(){

    std::ofstream datafile;
    datafile.open("plotting_full_ccbar.txt",std::ofstream::app);


    TFile* f1 = new TFile("obtaining_ccbar_from_mass_space.root","READ");
    TGraphErrors* graph_stat = (TGraphErrors*)f1->Get("graph_stat");
    TGraphErrors* graph_sys = (TGraphErrors*)f1->Get("graph_sys");

    graph_stat->SetMarkerStyle(20);
    graph_stat->SetMarkerSize(0.9);
    graph_stat->SetMarkerColor(kBlue-4);
    graph_stat->SetLineColor(kBlue-4);

    graph_sys->SetFillColorAlpha(kBlue-4, 0.3);
    graph_sys->SetLineColor(kBlue-4);

    for(int j = 0; j < (int)graph_stat->GetN(); j++){

        double x_graph, y_graph;
        graph_stat->GetPoint(j, x_graph, y_graph);

    }


    TFile* f2 = new TFile("plotting_ccbar_from_DCA_Fits_with_systematics.root","READ");
    TGraphErrors* graph_stat_dca = (TGraphErrors*)f2->Get("graph_stat");
    graph_stat_dca->SetName("graph_stat_dca");
    TGraphErrors* graph_sys_dca = (TGraphErrors*)f2->Get("graph_sys");
    graph_sys_dca->SetName("graph_sys_dca");

    
    graph_stat_dca->SetMarkerStyle(20);
    graph_stat_dca->SetMarkerSize(0.9);
    graph_stat_dca->SetMarkerColor(kBlack);
    graph_stat_dca->SetLineColor(kBlack);

    graph_sys_dca->SetFillColorAlpha(kBlack, 0.3);
    graph_sys_dca->SetLineColor(kBlack);

    for(int j = 0; j < (int)graph_stat_dca->GetN(); j++){

        double x_graph, y_graph;
        graph_stat_dca->GetPoint(j, x_graph, y_graph);

        datafile << x_graph << "\t" << y_graph << "\t" << graph_stat_dca->GetErrorY(j) << "\t" << graph_sys_dca->GetErrorY(j) << endl;

    }

    TMultiGraph* mg = new TMultiGraph();
    mg->SetName("mg");
    
    mg->Add(graph_stat,"AP");
    mg->Add(graph_sys, "e2same");
    mg->Add(graph_stat_dca,"psame");
    mg->Add(graph_sys_dca, "e2same");

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

    TFile* f3 = new TFile("normalizing_two_TGraphs_ccbar_stat.root","READ");
    TGraphErrors* band_stat = (TGraphErrors*)f3->Get("band");
    band_stat->SetName("band_stat");

    TFile* f4 = new TFile("normalizing_two_TGraphs_ccbar_sys.root","READ"); 
    TGraphErrors* band_sys = (TGraphErrors*)f4->Get("band");
    band_sys->SetName("band_sys");

    band_stat->SetFillColorAlpha(kBlack, 0.3);
    band_stat->SetLineColor(kBlack);
    band_stat->SetMarkerSize(0.0);
    band_stat->SetLineWidth(0);

    band_sys->SetFillColorAlpha(kBlack, 0.3);
    band_sys->SetFillStyle(3022);
    band_sys->SetMarkerSize(0.0);
    band_sys->SetLineWidth(0);

    band_stat->Draw("e3same");
    band_sys->Draw("e3same");

    TLegend* leg = new TLegend(0.2923706,0.6873418,0.7717984,0.8873418);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05164557);
    leg->SetNColumns(1);
    leg->AddEntry(graph_stat_dca, "c#bar{c} #rightarrow e^{+}e^{-} (From DCA Fits)");
    leg->AddEntry(graph_stat, "c#bar{c} #rightarrow e^{+}e^{-} (From Unlike-Sign Mass Spectra)");
    leg->AddEntry(band_stat, "c#bar{c} #rightarrow e^{+}e^{-} (PYTHIA8) Stat.");
    leg->AddEntry(band_sys, "c#bar{c} #rightarrow e^{+}e^{-} (PYTHIA8) Sys.");
    leg->Draw("same");
    
}