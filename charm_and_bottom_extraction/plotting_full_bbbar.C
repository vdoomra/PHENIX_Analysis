
const int npts = 3;
double x_lohi[npts] = {-999};
double x[npts] = {-999};
double y[npts] = {-999};
double stat_err[npts] = {-999};
double sys_err[npts] = {-999};

void plotting_full_bbbar(){

    std::ofstream datafile;
    datafile.open("plotting_full_bbbar_corrected.txt",std::ofstream::app);

    TFile* f1 = new TFile("plotting_bbbar_from_DCA_Fits_with_systematics_corrected.root","READ");
    TGraphErrors* graph_stat = (TGraphErrors*)f1->Get("graph_stat");
    TGraphErrors* graph_sys = (TGraphErrors*)f1->Get("graph_sys");

    graph_stat->SetMarkerStyle(20);
    graph_stat->SetMarkerSize(0.9);
    graph_stat->SetMarkerColor(kBlack);
    graph_stat->SetLineColor(kBlack);

    graph_sys->SetFillColorAlpha(kBlack, 0.3);
    graph_sys->SetLineColor(kBlack);

    for(int j = 0; j < (int)graph_stat->GetN(); j++){

        double x_graph, y_graph;
        graph_stat->GetPoint(j, x_graph, y_graph);

        datafile << x_graph << "\t" << y_graph << "\t" << graph_stat->GetErrorY(j) << "\t" << graph_sys->GetErrorY(j) << endl;

    }


    TF1* func = new TF1("func","pol2", 3.0, 4.5);
    func->SetParameters(6.44, -2.36, 0.264);

    TF1* fcorrection = new TF1("fcorrection","pol3", 1.0, 4.5);
    fcorrection->SetParameters(0.0736545, -0.0338781, 0.0137945, -0.00120594);

    TFile* f2 = new TFile("likesign_analysis_alternate_approach_stat.root","READ");
    TH1D* h1d_signal_stat = (TH1D*)f2->Get("h1d_signal_stat");
    TH1D* h1d_signal_sys = (TH1D*)f2->Get("h1d_signal_sys");


    int counter = 0;
    double x_graph, y_graph;

    for(int i = 1; i < (int)h1d_signal_stat->GetNbinsX()+1; i++){

        x[counter] = h1d_signal_stat->GetBinCenter(i);
        y[counter] = h1d_signal_stat->GetBinContent(i)*func->Eval(x[counter])/fcorrection->Eval(x[counter]);
        x_lohi[counter] = 0.05;
        stat_err[counter] = h1d_signal_stat->GetBinError(i)*func->Eval(x[counter])/fcorrection->Eval(x[counter]);
        sys_err[counter] = h1d_signal_sys->GetBinError(i)*func->Eval(x[counter])/fcorrection->Eval(x[counter]);

        datafile << x[counter] << "\t" << y[counter] << "\t" << stat_err[counter] << "\t" << sys_err[counter] << endl;

        counter++;

    }

    TGraphErrors* graph_stat_likesign = new TGraphErrors(counter, x, y, 0, stat_err);
    graph_stat_likesign->SetMarkerStyle(20);
    graph_stat_likesign->SetMarkerSize(0.9);
    graph_stat_likesign->SetMarkerColor(kBlue-4);
    graph_stat_likesign->SetLineColor(kBlue-4);
    graph_stat_likesign->SetName("graph_stat_likesign");

    TGraphErrors* graph_sys_likesign = new TGraphErrors(counter, x, y, x_lohi, sys_err);
    graph_sys_likesign->SetFillColorAlpha(kBlue-4, 0.3);
    graph_sys_likesign->SetLineColor(kBlue-4);
    graph_sys_likesign->SetName("graph_sys_likesign");


    TMultiGraph* mg = new TMultiGraph();
    mg->SetName("mg");
    
    mg->Add(graph_stat,"AP");
    mg->Add(graph_sys, "e2same");
    mg->Add(graph_stat_likesign,"psame");
    mg->Add(graph_sys_likesign, "e2same");

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

    TLatex *   tex5 = new TLatex(2.6,7e-9, "p_{T}^{trig, e} > 1.56 GeV/c (PbSc)");
    tex5->SetTextAlign(13);
    tex5->SetTextSize(0.052);
    tex5->SetLineWidth(1);
    tex5->SetTextFont(42);
    tex5->Draw();

    TLatex *   tex6 = new TLatex(2.6,2e-9, "p_{T}^{trig, e} > 2.40 GeV/c (PbGl)");
    tex6->SetTextAlign(13);
    tex6->SetTextSize(0.052);
    tex6->SetLineWidth(1);
    tex6->SetTextFont(42);
    tex6->Draw();

    TFile* f3 = new TFile("normalizing_two_TGraphs_bbbar_stat.root","READ");
    TGraphErrors* band_stat = (TGraphErrors*)f3->Get("band");
    band_stat->SetName("band_stat");

    TFile* f4 = new TFile("normalizing_two_TGraphs_bbbar_sys.root","READ"); 
    TGraphErrors* band_sys = (TGraphErrors*)f4->Get("band");
    band_sys->SetName("band_sys");

    band_stat->SetMarkerStyle(20);
    band_stat->SetMarkerSize(0.0);
    band_stat->SetLineWidth(0);

    band_sys->SetFillColorAlpha(kBlack, 0.3);
    band_sys->SetFillStyle(3022);
    band_stat->SetMarkerSize(0.0);
    band_sys->SetLineWidth(0);

    band_stat->Draw("e3same");
    band_sys->Draw("e3same");

    TLegend* leg = new TLegend(0.2923706,0.6873418,0.7717984,0.8873418);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05164557);
    leg->SetNColumns(2);
    leg->AddEntry(graph_stat, "b#bar{b} #rightarrow e^{+}e^{-} (From DCA Fits)");
    leg->AddEntry(graph_stat_likesign, "b#bar{b} #rightarrow e^{+}e^{-} (From LikeSign Mass Spectra)");
    leg->AddEntry(band_stat, "b#bar{b} #rightarrow e^{+}e^{-} (PYTHIA8) Stat.");
    leg->AddEntry(band_sys, "b#bar{b} #rightarrow e^{+}e^{-} (PYTHIA8) Sys.");
    leg->Draw("same");

    
}