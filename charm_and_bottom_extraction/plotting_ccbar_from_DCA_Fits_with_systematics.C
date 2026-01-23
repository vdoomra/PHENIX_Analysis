
const int npts = 7;
double x_lohi[npts] = {-999};
double x[npts] = {-999};
double y[npts] = {-999};
double stat_err[npts] = {-999};
double sys_err[npts] = {-999};


void plotting_ccbar_from_DCA_Fits_with_systematics(){

    TFile* fdata_central = new TFile("plotting_background_bbbar_central.root","READ");
    TGraphErrors* graph_central = (TGraphErrors*)fdata_central->Get("g2");
    graph_central->SetName("graph_central");

    TF1* func = new TF1("func","pol1", 1.0, 4.5);
    //func->SetParameters(0.2657, -0.0786);  //ccbar systematic parameters
    func->SetParameters(0.4891, -0.1212); //bbbar systematic parameters

    TF1* fcorrection = new TF1("fcorrection","pol3", 1.0, 4.5);
    //fcorrection->SetParameters(0.0591339, -0.0391315, 0.0196653, -0.00210255); // charm correction
    fcorrection->SetParameters(0.0736545, -0.0338781, 0.0137945, -0.00120594);

    int nbins = graph_central->GetN();
    int counter = 0;

    for(int ibin=0; ibin < nbins; ibin++){

        double x_central, y_central, err_y_central;
        graph_central->GetPoint(counter, x_central, y_central);
        err_y_central = graph_central->GetErrorY(counter)/fcorrection->Eval(x_central);

        x[counter] = x_central;
        y[counter] = y_central/fcorrection->Eval(x_central);
        x_lohi[counter] = 0.05;

        stat_err[counter] = err_y_central;
        sys_err[counter] = func->Eval(x_central)*y[counter];

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
    mg->SetName("mg");
    
    mg->Add(graph_stat,"AP");
    mg->Add(graph_sys, "e2same");

    TCanvas* c1 = new TCanvas();
    mg->Draw("a");

    mg->GetXaxis()->SetTitle("m_{ee}[GeV/c^{2}]");
    mg->GetYaxis()->SetTitle("#frac{1}{N_{events}}#frac{#epsilon_{BBC}#sigma_{pp}}{#epsilon_{bias}}#frac{dN}{dm_{ee}} [GeV/c^{2}]^{-1} [In PHENIX Acceptance]");

    TLatex *   tex1 = new TLatex(1.6,7e-8,"Run 2015 p+p, #sqrt{s} = 200 GeV, |#eta| < 0.35");
    tex1->SetTextAlign(13);
    tex1->SetTextSize(0.052);
    tex1->SetLineWidth(1);
    tex1->SetTextFont(42);
    tex1->Draw();

    TLatex *   tex4 = new TLatex(2.6,3e-9, "p_{T}^{pair} < 5 GeV/c");
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

    TLegend* leg = new TLegend(0.2923706,0.6873418,0.7717984,0.8873418);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05164557);
    leg->SetNColumns(2);
    leg->AddEntry(mg, "b#bar{b} #rightarrow e^{+}e^{-}");
    leg->Draw("same");

    TFile* fout = new TFile("plotting_bbbar_from_DCA_Fits_with_systematics_corrected.root","RECREATE");
    fout->cd();
    graph_stat->Write();
    graph_sys->Write();
    mg->Write();
    fout->Close();



    
}