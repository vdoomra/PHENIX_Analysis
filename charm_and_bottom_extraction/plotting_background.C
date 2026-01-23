
//These regions are used for ccbar extraction when I used SoftQCD for bbbar

//const int mass_regions = 6;
//double mass_lo[mass_regions] = {  1.30, 1.75, 2.20, 2.40, 1.60, 2.40};
//double mass_hi[mass_regions] = {  2.05, 2.40, 2.80, 2.95, 2.80, 3.45};

//const char* mass_lo_name[mass_regions] = { "130", "175", "220", "240", "160", "240"}; 
//const char* mass_hi_name[mass_regions] = { "205", "240", "280", "295", "280", "345"};
//const char* region_name[mass_regions] = {"Region4", "Region5", "Region6", "Region7", "Region8", "Region12"};

//These regions are used for bbbar extraction when I use ccbar from SoftQCD

const int mass_regions = 7;
double mass_lo[mass_regions] = { 1.75, 2.40, 1.60, 2.40, 2.60, 3.10, 2.80};
double mass_hi[mass_regions] = { 2.40, 2.95, 2.80, 3.45, 4.00, 4.25, 5.00};

const char* mass_lo_name[mass_regions] = { "175", "240", "160", "240", "260", "310", "280"}; 
const char* mass_hi_name[mass_regions] = { "240", "295", "280", "345", "400", "425", "500"};
const char* region_name[mass_regions] = { "Region5", "Region7", "Region8", "Region12", "Region13", "Region14", "Region15"};

//The following mass regions are used for BG extraction

//const int mass_regions = 13;
//double mass_lo[mass_regions] = { 0.70, 0.87, 1.30, 1.75, 2.20, 2.40, 1.60, 2.40, 3.30, 2.60, 3.10, 3.30, 1.15};
//double mass_hi[mass_regions] = { 1.15, 1.30, 2.05, 2.40, 2.80, 2.95, 2.80, 3.45, 4.50, 4.00, 4.25, 5.50, 1.75};

//const char* mass_lo_name[mass_regions] = { "070", "087", "130", "175", "220", "240", "160", "240", "330", "260", "310", "330", "115"}; 
//const char* mass_hi_name[mass_regions] = { "115", "130", "205", "240", "280", "295", "280", "345", "450", "400", "425", "550", "175"};
//const char* region_name[mass_regions] = {"Region1", "Region2", "Region3", "Region4", "Region5", "Region6", "Region7", "Region8", "Region9", "Region10", "Region11", "Region12", "Region13"};

double x[mass_regions] = {-999};
double y[mass_regions] = {-999};
double y_shape_correction[mass_regions] = {-999};
double err_x[mass_regions] = {-999};
double err_y[mass_regions] = {-999};
double err_y_shape_correction[mass_regions] = {-999};

TF1* func = NULL;

void plotting_background(){


    TFile* fdata = new TFile("/Users/vassu/Desktop/new_files/data_files/producing_templates.root","READ");
    TH1D* hdata = (TH1D*)fdata->Get("h1d_ee_FG_unlike_original");

    TFile* f1 = new TFile("template_fits_prompt_ccbar_fixed_upper.root","READ");
    TCanvas* c1 = new TCanvas();
    c1->Divide(4,2);

    TFile* f2 = new TFile("getting_shape_correction_bbbar.root","READ");
    TGraphErrors* graph_shape_correction = (TGraphErrors*)f2->Get("graph_shape_correction");

    int counter = 0;
    double x_graph, y_graph;

    for(int iregion = 0; iregion < mass_regions; iregion++){

        c1->cd(counter+1);

        TH1D* hist = (TH1D*)f1->Get(Form("h1d_para_bbbar_%s", region_name[iregion]));
        hist->SetName(Form("h1d_para_bbbar_%s", region_name[iregion]));

        func = new TF1("func","gaus", hist->GetMean()-1.6*hist->GetStdDev(), hist->GetMean()+1.6*hist->GetStdDev());
        //if(iregion == 0) func = new TF1("func","gaus", 0.10, 0.35);
        //if(iregion == 1) func = new TF1("func","gaus", 0.25, 0.50);
        //if(iregion == 5) func = new TF1("func","gaus", 0.12, 0.35);
        
        hist->Fit(func,"RQ+");

        hist->Draw();
        hist->GetXaxis()->SetTitle("Fraction");
        hist->GetYaxis()->SetTitle("Counts [Arbitrary Units]");

        TLatex *   tex1 = new TLatex(0.5,2000, Form("%.02f < m_{ee} < %.02f GeV/c^{2}", mass_lo[iregion], mass_hi[iregion]));
        tex1->SetTextAlign(13);
        tex1->SetTextSize(0.052);
        tex1->SetLineWidth(1);
        tex1->SetTextFont(42);
        tex1->Draw();

        int bin_lo = hdata->FindBin(mass_lo[iregion] + 0.01);
        int bin_hi = hdata->FindBin(mass_hi[iregion] - 0.01);

        double bin_content = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += hdata->GetBinContent(i)*hdata->GetBinWidth(i);

        }

        graph_shape_correction->GetPoint(counter, x_graph, y_graph);

        x[iregion] = (mass_lo[iregion] + mass_hi[iregion])/2;
        y[iregion] = (func->GetParameter(1)*bin_content)/(mass_hi[iregion] - mass_lo[iregion]);
        y_shape_correction[iregion] = y_graph*(func->GetParameter(1)*bin_content)/(mass_hi[iregion] - mass_lo[iregion]);
        err_x[iregion] = 0.0;
        err_y[iregion] = bin_content*func->GetParameter(2)/(mass_hi[iregion] - mass_lo[iregion]);
        err_y_shape_correction[iregion] = y_graph*bin_content*func->GetParameter(2)/(mass_hi[iregion] - mass_lo[iregion]);

        counter++;

    }


    std::ofstream datafile;
    datafile.open("dummy.txt",std::ofstream::app);

    for(int j = 0; j < mass_regions; j++){

        datafile << x[j] << "\t" << y_shape_correction[j] << "\t" << err_y_shape_correction[j] << endl;

    }


    TGraphErrors* g1 = new TGraphErrors(mass_regions, x, y, err_x, err_y);
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kRed);
    g1->SetName("g1");

    TGraphErrors* g2 = new TGraphErrors(mass_regions, x, y_shape_correction, err_x, err_y_shape_correction);
    g2->SetMarkerStyle(20);
    g2->SetMarkerColor(kRed);
    g2->SetName("g2");

    TCanvas* c2 = new TCanvas();
    c2->SetLogy();
    g1->Draw("AP");
    g2->Draw("same p");

    TCanvas* c3 = new TCanvas();
    c3->SetLogy();
    hdata->Draw();
    g2->Draw("samep");

    TFile* fout = new TFile("plotting_background_bbbar_upper.root","RECREATE");
    fout->cd();
    c1->Write();
    g1->Write();
    g2->Write();
    fout->Close();
    
}