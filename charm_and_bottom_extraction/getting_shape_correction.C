
//These regions are used for ccbar extraction using SoftQCD bbbar
//const int mass_regions = 6;
//double mass_lo[mass_regions] = {  1.30, 1.75, 2.20, 2.40, 1.60, 2.40};
//double mass_hi[mass_regions] = {  2.05, 2.40, 2.80, 2.95, 2.80, 3.45};

//These regions are used for bbbar extraction using SoftQCD ccbar as the intial guess
const int mass_regions = 7;
double mass_lo[mass_regions] = {1.75, 2.40, 1.60, 2.40, 2.60, 3.10, 2.80};
double mass_hi[mass_regions] = {2.40, 2.95, 2.80, 3.45, 4.00, 4.25, 5.00};

//These bins are used for BG extraction
//const int mass_regions = 13;
//double mass_lo[mass_regions] = { 0.70, 0.87, 1.30, 1.75, 2.20, 2.40, 1.60, 2.40, 3.30, 2.60, 3.10, 3.30, 1.15};
//double mass_hi[mass_regions] = { 1.15, 1.30, 2.05, 2.40, 2.80, 2.95, 2.80, 3.45, 4.50, 4.00, 4.25, 5.50, 1.75};

const int nbins_mass = 40;
const double mass_range[nbins_mass+1] = { 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75, 0.80, 0.87, 0.95, 1.05, 1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.05, 3.10, 3.15, 3.20, 3.30, 3.45, 3.60, 3.70, 3.80, 4.00, 4.25, 4.5, 5.0, 5.5};

double x[mass_regions] = {-999};
double y[mass_regions] = {-999};
double err_x[mass_regions] = {-999};
double err_y[mass_regions] = {-999};

double x_shape_correction[mass_regions] = {-999};
double y_shape_correction[mass_regions] = {-999};
double err_x_shape_correction[mass_regions] = {0};
double err_y_shape_correction[mass_regions] = {-999};

void getting_shape_correction(){

    TH1D* h1d_unlike_BG_mlp = new TH1D("h1d_unlike_BG_mlp", "h1d_unlike_BG_mlp", 200, 0.0, 5.5);

    std::ifstream inFile_BG("/Users/vassu/Desktop/new_files/data_files/MLP_files/mlp_comb_BG.txt");
    double x1, y1, err_y1;

    while(!inFile_BG.eof()){

        inFile_BG >> x1 >> y1 >> err_y1;

        int binX = h1d_unlike_BG_mlp->FindBin(x1);
        h1d_unlike_BG_mlp->SetBinContent(binX, y1);
        h1d_unlike_BG_mlp->SetBinError(binX, 0.0);

    }

    TH1D* h1d_unlike_BG = new TH1D("h1d_unlike_BG", "Unlike BG Distribution", nbins_mass, mass_range);
    h1d_unlike_BG->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        double bin_center = (mass_range[ibin] + mass_range[ibin+1])/2.;

        h1d_unlike_BG->SetBinContent(ibin+1, h1d_unlike_BG_mlp->GetBinContent(h1d_unlike_BG_mlp->FindBin(bin_center)));
        h1d_unlike_BG->SetBinError(ibin+1, 0.0);


    }


    for(int ibin = 0; ibin < mass_regions; ibin++){

        int bin_lo = h1d_unlike_BG->FindBin(mass_lo[ibin] + 0.01);
        int bin_hi = h1d_unlike_BG->FindBin(mass_hi[ibin] - 0.01);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_unlike_BG->GetBinContent(i)*h1d_unlike_BG->GetBinWidth(i);
            bin_error += pow(h1d_unlike_BG->GetBinError(i)*h1d_unlike_BG->GetBinWidth(i), 2);

        }

        x[ibin] = (mass_lo[ibin] + mass_hi[ibin])/2;
        y[ibin] = bin_content/(mass_hi[ibin] - mass_lo[ibin]);
        err_x[ibin] = (mass_hi[ibin] - mass_lo[ibin])/2;
        err_y[ibin] = sqrt(bin_error)/(mass_hi[ibin] - mass_lo[ibin]);


    }



    for(int ibin = 0; ibin < mass_regions; ibin++){


        x_shape_correction[ibin] = x[ibin];
        y_shape_correction[ibin] = (h1d_unlike_BG_mlp->GetBinContent(h1d_unlike_BG_mlp->FindBin(x[ibin]))/y[ibin]);
        err_y_shape_correction[ibin] = y_shape_correction[ibin]*sqrt( pow(err_y[ibin]/y[ibin], 2) + pow(h1d_unlike_BG_mlp->GetBinError(h1d_unlike_BG_mlp->FindBin(x[ibin]))/h1d_unlike_BG_mlp->GetBinContent(h1d_unlike_BG_mlp->FindBin(x[ibin])), 2) );

    }


    TGraphErrors* graph = new TGraphErrors(mass_regions, x, y, err_x, err_y);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kRed);
    graph->SetName("graph");

    TGraphErrors* graph_shape_correction = new TGraphErrors(mass_regions, x_shape_correction, y_shape_correction, 0, err_y_shape_correction);
    graph_shape_correction->SetMarkerStyle(20);
    graph_shape_correction->SetMarkerColor(kRed);
    graph_shape_correction->SetName("graph_shape_correction");

    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
    h1d_unlike_BG->Draw();
    graph->Draw("same");

    TCanvas* c2 = new TCanvas();
    graph_shape_correction->Draw("AP");


    TFile* fout = new TFile("getting_shape_correction_bbbar.root","RECREATE");
    fout->cd();
    graph->Write();
    h1d_unlike_BG->Write();
    h1d_unlike_BG_mlp->Write();
    graph_shape_correction->Write();
    fout->Close();


}