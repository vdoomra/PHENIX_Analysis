
const int nbins_mass = 26;
const double mass_range[nbins_mass+1] = {0.87, 0.95, 1.05, 1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.05, 3.10, 3.15, 3.20, 3.30, 3.45, 3.60, 3.70, 3.80, 4.00, 4.25, 4.5};

const double nevents = 4.216682807*pow(10,11)*(0.79);

void getting_final_background(){

    TFile* f1 = new TFile("/Users/vassu/Desktop/new_files/data_files/h1d_ee_FG_unlike.root","READ");
    TH1D* h1d_data_temp = (TH1D*)f1->Get("h1d_ee_FG_cent0");

    TH1D* h1d_data = new TH1D("h1d_data","h1d_data", nbins_mass, mass_range);
    h1d_data->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        int bin_lo = h1d_data_temp->FindBin(mass_range[ibin] + 0.001);
        int bin_hi = h1d_data_temp->FindBin(mass_range[ibin+1] - 0.001);

        double bin_content = 0; double bin_error = 0;

        for(int i = bin_lo; i <= bin_hi; i++){

            bin_content += h1d_data_temp->GetBinContent(i)*h1d_data_temp->GetBinWidth(i);
            bin_error += pow(h1d_data_temp->GetBinError(i)*h1d_data_temp->GetBinWidth(i), 2);

        }

        h1d_data->SetBinContent(ibin+1, bin_content/h1d_data->GetBinWidth(ibin+1));
        h1d_data->SetBinError(ibin+1, sqrt(bin_error)/h1d_data->GetBinWidth(ibin+1));


    }

    TH1D* h1d_mlp_dataBG_temp = new TH1D("h1d_mlp_dataBG_temp", "DCA Distribution", 200, 0.9, 4.50);

    std::ifstream inFile_data("mlp_DCA_Fits_BG.txt");
    double x1, y1, err_y1;

    while(!inFile_data.eof()){

        inFile_data >> x1 >> y1 >> err_y1;

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

    TH1D* h1d_mlp_like = new TH1D("h1d_mlp_like","Like Sign Distribution", 200, 0.5, 5.5);
    h1d_mlp_like->Sumw2();

    std::ifstream inFile("/Users/vassu/Desktop/new_files/data_files/mlp_likesign_mass.txt");
    double x_mlp, y_mlp, err_y_mlp;

    while(!inFile.eof()){

        inFile >> x_mlp >> y_mlp >> err_y_mlp;

        int binX = h1d_mlp_like->FindBin(x_mlp);
        h1d_mlp_like->SetBinContent(binX, y_mlp);
        h1d_mlp_like->SetBinError(binX, err_y_mlp);

    }

    TH1D* h1d_like_combined = new TH1D("h1d_like_combined","Like Sign Distribution", nbins_mass, mass_range);
    h1d_like_combined->Sumw2();

    for(int ibin=0; ibin < nbins_mass; ibin++){

        double bin_center = (mass_range[ibin] + mass_range[ibin+1])/2.;

        double bin_content = h1d_mlp_like->GetBinContent(h1d_mlp_like->FindBin(bin_center));
        double bin_error = h1d_mlp_like->GetBinError(h1d_mlp_like->FindBin(bin_center));

        h1d_like_combined->SetBinContent(ibin+1, bin_content);
        h1d_like_combined->SetBinError(ibin+1, bin_error);


    }

    h1d_like_combined->Scale(1./nevents);
    h1d_like_combined->Scale(1.0726*0.545*42.);

    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
    h1d_data->Draw();
    h1d_mlp_dataBG->Draw("same");
    h1d_like_combined->Draw("same");

    const int nBins = h1d_mlp_dataBG_temp->GetNbinsX();
    std::vector<double> x(nBins), y(nBins), ex(nBins), ey(nBins);

    int counter = 0;

    for (int i = 0; i < nBins; ++i) {
        int bin = i + 1;
        x[i] = h1d_mlp_dataBG_temp->GetBinCenter(bin);
        y[i] = h1d_mlp_dataBG_temp->GetBinContent(bin);
        if(y[i] == 0) continue;
        ex[i] = h1d_mlp_dataBG_temp->GetBinWidth(bin) / 2.0;
        ey[i] = h1d_mlp_dataBG_temp->GetBinError(bin);

        counter++;
   
    }

    TGraphErrors* g_band = new TGraphErrors(counter, &x[0], &y[0], &ex[0], &ey[0]);
    g_band->SetFillColorAlpha(kBlue, 0.3); 
    g_band->SetLineColor(kBlue);
    g_band->SetLineWidth(2);
    g_band->SetMarkerStyle(20);

    g_band->GetXaxis()->SetTitle("m_{ee}[GeV/c^{2}]");
    g_band->GetYaxis()->SetTitle("#frac{1}{N_{events}}#frac{#epsilon_{BBC}#sigma_{pp}}{#epsilon_{bias}}#frac{dN}{dm_{ee}} [GeV/c^{2}]^{-1} [In PHENIX Acceptance]");

    TCanvas* c2 = new TCanvas();
    g_band->Draw("A3");  // A3 = axis + filled band
    g_band->Draw("LX SAME"); 

    TFile* fout = new TFile("getting_final_background_BG.root","RECREATE");
    fout->cd();
    h1d_mlp_dataBG->Write();
    h1d_like_combined->Write();
    g_band->Write();
    fout->Close();



}