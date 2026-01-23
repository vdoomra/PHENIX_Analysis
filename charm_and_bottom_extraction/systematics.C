
const int npts = 7;
double x[npts] = {-999};
double y[npts] = {-999};
double y_err[npts] = {-999};

double ylo[npts] = {-999};
double ylo_err[npts] = {-999};

double yhi[npts] = {-999};
double yhi_err[npts] = {-999};

void systematics(){

    TFile* fdata_lower = new TFile("plotting_background_bbbar_lower.root","READ");
    TGraphErrors* graph_lower = (TGraphErrors*)fdata_lower->Get("g2");
    graph_lower->SetName("graph_lower");

    TFile* fdata_central = new TFile("plotting_background_bbbar_central.root","READ");
    TGraphErrors* graph_central = (TGraphErrors*)fdata_central->Get("g2");
    graph_central->SetName("graph_central");

    TFile* fdata_upper = new TFile("plotting_background_bbbar_upper.root","READ");
    TGraphErrors* graph_upper = (TGraphErrors*)fdata_upper->Get("g2");
    graph_upper->SetName("graph_upper");

    int nbins = graph_central->GetN();
    int counter = 0;

    for(int ibin=0; ibin < nbins; ibin++){

        double x_central, x_lower, x_upper;
        double y_central, y_lower, y_upper;

        graph_central->GetPoint(counter, x_central, y_central);
        graph_lower->GetPoint(counter, x_lower, y_lower);
        graph_upper->GetPoint(counter, x_upper, y_upper);

        double err_y_central = graph_central->GetErrorY(counter);
        double err_y_lower = graph_lower->GetErrorY(counter);
        double err_y_upper = graph_upper->GetErrorY(counter);

        x[counter] = x_central;

        ylo[counter] = fabs(y_central - y_lower)/y_central;
        yhi[counter] = fabs(y_central - y_upper)/y_central;

        y[counter] = (ylo[counter] + yhi[counter])/2;

        ylo_err[counter] = (y_lower/y_central)*sqrt( pow(err_y_central/y_central, 2)  + pow(err_y_lower/y_lower, 2) );

        yhi_err[counter] = (y_upper/y_central)*sqrt( pow(err_y_central/y_central, 2)  + pow(err_y_upper/y_upper, 2) );

        y_err[counter] = 0.5*sqrt(ylo_err[counter]*ylo_err[counter] + yhi_err[counter]*yhi_err[counter]); 

        counter++;
        

    }

    TGraphErrors* graph_lo = new TGraphErrors(counter, x, ylo, 0, ylo_err);
    graph_lo->SetName("graph_lo");

    TGraphErrors* graph_hi = new TGraphErrors(counter, x, yhi, 0, yhi_err);
    graph_hi->SetName("graph_hi");

    TGraphErrors* graph_average = new TGraphErrors(counter, x, y, 0, y_err);
    graph_average->SetName("graph_average");

    TGraph* graph_lo_no_err = new TGraph(counter, x, ylo);
    graph_lo_no_err->SetName("graph_lo_no_err");

    TGraph* graph_hi_no_err = new TGraph(counter, x, yhi);
    graph_hi_no_err->SetName("graph_hi_no_err");


    TFile* fout = new TFile("bbbar_DCA_Fits_systematics.root","RECREATE");
    fout->cd();
    graph_lo->Write();
    graph_average->Write();
    graph_lo_no_err->Write();
    graph_hi->Write();
    graph_hi_no_err->Write();
    fout->Close();

    
}