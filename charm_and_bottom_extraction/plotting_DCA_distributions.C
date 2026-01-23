const int mass_regions = 13;
double mass_lo[mass_regions] = { 0.70, 0.87, 0.87, 1.30, 1.75, 2.20, 2.40, 1.60, 2.80, 3.30, 3.30, 2.40, 2.60};
double mass_hi[mass_regions] = { 1.15, 1.30, 1.60, 2.05, 2.40, 2.80, 2.95, 2.80, 3.15, 4.50, 5.50, 3.45, 4.00};

const char* mass_lo_name[mass_regions] = { "070", "087", "087", "130", "175", "220", "240", "160", "280", "330", "330", "240", "260"}; 
const char* mass_hi_name[mass_regions] = { "115", "130", "160", "205", "240", "280", "295", "280", "315", "450", "550", "345", "400"};

const char* region_name[mass_regions] = {"Region1", "Region2", "Region3", "Region4", "Region5", "Region6", "Region7", "Region8", "Region9", "Region10", "Region11", "Region12", "Region13"};

const int num_templates = 6;
const char* templates[num_templates] = {"dataFG", "dataBG", "prompt", "ccbar", "bbbar", "sum"};
const char* naming[num_templates] = {"FG_{e^{+}e^{-}}", "BG_{e^{+}e^{-}}", "prompt", "c#bar{c}", "b#bar{b}", "sum"};
const int template_color[num_templates] = {kBlack, kRed, kGreen-4, kViolet, kOrange, kRed};

TH1D* hist[num_templates] = {NULL};

TCanvas* c1[mass_regions] = {NULL};

void plotting_DCA_distributions(){

    TFile* f1 = new TFile("template_fits_prompt_bbbar_fixed_stat.root","READ");

    TFile* fout = new TFile("plotting_DCA_distributions_prompt_bbbar_fixed_stat.root","RECREATE");

    for(int iregion = 0; iregion < mass_regions; iregion++){

        c1[iregion] = new TCanvas(Form("DCA_%s" ,region_name[iregion]));
        c1[iregion]->cd();
        c1[iregion]->SetLogy();

        TLegend* leg = new TLegend(0.2923706,0.6873418,0.7717984,0.8873418);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05164557);
        leg->SetNColumns(2);

        for(int itemp = 0; itemp < num_templates; itemp++){

            if(itemp==0) hist[itemp] = (TH1D*)f1->Get(Form("DCA_%s_%s_original", templates[itemp], region_name[iregion]));
            else if(itemp>0 && itemp < num_templates-1) hist[itemp] = (TH1D*)f1->Get(Form("DCA_%s_%s", templates[itemp], region_name[iregion]));
            else hist[itemp] = (TH1D*)f1->Get(Form("hDCA_%s_%s", templates[itemp], region_name[iregion]));

            hist[itemp]->SetName(Form("DCA_%s_%s", templates[itemp], region_name[iregion]));

            if(itemp == 0) { 
                hist[itemp]->Draw(); hist[itemp]->SetMarkerStyle(20); hist[itemp]->SetMarkerColor(template_color[itemp]); 
                hist[itemp]->GetXaxis()->SetTitle("DCA_{T}[#mum]");
                hist[itemp]->GetYaxis()->SetTitle("#frac{dN}{d(DCA_{T})} [#mum]^{-1}");
                
            }

            else if(itemp == num_templates-1) { hist[itemp]->Draw("same hist"); hist[itemp]->SetLineStyle(1); hist[itemp]->SetLineColor(template_color[itemp]); }
            else{ hist[itemp]->Draw("same hist[itemp]"); hist[itemp]->SetLineStyle(2); hist[itemp]->SetLineColor(template_color[itemp]); }

            leg->AddEntry(hist[itemp], Form("%s", naming[itemp]));

        }

        leg->Draw("same");

        TLatex *   tex1 = new TLatex(500,0.01, Form("%.02f < m_{ee} < %.02f GeV/c^{2}", mass_lo[iregion], mass_hi[iregion]));
        tex1->SetTextAlign(13);
        tex1->SetTextSize(0.052);
        tex1->SetLineWidth(1);
        tex1->SetTextFont(42);
        tex1->Draw();

        c1[iregion]->Write();

    }

}