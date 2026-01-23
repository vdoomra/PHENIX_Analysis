#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TRandom.h>
#include "/direct/phenix+u/vdoomra/install/include/DileptonAnalysis/MyEvent.h"

using namespace std;
using namespace DileptonAnalysis;
const double pi = TMath::ACos(-1);

const double DCENTERCUT = 0.05; // RICH ghost cut
const double PC1_DPHI_CUT = 0.02, PC1_DZ_CUT = 0.5; // PC1 ghost cut 

const int nlayers = 4;
const int ncharge = 2;
const int nsides = 2;
const int total_nlayers = 8;
const int counter = 3;

static const double Me = 0.000510998918;
static const double Me2 = Me*Me;

const int nbins_mass = 40;
const double mass_range[nbins_mass+1] = { 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75, 0.80, 0.87, 0.95, 1.05, 1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.05, 3.10, 3.15, 3.20, 3.30, 3.45, 3.60, 3.70, 3.80, 4.00, 4.25, 4.5, 5.0, 5.5};

const int nbins_pt = 33;
const double pt_range[nbins_pt+1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.50, 6.00, 7.00, 8.00, 9.00, 10.00, 12.00, 14.00, 16.00, 18.00};

const int nbins_dca = 38;
const double dca_range[nbins_dca+1] = {0,40,80,120,160,200,240,280,320,360,400,440,480,520,560,600,640,680,720,760,800,840,880,920,960,1000,1040,1080,1120,1160,1200,1240,1280,1320,1360,1400,1440,1480,1520};

const int nbins_mass_prof = 10;
const double mass_range_prof[nbins_mass_prof+1] = { 1.35, 1.65, 1.95, 2.25, 2.55, 2.85, 3.05, 3.25, 3.65, 3.85, 4.2 };

const int MaxTracks = 5;
const int ntype = 3;

const char* charge_type[ncharge] = {"like", "unlike"};
const char* type_name[ntype] = {"total", "ccbar", "bbbar"};
const char* process_name[ntype] = {"FE", "PC", "GS"};

//Pair Mass vs Pair DCA

TH2D* hist_temp[ncharge][ntype][ntype]= {NULL};
TH2D* hist[ncharge][ntype][ntype] = {NULL};
TProfile* hist_prof[ncharge][ntype][ntype] = {NULL};

// Pair mass vs pair pT
TH2D* h2d_ee_FG_temp[ncharge][ntype][ntype]= {NULL};
TH2D* h2d_ee_FG[ncharge][ntype][ntype]= {NULL};
TH1D* h1d_ee_FG[ncharge][ntype][ntype]= {NULL};

TH1D* h1d_delta_phi_temp[ncharge][ntype][ntype]= {NULL};
TH1D* h1d_delta_phi[ncharge][ntype][ntype]= {NULL};

TF1* func = NULL;

TH1D* pt_spectra_single_electron_temp[ntype][ntype] = {NULL};
TH1D* pt_spectra_single_positron_temp[ntype][ntype] = {NULL};
TH1D* pt_spectra_temp[ntype][ntype] = {NULL};

TH1D* pt_spectra_single_electron[ntype][ntype] = {NULL};
TH1D* pt_spectra_single_positron[ntype][ntype] = {NULL};
TH1D* pt_spectra[ntype][ntype] = {NULL};

TH2D* dca_dist_single_electron_temp[ntype][ntype] = {NULL};
TH2D* dca_dist_single_electron[ntype][ntype] = {NULL};

TH2D* h2d_delta_phi_pair_mass_temp[ncharge][ntype][ntype] = {NULL};
TH2D* h2d_delta_phi_pair_mass[ncharge][ntype][ntype] = {NULL};

TH1D* h_debug = NULL;

double veto_window_mean(int layer, int charge, int delphi_side, double pt){

  double veto_window_mean_par0[nlayers-1][ncharge][nsides] = {-999};
  double veto_window_mean_par1[nlayers-1][ncharge][nsides] = {-999};
  double veto_window_mean_par2[nlayers-1][ncharge][nsides] = {-999};

  // Parameter 0
  veto_window_mean_par0[0][0][0] = -0.001876;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE          
  veto_window_mean_par0[0][0][1] =  0.003013;
  veto_window_mean_par0[0][1][0] = -0.003025;
  veto_window_mean_par0[0][1][1] =  0.001911;

  veto_window_mean_par0[1][0][0] = -0.004279;
  veto_window_mean_par0[1][0][1] =  0.006654;
  veto_window_mean_par0[1][1][0] = -0.006633;
  veto_window_mean_par0[1][1][1] =  0.004403;

  veto_window_mean_par0[2][0][0] = -0.006293;
  veto_window_mean_par0[2][0][1] =  0.009369;
  veto_window_mean_par0[2][1][0] = -0.01025;
  veto_window_mean_par0[2][1][1] =  0.00720;


  // Parameter 1
  veto_window_mean_par1[0][0][0] = -0.20232;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE          
  veto_window_mean_par1[0][0][1] = 0.06833;
  veto_window_mean_par1[0][1][0] = -0.10529;
  veto_window_mean_par1[0][1][1] = 0.03370;

  veto_window_mean_par1[1][0][0] = -0.0252;
  veto_window_mean_par1[1][0][1] = 0.0338;
  veto_window_mean_par1[1][1][0] = -0.03311;
  veto_window_mean_par1[1][1][1] = 0.02312;

  veto_window_mean_par1[2][0][0] = -0.02067;
  veto_window_mean_par1[2][0][1] = 0.05404;
  veto_window_mean_par1[2][1][0] = -0.05873;
  veto_window_mean_par1[2][1][1] = 0.02255;

  // Parameter 2
  veto_window_mean_par2[0][0][0] = -9.1789;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE                                          
  veto_window_mean_par2[0][0][1] = -8.58896;
  veto_window_mean_par2[0][1][0] = -9.4364;
  veto_window_mean_par2[0][1][1] = -6.8837;

  veto_window_mean_par2[1][0][0] = -2.1499;
  veto_window_mean_par2[1][0][1] = -1.9659;
  veto_window_mean_par2[1][1][0] = -1.8537;
  veto_window_mean_par2[1][1][1] = -1.7946;

  veto_window_mean_par2[2][0][0] = -1.5888;
  veto_window_mean_par2[2][0][1] = -1.7085;
  veto_window_mean_par2[2][1][0] = -1.8974;
  veto_window_mean_par2[2][1][1] = -1.9511;

  return veto_window_mean_par0[layer][charge][delphi_side] + veto_window_mean_par1[layer][charge][delphi_side]*exp(veto_window_mean_par2[layer][charge][delphi_side] * pt);
}

double veto_window_sigma(int layer, int charge, int delphi_side, double pt){

  double veto_window_sigma_par0[nlayers-1][ncharge][nsides] = {-999};
  double veto_window_sigma_par1[nlayers-1][ncharge][nsides] = {-999};
  double veto_window_sigma_par2[nlayers-1][ncharge][nsides] = {-999};

  // Parameter 0
  veto_window_sigma_par0[0][0][0] = 0.0002732;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE         
  veto_window_sigma_par0[0][0][1] = 0.0008276;
  veto_window_sigma_par0[0][1][0] = 0.0008218;
  veto_window_sigma_par0[0][1][1] = 0.0002933;

  veto_window_sigma_par0[1][0][0] = 0.00041;
  veto_window_sigma_par0[1][0][1] = 0.001677;
  veto_window_sigma_par0[1][1][0] = 0.001626;
  veto_window_sigma_par0[1][1][1] = 0.000314;

  veto_window_sigma_par0[2][0][0] = 0.001299;
  veto_window_sigma_par0[2][0][1] = 0.001982;
  veto_window_sigma_par0[2][1][0] = 0.002024;
  veto_window_sigma_par0[2][1][1] = 0.000925;


  // Parameter 1
  veto_window_sigma_par1[0][0][0] = 0.008219;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE          
  veto_window_sigma_par1[0][0][1] = 0.11697;
  veto_window_sigma_par1[0][1][0] = 0.05595;
  veto_window_sigma_par1[0][1][1] = 0.01959;

  veto_window_sigma_par1[1][0][0] = 0.00903;
  veto_window_sigma_par1[1][0][1] = 0.005199;
  veto_window_sigma_par1[1][1][0] = 0.004186;
  veto_window_sigma_par1[1][1][1] = 0.01498;
  veto_window_sigma_par1[2][0][0] = 0.01592;
  veto_window_sigma_par1[2][0][1] = 0.006797;
  veto_window_sigma_par1[2][1][0] = 0.007108;
  veto_window_sigma_par1[2][1][1] = 0.009907;

  // Parameter 2           
  veto_window_sigma_par2[0][0][0] = -6.0461;        // B1 LAYER, -VE CHARGE, -VE DELTA PHI SIDE  
  veto_window_sigma_par2[0][0][1] = -11.8997;
  veto_window_sigma_par2[0][1][0] = -10.3209;
  veto_window_sigma_par2[0][1][1] = -7.5985;

  veto_window_sigma_par2[1][0][0] = -2.3919;
  veto_window_sigma_par2[1][0][1] = -2.1199;
  veto_window_sigma_par2[1][1][0] = -1.6409;
  veto_window_sigma_par2[1][1][1] = -2.8119;

  veto_window_sigma_par2[2][0][0] = -3.5943;
  veto_window_sigma_par2[2][0][1] = -1.9312;
  veto_window_sigma_par2[2][1][0] = -1.9647;
  veto_window_sigma_par2[2][1][1] = -2.5570;

  return veto_window_sigma_par0[layer][charge][delphi_side] + veto_window_sigma_par1[layer][charge][delphi_side]*exp(veto_window_sigma_par2[layer][charge][delphi_side] * pt);
}


bool single_cut(MyTrack* mytrk)
{

  float px = mytrk->GetPtPrime()*TMath::Cos(mytrk->GetPhi0Prime());
  float py = mytrk->GetPtPrime()*TMath::Sin(mytrk->GetPhi0Prime());
  float pz = mytrk->GetPtPrime()*(TMath::Cos(mytrk->GetThe0Prime()))/(TMath::Sin(mytrk->GetThe0Prime()));
  float p = sqrt(px*px + py*py + pz*pz);
  float pt = sqrt( px*px + py*py );

  if(pt < 0.4) return false;
  if(mytrk->GetN0() < 2) return false;
  if((mytrk->GetEcore()/p) < 0.7) return false;

  int indexL1 = mytrk->GetHitIndexL1();
  int indexL2 = mytrk->GetHitIndexL2();
  int indexL3 = mytrk->GetHitIndexL3();
  int indexL4 = mytrk->GetHitIndexL4();
  int indexL5 = mytrk->GetHitIndexL5();
  int indexL6 = mytrk->GetHitIndexL6();
  int indexL7 = mytrk->GetHitIndexL7();
  int indexL8 = mytrk->GetHitIndexL8();

  if( !((indexL1 != -999) && (indexL2 != -999 || indexL3 != -999 || indexL4 != -999 || indexL5 != -999 || indexL6 != -999 || indexL7 != -999 || indexL8 != -999 )) ) return false;

  return true;
}

double redefine_phi0(double vx, double vy, double x1, double y1, double x2, double y2, double err){


  double x[counter] = {vx, x1, x2};
  double y[counter] = {vy, y1, y2};

  double err_y[counter] = {0.0117, 0.0050, err};

  TGraphErrors* graph = new TGraphErrors(counter, x, y, 0, err_y);

  if(x[1] > 0) func = new TF1("func","pol2", -1.0, 20.0);
  else func = new TF1("func","pol2", -20.0, 1.0);
  graph->Fit(func,"QR");

  double par1 = func->GetParameter(1);
  double par2 = func->GetParameter(2);

  double slope = par1 + 2*par2*x[0];

  double phi0_new_method = TMath::ATan(slope);

  if((x[1] < 0 && y[1] > 0) || (x[1]<0 && y[1]<0)) phi0_new_method += pi;

  return phi0_new_method;

}

void DCA_ccbar_SoftQCD(const char* inFile, const char* outFile, const int is_ert, const int additional_rejection, const int gen_cut, const char* simfile)
{

  h_debug = new TH1D("h_debug","Debug Histogram", 5, -0.5, 4.5);

  for(int i = 0; i < ntype; i++){

    for(int k = 0; k < ntype; k++){

      for(int j=0; j<ncharge; j++){

        hist_temp[j][i][k] = new TH2D(Form("hist_%s_%s_%s_temp", charge_type[j], type_name[i], process_name[k]),"DCA Distribution", nbins_mass, mass_range, nbins_dca, dca_range);
        hist_temp[j][i][k]->Sumw2();

        hist_prof[j][i][k] = new TProfile(Form("hist_prof_%s_%s_%s", charge_type[j], type_name[i], process_name[k]),"DCA Distribution Profile Histogram", nbins_mass_prof, mass_range_prof);
        hist_prof[j][i][k]->Sumw2();

        h2d_ee_FG_temp[j][i][k] = new TH2D(Form("h2d_ee_FG_%s_%s_%s_temp", charge_type[j], type_name[i], process_name[k]), "ee mass distribution FG", nbins_mass, mass_range ,100, 0, 10);
        h2d_ee_FG_temp[j][i][k]->Sumw2();

        h1d_delta_phi_temp[j][i][k] = new TH1D(Form("h1d_delta_phi_%s_%s_%s_temp",charge_type[j], type_name[i], process_name[k]), "#Delta Phi Distribution", 10, 0, pi);
        h1d_delta_phi_temp[j][i][k]->Sumw2();

        h2d_delta_phi_pair_mass_temp[j][i][k] = new TH2D(Form("h2d_delta_phi_pair_mass_%s_%s_%s_temp", charge_type[j], type_name[i], process_name[k]), "#Delta #phi vs Pair Mass Distribution", nbins_mass, mass_range, 10, 0, TMath::Pi());
        h2d_delta_phi_pair_mass_temp[j][i][k]->Sumw2();

      }

      pt_spectra_single_electron_temp[i][k] = new TH1D(Form("pt_spectra_single_electron_%s_%s_temp", type_name[i], process_name[k]), "pT spectra", nbins_pt, pt_range);
      pt_spectra_single_electron_temp[i][k]->Sumw2();

      pt_spectra_single_positron_temp[i][k] = new TH1D(Form("pt_spectra_single_positron_%s_%s_temp", type_name[i], process_name[k]), "pT spectra", nbins_pt, pt_range);
      pt_spectra_single_positron_temp[i][k]->Sumw2();

      pt_spectra_temp[i][k] = new TH1D(Form("pt_spectra_%s_%s_temp", type_name[i], process_name[k]), "pT spectra", nbins_pt, pt_range);
      pt_spectra_temp[i][k]->Sumw2();

      dca_dist_single_electron_temp[i][k] = new TH2D(Form("dca_dist_single_electron_%s_%s_temp", type_name[i], process_name[k]), "DCA_distribution", nbins_pt, pt_range, 100, -1600, 1600);
      dca_dist_single_electron_temp[i][k]->Sumw2();


    }
  }

  TFile* input = new TFile(inFile,"READ");
  if(!(input))
  {
    cout << "no input file" << endl;
    exit(1);
  }

  TTree* T = (TTree*)input->Get("T");
  TBranch* br = T->GetBranch("MyEvent");
  MyEvent* event = 0;
  br->SetAddress(&event);

  cout << "Trees read!" << endl; 

  TFile* fsim = new TFile(simfile,"READ");
  TH1D* hnorm = (TH1D*)fsim->Get("hnorm");
  TTree* tree = (TTree*)fsim->Get("T");
  int ntracks = -999;
  int pid[MaxTracks] = {-999};
  double px[MaxTracks] = {-999};
  double py[MaxTracks] = {-999};
  double weight[MaxTracks] = {-999};
  int isCharm[MaxTracks] = {-999};
  int isBottom[MaxTracks] = {-999};
  int isFE[MaxTracks] = {-999};
  int isPC[MaxTracks] = {-999};
  int isGS[MaxTracks] = {-999};

  tree->SetBranchAddress("ntracks", &ntracks);
  tree->SetBranchAddress("pid",pid);
  tree->SetBranchAddress("px",px);
  tree->SetBranchAddress("py",py);
  tree->SetBranchAddress("weight",weight);
  tree->SetBranchAddress("isCharm",isCharm);
  tree->SetBranchAddress("isBottom",isBottom);
  tree->SetBranchAddress("isFE",isFE);
  tree->SetBranchAddress("isPC",isPC);
  tree->SetBranchAddress("isGS",isGS);

  double normalization_factor = hnorm->GetBinContent(hnorm->GetXaxis()->FindBin(0.0));

  TFile* frun_weights = new TFile("/phenix/plhf/vdoomra/DCA_analysis/runbyrun_weighting_output.root","READ");
  TH1D* hrunbyrun_weights = (TH1D*)frun_weights->Get("hist");
  hrunbyrun_weights->SetName("hrunbyrun_weights");

  TFile* output = new TFile(outFile,"RECREATE");
 
  int nevt = T->GetEntries();

  TF1* ftrigger = new TF1("ftrigger","[3] + ([0]-[3])/(1 + (x/[2])^[1])", 0, 10);

  int accepted_events = 0;
  int analysis_events = 0;

  double run_weight = 1.0;

  for (int ievent = 0; ievent < nevt; ievent++)
  {
    if (ievent%50000==0) cout << "Event: " << ievent << " / " << nevt << endl;
    gRandom->SetSeed(0);

    event->ClearEvent();
    br->GetEntry(ievent);
    
    int ntrack = event->GetNtrack();
    int nvtxhits = event->GetNVTXhit();
    float vtxz = event->GetPreciseZ();
    int run_number = event->GetRunNumber();
    run_weight = hrunbyrun_weights->GetBinContent(hrunbyrun_weights->FindBin(run_number));
    int evtno = event->GetEvtNo();

    if(fabs(vtxz) > 10) continue;

    accepted_events++;

    int ghost[40], real[40];
    int all[40];
    int nghost = 0, nreal = 0, nall = 0;
    
    double vx = event->GetPreciseX();
    double vy = event->GetPreciseY();

    for (int k1 = 0; k1 < ntrack; ++k1)
    {
      
      all[nall] = k1;
      ++nall;        
      int ghost_flag = 0;

      for (int k2 = 0; k2 < ntrack; ++k2)
      {
        if ( k1==k2 ) continue;
        if ( fabs((event->GetEntry(k1)).GetCrkphi() - (event->GetEntry(k2)).GetCrkphi() ) < DCENTERCUT ) { ghost_flag = 1; break; }
      }

      if ( ghost_flag ) { ghost[nghost] = k1; ++nghost; }
      else { real[nreal] = k1; ++nreal; }
    }

    if(nreal<1) continue;

    int nprimary = 0;
    int primary_id[10] = {-999};

    if(additional_rejection){

      for(int ireal=0; ireal<nreal; ireal++){

        MyTrack mytrk = event->GetEntry(real[ireal]);
        double px = mytrk.GetPtPrime()*TMath::Cos(mytrk.GetPhi0Prime());
        double py = mytrk.GetPtPrime()*TMath::Sin(mytrk.GetPhi0Prime());
        double pt = sqrt( px*px + py*py );

        int charge = mytrk.GetChargePrime();

        int charge_index = -999;
        if(charge<0) charge_index = 0;
        else charge_index = 1;

        int hit_index[total_nlayers] = {-999};

        hit_index[0] = mytrk.GetHitIndexL1();
        hit_index[1] = mytrk.GetHitIndexL2();
        hit_index[2] = mytrk.GetHitIndexL3();
        hit_index[3] = mytrk.GetHitIndexL4();
        hit_index[4] = mytrk.GetHitIndexL5();
        hit_index[5] = mytrk.GetHitIndexL6();
        hit_index[6] = mytrk.GetHitIndexL7();
        hit_index[7] = mytrk.GetHitIndexL8();

	      bool veto_L1 = false; bool veto_L2 = false; bool veto_L3 = false;

        for(int ilayer=1; ilayer<total_nlayers; ilayer++){

          if(hit_index[ilayer] == -999) continue;

          int layer_index = -999;
          if(ilayer == 1) layer_index = 0;
          else if(ilayer == 2 || ilayer == 3 || ilayer == 4) layer_index = 1;
          else layer_index = 2;

          MyVTXHit vtxhit = event->GetVTXHitEntry(hit_index[ilayer]);
          double zhit = vtxhit.GetZHit();
          double phi_hit = TMath::ATan2(vtxhit.GetYHit(), vtxhit.GetXHit());
          if(phi_hit < -pi/2) phi_hit += 2*pi;

	        for(int ihit = 0; ihit < nvtxhits; ihit++){

            if(ihit == hit_index[ilayer]) continue;
            int dphi_index = -999;
            MyVTXHit vtx_near_hit = event->GetVTXHitEntry(ihit);
            int layer = vtx_near_hit.GetLayer();
            if( (ilayer==1) && ilayer != layer) continue;
            if( (ilayer==2 || ilayer==3 || ilayer==4) && layer !=2) continue;
            if( (ilayer==5 || ilayer==6 || ilayer==7) && layer !=3) continue;
            double phi_near_hit = TMath::ATan2(vtx_near_hit.GetYHit(), vtx_near_hit.GetXHit());
            if(phi_near_hit < -pi/2) phi_near_hit += 2*pi;
            double z_near_hit = vtx_near_hit.GetZHit();

            double dphi = phi_near_hit - phi_hit;
            double dz = z_near_hit - zhit;

	          if(layer==1 && fabs(dz) > 0.2) continue;
            else if((layer==2 || layer==3) && fabs(dz) > 0.4 ) continue;

            if(dphi<0) dphi_index = 0;
            else dphi_index = 1;

            double mean = veto_window_mean(layer_index, charge_index, dphi_index, pt);
            double sigma = veto_window_sigma(layer_index, charge_index, dphi_index, pt);

            if( fabs(dphi) > 0.001 && fabs(dphi) < fabs(mean)+15*sigma && layer==1) veto_L1 = true;
            if( fabs(dphi) > 0.001 && fabs(dphi) < fabs(mean)+15*sigma && layer==2) veto_L2 = true;
            if( fabs(dphi) > 0.001 && fabs(dphi) < fabs(mean)+15*sigma && layer==3) veto_L3 = true;

          }
        }

	      if( (veto_L1 || veto_L2 || veto_L3) ) continue;
        primary_id[nprimary] = real[ireal];
        nprimary++;

      }

    }

    else{

      nprimary = nreal;
      for(int ii=0; ii<nreal; ii++) { primary_id[ii] = real[ii];  }

    }

    for (int iprimary_A1 = 0; iprimary_A1 < nprimary; iprimary_A1++)
    {
      MyTrack mytrk_A1 = event->GetEntry(primary_id[iprimary_A1]);
      if(!(single_cut(&mytrk_A1))) continue;

      int hit_index_A1[total_nlayers] = {-999};

      hit_index_A1[0] = mytrk_A1.GetHitIndexL1();
      hit_index_A1[1] = mytrk_A1.GetHitIndexL2();
      hit_index_A1[2] = mytrk_A1.GetHitIndexL3();
      hit_index_A1[3] = mytrk_A1.GetHitIndexL4();
      hit_index_A1[4] = mytrk_A1.GetHitIndexL5();
      hit_index_A1[5] = mytrk_A1.GetHitIndexL6();
      hit_index_A1[6] = mytrk_A1.GetHitIndexL7();
      hit_index_A1[7] = mytrk_A1.GetHitIndexL8();

      double x1_A1 = -999; double y1_A1 = -999;
      double x2_A1 = -999; double y2_A1 = -999;
      int layer2_A1 = -999;

      for(int ilayer=0; ilayer<total_nlayers; ilayer++){

        if(hit_index_A1[ilayer] == -999) continue;

        MyVTXHit vtxhit = event->GetVTXHitEntry(hit_index_A1[ilayer]);
        if(x1_A1 == -999){ x1_A1 = vtxhit.GetXHit();  y1_A1 = vtxhit.GetYHit(); continue; }
        if(x2_A1 == -999){ x2_A1 = vtxhit.GetXHit();  y2_A1 = vtxhit.GetYHit(); layer2_A1 = vtxhit.GetLayer(); }

      }

      double err_A1 = -999;
      if( layer2_A1 == 1) err_A1 = 0.0050;
      else err_A1 = 0.0080;

      double phi0_A1 = redefine_phi0(0.0, 0.0, x1_A1, y1_A1, x2_A1, y2_A1, err_A1);

      double px_A1 = mytrk_A1.GetPtPrime()*TMath::Cos(phi0_A1);
      double py_A1 = mytrk_A1.GetPtPrime()*TMath::Sin(phi0_A1);
      double pz_A1 = mytrk_A1.GetPtPrime()*(TMath::Cos(mytrk_A1.GetThe0Prime()))/(TMath::Sin(mytrk_A1.GetThe0Prime()));

      double pt_A1 = sqrt(px_A1*px_A1 + py_A1*py_A1);
      double R_A1 = pt_A1/(0.003*0.90); // Radius in cm

      double phiDC_A1 = mytrk_A1.GetPhiDC();
      double zed_A1 = mytrk_A1.GetZDC();

      int charge_A1 = mytrk_A1.GetChargePrime();
      int mcid_A1 = mytrk_A1.GetMcId();
      bool found_A1 = false;
      bool is_charm_A1 = false;
      bool is_bottom_A1 = false;
      double true_pt_A1 = -999;
      double weight_A1 = 1.0;

      bool flavor_excitation_A1 = false;
      bool pair_creation_A1 = false;
      bool gluon_splitting_A1 = false;

      for(int k=1; k<(int)event->GetNMcTrack(); k++){

        McTrack mctrk = event->GetMcEntry(k);
        if(mctrk.GetMcIndex() != mcid_A1) continue;
        int gen = mctrk.GetGen();
        if(gen == 1) found_A1 = true;

	if(found_A1) {  
		
	  true_pt_A1 = sqrt(mctrk.GetPx()*mctrk.GetPx() + mctrk.GetPy()*mctrk.GetPy());
	  int pid_A1 = mctrk.GetParticleID();

	  tree->GetEntry(evtno-1);
	  
	  for(int jj=0 ; jj<(int)ntracks; jj++){

	    if( pid_A1 == 2 && pid[jj] != -11 ) continue;
	    if( pid_A1 == 3 && pid[jj] != 11 ) continue;
            if(weight[jj] < 0) weight[jj] = 1.0;
	    weight_A1 = weight[jj];
            if(isCharm[jj] > 0) is_charm_A1 = true;
            else is_bottom_A1 = true;

            if(isFE[jj] > 0) flavor_excitation_A1 = true;
            if(isPC[jj] > 0) pair_creation_A1 = true;
            if(isGS[jj] > 0) gluon_splitting_A1 = true;
	    
	  }
	  
	}

      }

      if(gen_cut){ if(!found_A1) continue; }

      TLorentzVector trk_A1;
      trk_A1.SetPxPyPzE(px_A1, py_A1, pz_A1, sqrt(px_A1*px_A1 + py_A1*py_A1 + pz_A1*pz_A1 + Me2));

      TVector3 trk3_A1;
      trk3_A1.SetXYZ(px_A1, py_A1, pz_A1);

      bool isERT_A1 = false;
      double rand_A1 = gRandom->Uniform(0,1);
      double ecore_A1 = sqrt(pt_A1*pt_A1 + pz_A1*pz_A1 + Me2);
      double energy_threshold_A1 = -999;

      if(mytrk_A1.GetArm() == 1){  // West Arm

	if (mytrk_A1.GetSect() == 0) { ftrigger->SetParameters(0.0, 7.74596, 1.80702, 1.0); energy_threshold_A1 = 1.652; }
        else if (mytrk_A1.GetSect() == 1) { ftrigger->SetParameters(0.0, 8.19638, 1.64472, 1.0); energy_threshold_A1 = 1.511; }
        else if (mytrk_A1.GetSect() == 2) { ftrigger->SetParameters(0.0, 7.56855, 1.68859, 1.0); energy_threshold_A1 = 1.541; }
        else { ftrigger->SetParameters(0.0, 10.7913, 1.59629, 1.0); energy_threshold_A1 = 1.497; }

      }


      else {  // East Arm                                                                                                                               

        if (mytrk_A1.GetSect() == 0) { ftrigger->SetParameters(0.0, 6.07937, 2.50881, 1.0); energy_threshold_A1 = 2.239; }
        else if (mytrk_A1.GetSect() == 1) { ftrigger->SetParameters(0.0, 6.87068, 2.83264, 1.0); energy_threshold_A1 = 2.561; }
        else if(mytrk_A1.GetSect() == 2) { ftrigger->SetParameters(0.0, 8.078, 1.82598, 1.0); energy_threshold_A1 = 1.676; }
        else { ftrigger->SetParameters(0.0, 7.05411, 1.61651, 1.0); energy_threshold_A1 = 1.465; }

      }

      if(rand_A1 < ftrigger->Eval(mytrk_A1.GetEcore()) && ecore_A1 > energy_threshold_A1) isERT_A1 = true;

      int heavy_flavor_index_A1 = -999;
      int process_index_A1 = -999;

      if(flavor_excitation_A1) process_index_A1 = 0;
      else if(pair_creation_A1) process_index_A1 = 1;
      else process_index_A1 = 2;

      if(is_charm_A1) heavy_flavor_index_A1 = 1;
      else if(is_bottom_A1) heavy_flavor_index_A1 = 2;

      if(is_ert){

        if(!isERT_A1) continue;

        pt_spectra_temp[0][process_index_A1]->Fill(pt_A1, run_weight*weight_A1);
        pt_spectra_temp[heavy_flavor_index_A1][process_index_A1]->Fill(pt_A1, run_weight*weight_A1);

        if(charge_A1 == -1) {
          
          pt_spectra_single_electron_temp[0][process_index_A1]->Fill(pt_A1, run_weight*weight_A1);
          pt_spectra_single_electron_temp[heavy_flavor_index_A1][process_index_A1]->Fill(pt_A1, run_weight*weight_A1); 
          
        }
  
        else { 
            
          pt_spectra_single_positron_temp[0][process_index_A1]->Fill(pt_A1, run_weight*weight_A1);
          pt_spectra_single_positron_temp[heavy_flavor_index_A1][process_index_A1]->Fill(pt_A1, run_weight*weight_A1);
  
        }

      }


      else{

        pt_spectra_temp[0][process_index_A1]->Fill(pt_A1, run_weight*weight_A1);
        pt_spectra_temp[heavy_flavor_index_A1][process_index_A1]->Fill(pt_A1, run_weight*weight_A1);

        if(charge_A1 == -1) {
          
          pt_spectra_single_electron_temp[0][process_index_A1]->Fill(pt_A1, run_weight*weight_A1);
          pt_spectra_single_electron_temp[heavy_flavor_index_A1][process_index_A1]->Fill(pt_A1, run_weight*weight_A1); 
          
        }
  
        else { 
            
          pt_spectra_single_positron_temp[0][process_index_A1]->Fill(pt_A1, run_weight*weight_A1);
          pt_spectra_single_positron_temp[heavy_flavor_index_A1][process_index_A1]->Fill(pt_A1, run_weight*weight_A1);
  
        }

      }

      double a_A1 = ((x1_A1-x2_A1) -(y1_A1*y1_A1 - y2_A1*y2_A1)/(x1_A1-x2_A1))*0.5;
      double b_A1 = (y1_A1-y2_A1)/(x1_A1-x2_A1);

      double yc1_A1 = (y1_A1-a_A1*b_A1 + sqrt(R_A1*R_A1*(1+b_A1*b_A1) - (a_A1 + b_A1*y1_A1)*(a_A1 + b_A1*y1_A1)))/(1+b_A1*b_A1);
      double yc2_A1 = (y1_A1-a_A1*b_A1 - sqrt(R_A1*R_A1*(1+b_A1*b_A1) - (a_A1 + b_A1*y1_A1)*(a_A1 + b_A1*y1_A1)))/(1+b_A1*b_A1);

      double xc1_A1 = x1_A1 - a_A1 - b_A1*yc1_A1;
      double xc2_A1 = x1_A1 - a_A1 - b_A1*yc2_A1;

      double xc_A1 = -999; double yc_A1 = -999;

      if(phi0_A1 > -pi/2 && phi0_A1 < pi/2 && charge_A1 == -1){ 
	
	      if(yc1_A1<0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
	      else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      if(phi0_A1 > -pi/2 && phi0_A1 < pi/2 && charge_A1 == 1){

              if(yc1_A1>0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
	      else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      if(phi0_A1 > pi/2 && phi0_A1 < 3*pi/2 && charge_A1 == -1){

              if(yc1_A1>0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
	      else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      if(phi0_A1 > pi/2 && phi0_A1 < 3*pi/2 && charge_A1 == 1){

        if(yc1_A1<0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
        else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      double L_A1 = sqrt((xc_A1-vx)*(xc_A1-vx) + (yc_A1-vy)*(yc_A1-vy));
      double dca_A1 = (L_A1-R_A1)*10000; // In Micro Meters

      if(is_ert) { 

        if(!isERT_A1) continue;
  
        dca_dist_single_electron_temp[0][process_index_A1]->Fill(pt_A1, dca_A1, run_weight*weight_A1); 
        dca_dist_single_electron_temp[heavy_flavor_index_A1][process_index_A1]->Fill(pt_A1, dca_A1, run_weight*weight_A1);  
      
      }

      else { 
        
        dca_dist_single_electron_temp[0][process_index_A1]->Fill(pt_A1, dca_A1, run_weight*weight_A1); 
        dca_dist_single_electron_temp[heavy_flavor_index_A1][process_index_A1]->Fill(pt_A1, dca_A1, run_weight*weight_A1); 
      
      }

      for (int iprimary_A2 = iprimary_A1 +1; iprimary_A2 < nprimary; iprimary_A2++)
      {
        MyTrack mytrk_A2 = event->GetEntry(primary_id[iprimary_A2]);
        if(!(single_cut(&mytrk_A2))) continue;

	      int hit_index_A2[total_nlayers] = {-999};

        hit_index_A2[0] = mytrk_A2.GetHitIndexL1();
        hit_index_A2[1] = mytrk_A2.GetHitIndexL2();
        hit_index_A2[2] = mytrk_A2.GetHitIndexL3();
        hit_index_A2[3] = mytrk_A2.GetHitIndexL4();
        hit_index_A2[4] = mytrk_A2.GetHitIndexL5();
        hit_index_A2[5] = mytrk_A2.GetHitIndexL6();
        hit_index_A2[6] = mytrk_A2.GetHitIndexL7();
        hit_index_A2[7] = mytrk_A2.GetHitIndexL8();

        double x1_A2 = -999; double y1_A2 = -999;
        double x2_A2 = -999; double y2_A2 = -999;
        int layer2_A2 = -999;

        for(int ilayer=0; ilayer<total_nlayers; ilayer++){

          if(hit_index_A2[ilayer] == -999) continue;

          MyVTXHit vtxhit = event->GetVTXHitEntry(hit_index_A2[ilayer]);
          if(x1_A2 == -999){ x1_A2 = vtxhit.GetXHit();  y1_A2 = vtxhit.GetYHit(); continue; }
          if(x2_A2 == -999){ x2_A2 = vtxhit.GetXHit();  y2_A2 = vtxhit.GetYHit(); layer2_A2 = vtxhit.GetLayer(); }

        }

	      double err_A2 = -999;
        if( layer2_A2 == 1) err_A2 = 0.0050;
        else err_A2 = 0.0080;

        double phi0_A2 = redefine_phi0(0.0, 0.0, x1_A2, y1_A2, x2_A2, y2_A2, err_A2);

        double px_A2 = mytrk_A2.GetPtPrime()*TMath::Cos(phi0_A2);
        double py_A2 = mytrk_A2.GetPtPrime()*TMath::Sin(phi0_A2);
        double pz_A2 = mytrk_A2.GetPtPrime()*(TMath::Cos(mytrk_A2.GetThe0Prime()))/(TMath::Sin(mytrk_A2.GetThe0Prime()));
	      double pt_A2 = sqrt(px_A2*px_A2 + py_A2*py_A2);
        double R_A2 = pt_A2/(0.003*0.90); // Radius in cm

	      int charge_A2 = mytrk_A2.GetChargePrime();
	      int mcid_A2 = mytrk_A2.GetMcId();
        bool found_A2 = false;

        double phiDC_A2 = mytrk_A2.GetPhiDC();
	      double zed_A2 = mytrk_A2.GetZDC();
	      double true_pt_A2 = -999;
	      double weight_A2 = 1.0;
        bool is_charm_A2 = false;
        bool is_bottom_A2 = false;
        bool flavor_excitation_A2 = false;
        bool pair_creation_A2 = false;
        bool gluon_splitting_A2 = false;

        for(int k=1; k<(int)event->GetNMcTrack(); k++){

          McTrack mctrk = event->GetMcEntry(k);
          if(mctrk.GetMcIndex() != mcid_A2) continue;
          int gen = mctrk.GetGen();
          if(gen == 1) found_A2 = true;
	  
	        if(found_A2) {

            true_pt_A2 = sqrt(mctrk.GetPx()*mctrk.GetPx() + mctrk.GetPy()*mctrk.GetPy());
            int pid_A2 = mctrk.GetParticleID();
        
            tree->GetEntry(evtno-1);

            for(int jj=0; jj<(int)ntracks; jj++){

              if( pid_A2 == 2 && pid[jj] != -11 ) continue;
              if( pid_A2 == 3 && pid[jj] != 11 ) continue;
              if(weight[jj] < 0) weight[jj] = 1.0;
              weight_A2 = weight[jj];
              if(isCharm[jj] > 0) is_charm_A2 = true;
              else is_bottom_A2 = true;

              if(isFE[jj] > 0) flavor_excitation_A2 = true;
              if(isPC[jj] > 0) pair_creation_A2 = true;
              if(isGS[jj] > 0) gluon_splitting_A2 = true;

            }

          }

        }

        if(gen_cut){ if (!found_A2) continue; }

        TLorentzVector trk_A2;
        trk_A2.SetPxPyPzE(px_A2, py_A2, pz_A2, sqrt(px_A2*px_A2 + py_A2*py_A2 + pz_A2*pz_A2 + Me2));

        TVector3 trk3_A2;
        trk3_A2.SetXYZ(px_A2, py_A2, pz_A2);

        bool isERT_A2 = false;

        double rand_A2 = gRandom->Uniform(0,1);
	      double ecore_A2 = sqrt(pt_A2*pt_A2 + pz_A2*pz_A2 + Me2);
	      double energy_threshold_A2 = -999;

        if(mytrk_A2.GetArm() == 1){  // West Arm                                                                                                                            

          if (mytrk_A2.GetSect() == 0) { ftrigger->SetParameters(0.0, 7.74596, 1.80702, 1.0); energy_threshold_A2 = 1.652; }
          else if (mytrk_A2.GetSect() == 1) { ftrigger->SetParameters(0.0, 8.19638, 1.64472, 1.0); energy_threshold_A2 = 1.511; }
          else if (mytrk_A2.GetSect() == 2) { ftrigger->SetParameters(0.0, 7.56855, 1.68859, 1.0); energy_threshold_A2 = 1.541; }
          else { ftrigger->SetParameters(0.0, 10.7913, 1.59629, 1.0); energy_threshold_A2 = 1.497; }

        }


        else {  // East Arm                                                                                                                                                 

          if (mytrk_A2.GetSect() == 0) { ftrigger->SetParameters(0.0, 6.07937, 2.50881, 1.0); energy_threshold_A2 = 2.239;}
          else if (mytrk_A2.GetSect() == 1) { ftrigger->SetParameters(0.0, 6.87068, 2.83264, 1.0); energy_threshold_A2 = 2.561; }
          else if(mytrk_A2.GetSect() == 2) { ftrigger->SetParameters(0.0, 8.078, 1.82598, 1.0); energy_threshold_A2 = 1.676; }
          else { ftrigger->SetParameters(0.0, 7.05411, 1.61651, 1.0); energy_threshold_A2 = 1.465; }

        }

	if(rand_A2 < ftrigger->Eval(mytrk_A2.GetEcore()) && ecore_A2 > energy_threshold_A2) isERT_A2 = true;
	      
	if ( is_ert ) {  if(!( isERT_A1 || isERT_A2 )){  continue; } }
	if ( fabs(phiDC_A1-phiDC_A2)<PC1_DPHI_CUT && fabs(zed_A1-zed_A2)<PC1_DZ_CUT) continue;
	      
        double a_A2 = ((x1_A2-x2_A2) -(y1_A2*y1_A2 - y2_A2*y2_A2)/(x1_A2-x2_A2))*0.5;
        double b_A2 = (y1_A2-y2_A2)/(x1_A2-x2_A2);

        double yc1_A2 = (y1_A2-a_A2*b_A2 + sqrt(R_A2*R_A2*(1+b_A2*b_A2) - (a_A2 + b_A2*y1_A2)*(a_A2 + b_A2*y1_A2)))/(1+b_A2*b_A2);
        double yc2_A2 = (y1_A2-a_A2*b_A2 - sqrt(R_A2*R_A2*(1+b_A2*b_A2) - (a_A2 + b_A2*y1_A2)*(a_A2 + b_A2*y1_A2)))/(1+b_A2*b_A2);

        double xc1_A2 = x1_A2 - a_A2 - b_A2*yc1_A2;
        double xc2_A2 = x1_A2 - a_A2 - b_A2*yc2_A2;

        double xc_A2 = -999; double yc_A2 = -999;

        if(phi0_A2 > -pi/2 && phi0_A2 < pi/2 && charge_A2 == -1){ 
	
	        if(yc1_A2<0){ yc_A2 = yc1_A2; xc_A2 = xc1_A2; }
	        else{ yc_A2 = yc2_A2; xc_A2 = xc2_A2; }

        }

        if(phi0_A2 > -pi/2 && phi0_A2 < pi/2 && charge_A2 == 1){

          if(yc1_A2>0){ yc_A2 = yc1_A2; xc_A2 = xc1_A2; }
	        else{ yc_A2 = yc2_A2; xc_A2 = xc2_A2; }

        }

        if(phi0_A2 > pi/2 && phi0_A2 < 3*pi/2 && charge_A2 == -1){

          if(yc1_A2>0){ yc_A2 = yc1_A2; xc_A2 = xc1_A2; }
	        else{ yc_A2 = yc2_A2; xc_A2 = xc2_A2; }

        }

        if(phi0_A2 > pi/2 && phi0_A2 < 3*pi/2 && charge_A2 == 1){

          if(yc1_A2<0){ yc_A2 = yc1_A2; xc_A2 = xc1_A2; }
          else{ yc_A2 = yc2_A2; xc_A2 = xc2_A2; }

        }

        double L_A2 = sqrt((xc_A2-vx)*(xc_A2-vx) + (yc_A2-vy)*(yc_A2-vy));
        double dca_A2 = (L_A2-R_A2)*10000; // In Micro Meters

        TLorentzVector pair = trk_A1 + trk_A2;
        double pair_mass = pair.M();
        double pair_pt = pair.Pt();

        double dca1 = sqrt( fabs(dca_A1*dca_A1 - dca_A2*dca_A2) );

        double angle = trk3_A1.Angle(trk3_A2);

        if(dca1 > 1520) continue;
        if(pair_pt > 5.0) continue;
	
	      analysis_events++;

        int process_index = -999;
        int q_index = -999;
        int heavy_flavor_index = -999;

	if(charge_A1 == charge_A2) q_index = 0;
        else q_index = 1;

        if(flavor_excitation_A1 && flavor_excitation_A2) process_index = 0;
        else if(pair_creation_A1 && pair_creation_A2) process_index = 1;
        else if(gluon_splitting_A1 && gluon_splitting_A2) process_index = 2;

        h2d_ee_FG_temp[q_index][0][process_index]->Fill(pair_mass, pair_pt, run_weight*weight_A1*weight_A2);
        hist_temp[q_index][0][process_index]->Fill(pair_mass, dca1, run_weight*weight_A1*weight_A2);
        hist_prof[q_index][0][process_index]->Fill(pair_mass, dca1, run_weight*weight_A1*weight_A2);	
        h1d_delta_phi_temp[q_index][0][process_index]->Fill(angle, run_weight*weight_A1*weight_A2);
        h2d_delta_phi_pair_mass_temp[q_index][0][process_index]->Fill(pair_mass, angle, run_weight*weight_A1*weight_A2);

        if(is_charm_A1 && is_charm_A2) heavy_flavor_index = 1;
        else if(is_bottom_A1 && is_bottom_A2) heavy_flavor_index = 2;

        h2d_ee_FG_temp[q_index][heavy_flavor_index][process_index]->Fill(pair_mass, pair_pt, run_weight*weight_A1*weight_A2);
        hist_temp[q_index][heavy_flavor_index][process_index]->Fill(pair_mass, dca1, run_weight*weight_A1*weight_A2);
        hist_prof[q_index][heavy_flavor_index][process_index]->Fill(pair_mass, dca1, run_weight*weight_A1*weight_A2);	
        h1d_delta_phi_temp[q_index][heavy_flavor_index][process_index]->Fill(angle, run_weight*weight_A1*weight_A2);
        h2d_delta_phi_pair_mass_temp[q_index][heavy_flavor_index][process_index]->Fill(pair_mass, angle, run_weight*weight_A1*weight_A2);
	      	
      } // End Track A2 Loop
    
    } // End Track A1 Loop
  } // End Event A loop

  cout << analysis_events << "\t" << accepted_events << "\t" << nevt << endl;

  h_debug->Fill(0.0);
  if(analysis_events > 0) h_debug->Fill(1.0);

  normalization_factor *= run_weight;


  for(int itype = 0; itype < ntype; itype++){

    for(int icharge = 0; icharge < ncharge; icharge++){

      for(int iprocess = 0; iprocess < ntype; iprocess++){

        h2d_ee_FG[icharge][itype][iprocess] = (TH2D*)h2d_ee_FG_temp[icharge][itype][iprocess]->Clone();
        h2d_ee_FG[icharge][itype][iprocess]->Reset("ICESM");
        h2d_ee_FG[icharge][itype][iprocess]->SetName(Form("h2d_ee_FG_%s_%s_%s", charge_type[icharge], type_name[itype], process_name[iprocess]));

        int nbinsX = h2d_ee_FG_temp[icharge][itype][iprocess]->GetNbinsX();
        int nbinsY = h2d_ee_FG_temp[icharge][itype][iprocess]->GetNbinsY();

        for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

          double binwidthX = h2d_ee_FG_temp[icharge][itype][iprocess]->GetXaxis()->GetBinWidth(ibinX);

          for(int ibinY=1; ibinY < nbinsY+1; ibinY++){
  
            double content = h2d_ee_FG_temp[icharge][itype][iprocess]->GetBinContent(ibinX, ibinY);
            double err = h2d_ee_FG_temp[icharge][itype][iprocess]->GetBinError(ibinX,ibinY);
            h2d_ee_FG[icharge][itype][iprocess]->SetBinContent(ibinX, ibinY, content/binwidthX);
            h2d_ee_FG[icharge][itype][iprocess]->SetBinError(ibinX, ibinY, err/binwidthX);

          }
        }

        h2d_ee_FG[icharge][itype][iprocess]->Scale(1./normalization_factor);

      }

    }
  }

  for(int itype = 0; itype < ntype; itype++){

    for(int icharge = 0; icharge < ncharge; icharge++){

      for(int iprocess = 0; iprocess < ntype; iprocess++){

        h1d_ee_FG[icharge][itype][iprocess] = (TH1D*)h2d_ee_FG[icharge][itype][iprocess]->ProjectionX();
        h1d_ee_FG[icharge][itype][iprocess]->SetName(Form("h1d_ee_FG_%s_%s_%s", charge_type[icharge], type_name[itype], process_name[iprocess]));

      }

    }

  }

  for(int itype = 0; itype < ntype; itype++){

    for(int icharge = 0; icharge < ncharge; icharge++){

      for(int iprocess = 0; iprocess < ntype; iprocess++){

        h1d_delta_phi[icharge][itype][iprocess] = (TH1D*)h1d_delta_phi_temp[icharge][itype][iprocess]->Clone();
        h1d_delta_phi[icharge][itype][iprocess]->Reset("ICESM");
        h1d_delta_phi[icharge][itype][iprocess]->SetName(Form("h1d_delta_phi_%s_%s_%s", charge_type[icharge], type_name[itype], process_name[iprocess]));

        int nbins = h1d_delta_phi_temp[icharge][itype][iprocess]->GetNbinsX();

        for(int ibin=1; ibin < nbins+1; ibin++){

          double binwidth = h1d_delta_phi_temp[icharge][itype][iprocess]->GetXaxis()->GetBinWidth(ibin);
          double content = h1d_delta_phi_temp[icharge][itype][iprocess]->GetBinContent(ibin);
          double err = h1d_delta_phi_temp[icharge][itype][iprocess]->GetBinError(ibin);
          h1d_delta_phi[icharge][itype][iprocess]->SetBinContent(ibin, content/binwidth);
          h1d_delta_phi[icharge][itype][iprocess]->SetBinError(ibin, err/binwidth);

        
        }

        h1d_delta_phi[icharge][itype][iprocess]->Scale(1./normalization_factor);

      }

    }

  }

  for(int itype = 0; itype < ntype; itype++){

    for(int icharge = 0; icharge < ncharge; icharge++){

      for(int iprocess = 0; iprocess < ntype; iprocess++){

        hist[icharge][itype][iprocess] = (TH2D*)hist_temp[icharge][itype][iprocess]->Clone();
        hist[icharge][itype][iprocess]->Reset("ICESM");
        hist[icharge][itype][iprocess]->SetName(Form("hist_%s_%s_%s", charge_type[icharge], type_name[itype], process_name[iprocess]));

        int nbinsX = hist_temp[icharge][itype][iprocess]->GetNbinsX();
        int nbinsY = hist_temp[icharge][itype][iprocess]->GetNbinsY();

        for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

          double binwidthX = hist_temp[icharge][itype][iprocess]->GetXaxis()->GetBinWidth(ibinX);

          for(int ibinY=1; ibinY < nbinsY+1; ibinY++){
    
            double content = hist_temp[icharge][itype][iprocess]->GetBinContent(ibinX, ibinY);
            double err = hist_temp[icharge][itype][iprocess]->GetBinError(ibinX,ibinY);
            hist[icharge][itype][iprocess]->SetBinContent(ibinX, ibinY, content/binwidthX);
            hist[icharge][itype][iprocess]->SetBinError(ibinX, ibinY, err/binwidthX);

          }
        }

        hist[icharge][itype][iprocess]->Scale(1./normalization_factor);
      }

    }

  }

  for(int itype = 0; itype < ntype; itype++){

    for(int iprocess = 0; iprocess < ntype; iprocess++){

    dca_dist_single_electron[itype][iprocess] = (TH2D*)dca_dist_single_electron_temp[itype][iprocess]->Clone();
    dca_dist_single_electron[itype][iprocess]->Reset("ICESM");
    dca_dist_single_electron[itype][iprocess]->SetName(Form("dca_dist_single_electron_%s_%s", type_name[itype], process_name[iprocess]));

    int nbinsX = dca_dist_single_electron_temp[itype][iprocess]->GetNbinsX();
    int nbinsY = dca_dist_single_electron_temp[itype][iprocess]->GetNbinsY();

    for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

      double binwidthX = dca_dist_single_electron_temp[itype][iprocess]->GetXaxis()->GetBinWidth(ibinX);

      for(int ibinY=1; ibinY < nbinsY+1; ibinY++){

        double content = dca_dist_single_electron_temp[itype][iprocess]->GetBinContent(ibinX, ibinY);
        double err = dca_dist_single_electron_temp[itype][iprocess]->GetBinError(ibinX,ibinY);
        dca_dist_single_electron[itype][iprocess]->SetBinContent(ibinX, ibinY, content/binwidthX);
        dca_dist_single_electron[itype][iprocess]->SetBinError(ibinX, ibinY, err/binwidthX);

      }
    }

    dca_dist_single_electron[itype][iprocess]->Scale(1./normalization_factor);
  }

  }

  for(int itype = 0; itype < ntype; itype++){

    for(int iprocess = 0; iprocess < ntype; iprocess++){
    
    pt_spectra[itype][iprocess] = (TH1D*)pt_spectra_temp[itype][iprocess]->Clone();
    pt_spectra[itype][iprocess]->Reset("ICESM");
    pt_spectra[itype][iprocess]->SetName(Form("pt_spectra_%s_%s", type_name[itype], process_name[iprocess]));

    int nbins = pt_spectra_temp[itype][iprocess]->GetNbinsX();

    for(int ibin=1; ibin < nbins+1; ibin++){

      double binwidth = pt_spectra_temp[itype][iprocess]->GetBinWidth(ibin);
      double content = pt_spectra_temp[itype][iprocess]->GetBinContent(ibin);
      double err = pt_spectra_temp[itype][iprocess]->GetBinError(ibin);
      pt_spectra[itype][iprocess]->SetBinContent(ibin, content/binwidth);
      pt_spectra[itype][iprocess]->SetBinError(ibin, err/binwidth);

    }

    pt_spectra[itype][iprocess]->Scale(1./normalization_factor);

  }

  }

  for(int itype = 0; itype < ntype; itype++){

    for(int iprocess = 0; iprocess < ntype; iprocess++){
    
    pt_spectra_single_electron[itype][iprocess] = (TH1D*)pt_spectra_single_electron_temp[itype][iprocess]->Clone();
    pt_spectra_single_electron[itype][iprocess]->Reset("ICESM");
    pt_spectra_single_electron[itype][iprocess]->SetName(Form("pt_spectra_single_electron_%s_%s", type_name[itype], process_name[iprocess]));

    int nbins = pt_spectra_single_electron_temp[itype][iprocess]->GetNbinsX();

    for(int ibin=1; ibin < nbins+1; ibin++){

      double binwidth = pt_spectra_single_electron_temp[itype][iprocess]->GetBinWidth(ibin);
      double content = pt_spectra_single_electron_temp[itype][iprocess]->GetBinContent(ibin);
      double err = pt_spectra_single_electron_temp[itype][iprocess]->GetBinError(ibin);
      pt_spectra_single_electron[itype][iprocess]->SetBinContent(ibin, content/binwidth);
      pt_spectra_single_electron[itype][iprocess]->SetBinError(ibin, err/binwidth);

    }

    pt_spectra_single_electron[itype][iprocess]->Scale(1./normalization_factor);
  }

  }

  for(int itype = 0; itype < ntype; itype++){

    for(int iprocess = 0; iprocess < ntype; iprocess++){
    
      pt_spectra_single_positron[itype][iprocess] = (TH1D*)pt_spectra_single_positron_temp[itype][iprocess]->Clone();
      pt_spectra_single_positron[itype][iprocess]->Reset("ICESM");
      pt_spectra_single_positron[itype][iprocess]->SetName(Form("pt_spectra_single_positron_%s_%s", type_name[itype], process_name[iprocess]));

      int nbins = pt_spectra_single_positron_temp[itype][iprocess]->GetNbinsX();

      for(int ibin=1; ibin < nbins+1; ibin++){

        double binwidth = pt_spectra_single_positron_temp[itype][iprocess]->GetBinWidth(ibin);
        double content = pt_spectra_single_positron_temp[itype][iprocess]->GetBinContent(ibin);
        double err = pt_spectra_single_positron_temp[itype][iprocess]->GetBinError(ibin);
        pt_spectra_single_positron[itype][iprocess]->SetBinContent(ibin, content/binwidth);
        pt_spectra_single_positron[itype][iprocess]->SetBinError(ibin, err/binwidth);

      }

      pt_spectra_single_positron[itype][iprocess]->Scale(1./normalization_factor);
    }

  }

  for(int itype = 0; itype < ntype; itype++){

    for(int icharge = 0; icharge < ncharge; icharge++){

      for(int iprocess = 0; iprocess < ntype; iprocess++){

        h2d_delta_phi_pair_mass[icharge][itype][iprocess] = (TH2D*)h2d_delta_phi_pair_mass_temp[icharge][itype][iprocess]->Clone();
        h2d_delta_phi_pair_mass[icharge][itype][iprocess]->Reset("ICESM");
        h2d_delta_phi_pair_mass[icharge][itype][iprocess]->SetName(Form("h2d_delta_phi_pair_mass_%s_%s_%s", charge_type[icharge], type_name[itype], process_name[iprocess]));

        int nbinsX = h2d_delta_phi_pair_mass_temp[icharge][itype][iprocess]->GetNbinsX();
        int nbinsY = h2d_delta_phi_pair_mass_temp[icharge][itype][iprocess]->GetNbinsY();

        for(int ibinY=1; ibinY < nbinsY+1; ibinY++){

          double binwidthY = h2d_delta_phi_pair_mass_temp[icharge][itype][iprocess]->GetYaxis()->GetBinWidth(ibinY);

          for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

            double content = h2d_delta_phi_pair_mass_temp[icharge][itype][iprocess]->GetBinContent(ibinX, ibinY);
            double err = h2d_delta_phi_pair_mass_temp[icharge][itype][iprocess]->GetBinError(ibinX,ibinY);
            h2d_delta_phi_pair_mass[icharge][itype][iprocess]->SetBinContent(ibinX, ibinY, content/binwidthY);
            h2d_delta_phi_pair_mass[icharge][itype][iprocess]->SetBinError(ibinX, ibinY, err/binwidthY);

          }
        }

        h2d_delta_phi_pair_mass[icharge][itype][iprocess]->Scale(1./normalization_factor);
      }

    }

  }

  output->cd();

  h_debug->Write();

  for(int itype = 0; itype < ntype; itype++){

    for(int iprocess = 0; iprocess < ntype; iprocess++){

      pt_spectra_single_electron[itype][iprocess]->Write();
      pt_spectra_single_positron[itype][iprocess]->Write();
      pt_spectra[itype][iprocess]->Write();
      dca_dist_single_electron[itype][iprocess]->Write();
      
      for(int j=0; j<ncharge; j++){

        hist[j][itype][iprocess]->Write();
        hist_prof[j][itype][iprocess]->Write();
        h2d_ee_FG[j][itype][iprocess]->Write();
        h1d_ee_FG[j][itype][iprocess]->Write();
        h1d_delta_phi[j][itype][iprocess]->Write();
        h2d_delta_phi_pair_mass[j][itype][iprocess]->Write();

      }

    }


  }

  output->Close();

}
