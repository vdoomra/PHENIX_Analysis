#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
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

const int nbins_dca = 38;
const double dca_range[nbins_dca+1] = {0,40,80,120,160,200,240,280,320,360,400,440,480,520,560,600,640,680,720,760,800,840,880,920,960,1000,1040,1080,1120,1160,1200,1240,1280,1320,1360,1400,1440,1480,1520};

const int nbins_pt = 33;
const double pt_range[nbins_pt+1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.50, 6.00, 7.00, 8.00, 9.00, 10.00, 12.00, 14.00, 16.00, 18.00};

const int nbins_mass_prof = 10;
const double mass_range_prof[nbins_mass_prof+1] = { 1.35, 1.65, 1.95, 2.25, 2.55, 2.85, 3.05, 3.25, 3.65, 3.85, 4.2 };

const char* charge_type[ncharge] = {"like", "unlike"};

//Pair Mass vs Pair DCA

TH2D* hist_temp[ncharge]= {NULL};
TH2D* hist[ncharge] = {NULL};
TProfile* hist_prof[ncharge] = {NULL};

// Pair mass vs pair pT
TH2D* h2d_ee_FG_temp[ncharge]= {NULL};
TH2D* h2d_ee_FG[ncharge]= {NULL};
TH1D* h1d_ee_FG[ncharge]= {NULL};

TH1D* h1d_delta_phi_temp[ncharge]= {NULL};
TH1D* h1d_delta_phi[ncharge]= {NULL};

TF1* func = NULL;

TH1D* pt_spectra_single_electron_temp = NULL;
TH1D* pt_spectra_single_positron_temp = NULL;
TH1D* pt_spectra_temp = NULL;

TH1D* pt_spectra_single_electron = NULL;
TH1D* pt_spectra_single_positron = NULL;
TH1D* pt_spectra = NULL;

TH2D* dca_dist_single_electron_temp = NULL;
TH2D* dca_dist_single_electron = NULL;

TH2D* h2d_delta_phi_pair_mass_temp[ncharge] = {NULL};
TH2D* h2d_delta_phi_pair_mass[ncharge] = {NULL};


const double trigger_eff_W1 = 0.9412;
const double trigger_eff_W2 = 0.9884;
const double trigger_eff_W3 = 0.9878;
const double trigger_eff_W4 = 0.9893;

const double trigger_eff_E1 = 0.9886;
const double trigger_eff_E2 = 0.9157;
const double trigger_eff_E3 = 0.9345;
const double trigger_eff_E4 = 0.9477;

double trigger_eff_W[4] = {trigger_eff_W1, trigger_eff_W2, trigger_eff_W3, trigger_eff_W4};
double trigger_eff_E[4] = {trigger_eff_E1, trigger_eff_E2, trigger_eff_E3, trigger_eff_E4};

const double energy_threshold_W1 = 1.652;
const double energy_threshold_W2 = 1.511;
const double energy_threshold_W3 = 1.541;
const double energy_threshold_W4 = 1.497;

const double energy_threshold_E1 = 2.239;
const double energy_threshold_E2 = 2.561;
const double energy_threshold_E3 = 1.676;
const double energy_threshold_E4 = 1.465;


double offset = 0.0;

double get_trigger_efficiency(int arm, int sect) {

  if (sect < 0 || sect > 3) return 1.0;
  return (arm == 1) ? trigger_eff_W[sect] : trigger_eff_E[sect];

}

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

  if(mytrk->GetN0() < 2) return false;
  if((mytrk->GetEcore()/p) < 0.7) return false;
  if(pt < 0.4) return false;
  if( fabs(mytrk->GetEmcTOF()) > 5 ) return false;

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

void DCA_data_phi_redef(const char* inFile, const char* outFile, const int isert, const int additional_rejection)
{

  
  for(int j=0; j<ncharge; j++){

    hist_temp[j] = new TH2D(Form("hist_%s_temp", charge_type[j]),"DCA Distribution", nbins_mass, mass_range, nbins_dca, dca_range);
    hist_temp[j]->Sumw2();

    hist_prof[j] = new TProfile(Form("hist_prof_%s", charge_type[j]),"DCA Distribution Profile Histogram", nbins_mass_prof, mass_range_prof);
    hist_prof[j]->Sumw2();

    h2d_ee_FG_temp[j] = new TH2D(Form("h2d_ee_FG_%s_temp", charge_type[j]), "ee mass distribution FG", nbins_mass, mass_range ,100, 0, 10);
    h2d_ee_FG_temp[j]->Sumw2();

    h1d_delta_phi_temp[j] = new TH1D(Form("h1d_delta_phi_%s_temp",charge_type[j]), "#Delta Phi Distribution", 50, 0, pi);
    h1d_delta_phi_temp[j]->Sumw2();

    h2d_delta_phi_pair_mass_temp[j] = new TH2D(Form("h2d_delta_phi_pair_mass_%s_temp", charge_type[j]), "#Delta #phi vs Pair Mass Distribution", nbins_mass, mass_range, 10, 0, TMath::Pi());
    h2d_delta_phi_pair_mass_temp[j]->Sumw2();


  }

  pt_spectra_single_electron_temp = new TH1D("pt_spectra_single_electron_temp", "pT spectra", nbins_pt, pt_range);
  pt_spectra_single_electron_temp->Sumw2();

  pt_spectra_single_positron_temp = new TH1D("pt_spectra_single_positron_temp", "pT spectra", nbins_pt, pt_range);
  pt_spectra_single_positron_temp->Sumw2();

  pt_spectra_temp = new TH1D("pt_spectra_temp", "pT spectra", nbins_pt, pt_range);
  pt_spectra_temp->Sumw2();

  dca_dist_single_electron_temp = new TH2D("dca_dist_single_electron_temp", "DCA_distribution", nbins_pt, pt_range, 100, -1600, 1600);
  dca_dist_single_electron_temp->Sumw2();

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

  TFile* output = new TFile(outFile,"RECREATE");
 
  int nevt = T->GetEntries();

  for (int ievent = 0; ievent < nevt; ievent++)
  {
    if (ievent%100000==0) cout << "Event: " << ievent << " / " << nevt << endl;

    event->ClearEvent();
    br->GetEntry(ievent);
    
    int ntrack = event->GetNtrack();
    int nvtxhits = event->GetNVTXhit();
    int nclust = event->GetNcluster();
    float vtxz = event->GetPreciseZ();
    float vtxx = event->GetPreciseX();
    float vtxy = event->GetPreciseY();
    if( fabs(vtxz) > 10 ) continue;

    int ghost[40], real[40];
    int all[40];
    int nghost = 0, nreal = 0, nall = 0;
   
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
      float px = mytrk.GetPtPrime()*TMath::Cos(mytrk.GetPhi0Prime());
      float py = mytrk.GetPtPrime()*TMath::Sin(mytrk.GetPhi0Prime());
      float pz = mytrk.GetPtPrime()*(TMath::Cos(mytrk.GetThe0Prime()))/(TMath::Sin(mytrk.GetThe0Prime()));
      float pt = sqrt( px*px + py*py );
      
      
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

    if(nprimary<1) continue;

    }

    else{

      nprimary = nreal;
      for(int ii=0; ii<nreal; ii++) { primary_id[ii] = real[ii];  }

    }

    double weight = 1.0;

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

      double phi0_A1 = redefine_phi0(vtxx, vtxy, x1_A1, y1_A1, x2_A1, y2_A1, err_A1);
      double px_A1 = mytrk_A1.GetPtPrime()*TMath::Cos(phi0_A1);
      double py_A1 = mytrk_A1.GetPtPrime()*TMath::Sin(phi0_A1);
      double pz_A1 = mytrk_A1.GetPtPrime()*(TMath::Cos(mytrk_A1.GetThe0Prime()))/(TMath::Sin(mytrk_A1.GetThe0Prime()));
      double pt_A1 = sqrt(px_A1*px_A1 + py_A1*py_A1);

      double p_A1 = sqrt(pt_A1*pt_A1 + pz_A1*pz_A1);
      double x_A1 = 220*sin(mytrk_A1.GetThe0Prime())*cos(mytrk_A1.GetPhi0Prime());
      double y_A1 = 220*sin(mytrk_A1.GetThe0Prime())*sin(mytrk_A1.GetPhi0Prime());
      double z_A1 = 220*cos(mytrk_A1.GetThe0Prime());

      double phi_A1 = mytrk_A1.GetPhiDC();
      double zed_A1 = mytrk_A1.GetZDC();

      int charge_A1 = mytrk_A1.GetChargePrime();
      int emcid_A1 = mytrk_A1.GetEmcId();
      int arm_A1 = mytrk_A1.GetArm();
      int isERT_A1 = mytrk_A1.GetisERT();
      double ecore_A1 = sqrt(pt_A1*pt_A1 + pz_A1*pz_A1 + Me2);
      bool isert_fired_A1 = false;

      double energy_threshold_A1 = -999;

      if(mytrk_A1.GetArm() == 1){  // West Arm

	      if (mytrk_A1.GetSect() == 0) { energy_threshold_A1 = energy_threshold_W1 + offset; }  // 1.652
	      else if (mytrk_A1.GetSect() == 1) { energy_threshold_A1 = energy_threshold_W2 + offset; } //1.511
	      else if (mytrk_A1.GetSect() == 2) { energy_threshold_A1 = energy_threshold_W3 + offset; } // 1.541
	      else { energy_threshold_A1 = energy_threshold_W4 + offset; } //1.497

      }

      else { // East Arm

	      if (mytrk_A1.GetSect() == 0) { energy_threshold_A1 = energy_threshold_E1 + offset; } // 2.239
	      else if (mytrk_A1.GetSect() == 1) { energy_threshold_A1 = energy_threshold_E2 + offset; } // 2.561
	      else if(mytrk_A1.GetSect() == 2) { energy_threshold_A1 = energy_threshold_E3 + offset; } // 1.676
	      else { energy_threshold_A1 = energy_threshold_E4 + offset; } // 1.465

      }

      isert_fired_A1 = isERT_A1 && ecore_A1 > energy_threshold_A1;

      if(isert){

        if(isert_fired_A1 && charge_A1 == -1) pt_spectra_single_electron_temp->Fill(pt_A1);
        if(isert_fired_A1 && charge_A1 == 1) pt_spectra_single_positron_temp->Fill(pt_A1);
        if(isert_fired_A1) pt_spectra_temp->Fill(pt_A1);

      }

      else{

        if( charge_A1 == -1) pt_spectra_single_electron_temp->Fill(pt_A1);
        if( charge_A1 == 1) pt_spectra_single_positron_temp->Fill(pt_A1);
        pt_spectra_temp->Fill(pt_A1);

      }

      TLorentzVector trk_A1;
      trk_A1.SetPxPyPzE(px_A1, py_A1, pz_A1, ecore_A1);

      TVector3 trk3_A1;
      trk3_A1.SetXYZ(px_A1, py_A1, pz_A1);

      double R_A1 = pt_A1/(0.003*0.90); // Radius in cm

      double a_A1 = ((x1_A1-x2_A1) -(y1_A1*y1_A1 - y2_A1*y2_A1)/(x1_A1-x2_A1))*0.5;
      double b_A1 = (y1_A1-y2_A1)/(x1_A1-x2_A1);

      double yc1_A1 = (y1_A1-a_A1*b_A1 + sqrt(R_A1*R_A1*(1+b_A1*b_A1) - (a_A1 + b_A1*y1_A1)*(a_A1 + b_A1*y1_A1)))/(1+b_A1*b_A1);
      double yc2_A1 = (y1_A1-a_A1*b_A1 - sqrt(R_A1*R_A1*(1+b_A1*b_A1) - (a_A1 + b_A1*y1_A1)*(a_A1 + b_A1*y1_A1)))/(1+b_A1*b_A1);

      double xc1_A1 = x1_A1 - a_A1 - b_A1*yc1_A1;
      double xc2_A1 = x1_A1 - a_A1 - b_A1*yc2_A1;

      double xc_A1 = -999; double yc_A1 = -999;

      if( phi0_A1 > -pi/2 && phi0_A1 < pi/2 && charge_A1 == -1 ){ 
	
	      if(yc1_A1<0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
	      else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      if( phi0_A1 > -pi/2 && phi0_A1 < pi/2 && charge_A1 == 1 ){

        if(yc1_A1>0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
	      else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      if( phi0_A1 > pi/2 && phi0_A1 < 3*pi/2 && charge_A1 == -1 ){

        if(yc1_A1>0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
	      else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      if( phi0_A1 > pi/2 && phi0_A1 < 3*pi/2 && charge_A1 == 1 ){

        if(yc1_A1<0){ yc_A1 = yc1_A1; xc_A1 = xc1_A1; }
        else{ yc_A1 = yc2_A1; xc_A1 = xc2_A1; }

      }

      double L_A1 = sqrt((xc_A1-vtxx)*(xc_A1-vtxx) + (yc_A1-vtxy)*(yc_A1-vtxy));
      double dca_A1 = (L_A1-R_A1)*10000; // In Micro Meters

      if(isert) { if(isert_fired_A1) dca_dist_single_electron_temp->Fill(pt_A1, dca_A1); }
      else dca_dist_single_electron_temp->Fill(pt_A1, dca_A1);

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

        double phi0_A2 = redefine_phi0(vtxx, vtxy, x1_A2, y1_A2, x2_A2, y2_A2, err_A2);
        double px_A2 = mytrk_A2.GetPtPrime()*TMath::Cos(phi0_A2);
        double py_A2 = mytrk_A2.GetPtPrime()*TMath::Sin(phi0_A2);
        double pz_A2 = mytrk_A2.GetPtPrime()*(TMath::Cos(mytrk_A2.GetThe0Prime()))/(TMath::Sin(mytrk_A2.GetThe0Prime()));
        double pt_A2 = sqrt(px_A2*px_A2 + py_A2*py_A2);

        int charge_A2 = mytrk_A2.GetChargePrime();
        int emcid_A2 = mytrk_A2.GetEmcId();
        int arm_A2 = mytrk_A2.GetArm();
        int isERT_A2 = mytrk_A2.GetisERT();
        double ecore_A2 = sqrt(pt_A2*pt_A2 + pz_A2*pz_A2 + Me2);
        
        double phi_A2 = mytrk_A2.GetPhiDC();
        double zed_A2 = mytrk_A2.GetZDC();

        bool isert_fired_A2 = false;
        double energy_threshold_A2 = -999;

        if(mytrk_A2.GetArm() == 1){  // West Arm
	  
          if (mytrk_A2.GetSect() == 0) { energy_threshold_A2 = energy_threshold_W1 + offset; } //1.652
          else if (mytrk_A2.GetSect() == 1) { energy_threshold_A2 = energy_threshold_W2 + offset; } // 1.511
          else if (mytrk_A2.GetSect() == 2) { energy_threshold_A2 = energy_threshold_W3 + offset; } // 1.541
          else { energy_threshold_A2 = energy_threshold_W4 + offset; } // 1.497
          
        }


        else {  // East Arm                                                                                                                               

          if (mytrk_A2.GetSect() == 0) { energy_threshold_A2 = energy_threshold_E1 + offset; } // 2.239
          else if (mytrk_A2.GetSect() == 1) { energy_threshold_A2 = energy_threshold_E2 + offset; } // 2.561
          else if(mytrk_A2.GetSect() == 2) { energy_threshold_A2 = energy_threshold_E3 + offset; } // 1.676
          else { energy_threshold_A2 = energy_threshold_E4 + offset; } // 1.465

        }

        isert_fired_A2 = isERT_A2 && ecore_A2 > energy_threshold_A2;

        if (isert && !(isert_fired_A1 || isert_fired_A2)) continue;
	      if (fabs(phi_A1-phi_A2)<PC1_DPHI_CUT && fabs(zed_A1-zed_A2)<PC1_DZ_CUT) continue;

        if(isert){

          if (isert_fired_A1 && !isert_fired_A2) weight /= get_trigger_efficiency(mytrk_A1.GetArm(), mytrk_A1.GetSect());
          else if (isert_fired_A2 && !isert_fired_A1) weight /= get_trigger_efficiency(mytrk_A2.GetArm(), mytrk_A2.GetSect());
          else { 

            if(gRandom->Uniform(0,1) > 0.5) weight /= get_trigger_efficiency(mytrk_A1.GetArm(), mytrk_A1.GetSect()); 
            else weight /= get_trigger_efficiency(mytrk_A2.GetArm(), mytrk_A2.GetSect()); 
          
          }

        }

        TLorentzVector trk_A2;
        trk_A2.SetPxPyPzE(px_A2, py_A2, pz_A2, ecore_A2);

	      TVector3 trk3_A2;
	      trk3_A2.SetXYZ(px_A2, py_A2, pz_A2);
	      
        double R_A2 = pt_A2/(0.003*0.90); // Radius in cm
     
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

        double L_A2 = sqrt((xc_A2-vtxx)*(xc_A2-vtxx) + (yc_A2-vtxy)*(yc_A2-vtxy));
        double dca_A2 = (L_A2-R_A2)*10000; // In Micro Meters

        TLorentzVector pair = trk_A1 + trk_A2;
        double pair_mass = pair.M();
        double pair_pt = pair.Pt();
	      double dca1 = sqrt( fabs(dca_A1*dca_A1 - dca_A2*dca_A2) );
	      
        double angle = trk3_A1.Angle(trk3_A2);

        if(dca1 > 1520) continue;
        if(pair_pt > 5.0) continue;

	      
	      if(charge_A1 == charge_A2){

          h2d_ee_FG_temp[0]->Fill(pair_mass, pair_pt, weight);
          hist_temp[0]->Fill(pair_mass, dca1, weight);
          hist_prof[0]->Fill(pair_mass, dca1, weight);	
          h1d_delta_phi_temp[0]->Fill(angle, weight);
	        h2d_delta_phi_pair_mass_temp[0]->Fill(pair_mass, angle, weight);

	      }


	      else{ 

          h2d_ee_FG_temp[1]->Fill(pair_mass, pair_pt, weight); 
          hist_temp[1]->Fill(pair_mass, dca1, weight);
          hist_prof[1]->Fill(pair_mass, dca1, weight);
	        h2d_delta_phi_pair_mass_temp[1]->Fill(pair_mass, angle, weight);

          if( (pair_mass > 1.15 && pair_mass < 4.50) ) h1d_delta_phi_temp[1]->Fill(angle, weight);

	      }

      } // End Track A2 Loop
    
    } // End Track A1 Loop
  } // End Event A loop

  for(int icharge = 0; icharge < ncharge; icharge++){

    h2d_ee_FG[icharge] = (TH2D*)h2d_ee_FG_temp[icharge]->Clone();
    h2d_ee_FG[icharge]->Reset("ICESM");
    h2d_ee_FG[icharge]->SetName(Form("h2d_ee_FG_%s", charge_type[icharge]));

    int nbinsX = h2d_ee_FG_temp[icharge]->GetNbinsX();
    int nbinsY = h2d_ee_FG_temp[icharge]->GetNbinsY();

    for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

      double binwidthX = h2d_ee_FG_temp[icharge]->GetXaxis()->GetBinWidth(ibinX);

      for(int ibinY=1; ibinY < nbinsY+1; ibinY++){
 
        double content = h2d_ee_FG_temp[icharge]->GetBinContent(ibinX, ibinY);
        double err = h2d_ee_FG_temp[icharge]->GetBinError(ibinX,ibinY);
        h2d_ee_FG[icharge]->SetBinContent(ibinX, ibinY, content/binwidthX);
        h2d_ee_FG[icharge]->SetBinError(ibinX, ibinY, err/binwidthX);

      }
    }
  }

  for(int icharge = 0; icharge < ncharge; icharge++){

    h1d_ee_FG[icharge] = (TH1D*)h2d_ee_FG[icharge]->ProjectionX();
    h1d_ee_FG[icharge]->SetName(Form("h1d_ee_FG_%s", charge_type[icharge]));

  }


  for(int icharge = 0; icharge < ncharge; icharge++){

    h1d_delta_phi[icharge] = (TH1D*)h1d_delta_phi_temp[icharge]->Clone();
    h1d_delta_phi[icharge]->Reset("ICESM");
    h1d_delta_phi[icharge]->SetName(Form("h1d_delta_phi_%s", charge_type[icharge]));

    int nbins = h1d_delta_phi_temp[icharge]->GetNbinsX();

    for(int ibin=1; ibin < nbins+1; ibin++){

      double binwidth = h1d_delta_phi_temp[icharge]->GetXaxis()->GetBinWidth(ibin);
      double content = h1d_delta_phi_temp[icharge]->GetBinContent(ibin);
      double err = h1d_delta_phi_temp[icharge]->GetBinError(ibin);
      h1d_delta_phi[icharge]->SetBinContent(ibin, content/binwidth);
      h1d_delta_phi[icharge]->SetBinError(ibin, err/binwidth);

      
    }
  }

  for(int icharge = 0; icharge < ncharge; icharge++){

    hist[icharge] = (TH2D*)hist_temp[icharge]->Clone();
    hist[icharge]->Reset("ICESM");
    hist[icharge]->SetName(Form("hist_%s", charge_type[icharge]));

    int nbinsX = hist_temp[icharge]->GetNbinsX();
    int nbinsY = hist_temp[icharge]->GetNbinsY();

    for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

      double binwidthX = hist_temp[icharge]->GetXaxis()->GetBinWidth(ibinX);

      for(int ibinY=1; ibinY < nbinsY+1; ibinY++){
 
        double content = hist_temp[icharge]->GetBinContent(ibinX, ibinY);
        double err = hist_temp[icharge]->GetBinError(ibinX,ibinY);
        hist[icharge]->SetBinContent(ibinX, ibinY, content/binwidthX);
        hist[icharge]->SetBinError(ibinX, ibinY, err/binwidthX);

      }
    }
  }

  dca_dist_single_electron = (TH2D*)dca_dist_single_electron_temp->Clone();
  dca_dist_single_electron->Reset("ICESM");
  dca_dist_single_electron->SetName("dca_dist_single_electron");

  int nbinsX = dca_dist_single_electron_temp->GetNbinsX();
  int nbinsY = dca_dist_single_electron_temp->GetNbinsY();

  for(int ibinX=1; ibinX < nbinsX+1; ibinX++){

    double binwidthX = dca_dist_single_electron_temp->GetXaxis()->GetBinWidth(ibinX);

    for(int ibinY=1; ibinY < nbinsY+1; ibinY++){

      double content = dca_dist_single_electron_temp->GetBinContent(ibinX, ibinY);
      double err = dca_dist_single_electron_temp->GetBinError(ibinX,ibinY);
      dca_dist_single_electron->SetBinContent(ibinX, ibinY, content/binwidthX);
      dca_dist_single_electron->SetBinError(ibinX, ibinY, err/binwidthX);

    }
  }
    
  pt_spectra = (TH1D*)pt_spectra_temp->Clone();
  pt_spectra->Reset("ICESM");
  pt_spectra->SetName("pt_spectra");

  int nbins = pt_spectra_temp->GetNbinsX();

  for(int ibin=1; ibin < nbins+1; ibin++){

    double binwidth = pt_spectra_temp->GetBinWidth(ibin);
    double content = pt_spectra_temp->GetBinContent(ibin);
    double err = pt_spectra_temp->GetBinError(ibin);
    pt_spectra->SetBinContent(ibin, content/binwidth);
    pt_spectra->SetBinError(ibin, err/binwidth);

  }

  pt_spectra_single_electron = (TH1D*)pt_spectra_single_electron_temp->Clone();
  pt_spectra_single_electron->Reset("ICESM");
  pt_spectra_single_electron->SetName("pt_spectra_single_electron");

  nbins = pt_spectra_single_electron_temp->GetNbinsX();

  for(int ibin=1; ibin < nbins+1; ibin++){

    double binwidth = pt_spectra_single_electron_temp->GetBinWidth(ibin);
    double content = pt_spectra_single_electron_temp->GetBinContent(ibin);
    double err = pt_spectra_single_electron_temp->GetBinError(ibin);
    pt_spectra_single_electron->SetBinContent(ibin, content/binwidth);
    pt_spectra_single_electron->SetBinError(ibin, err/binwidth);

  }

  pt_spectra_single_positron = (TH1D*)pt_spectra_single_positron_temp->Clone();
  pt_spectra_single_positron->Reset("ICESM");
  pt_spectra_single_positron->SetName("pt_spectra_single_positron");

  nbins = pt_spectra_single_positron_temp->GetNbinsX();

  for(int ibin=1; ibin < nbins+1; ibin++){

    double binwidth = pt_spectra_single_positron_temp->GetBinWidth(ibin);
    double content = pt_spectra_single_positron_temp->GetBinContent(ibin);
    double err = pt_spectra_single_positron_temp->GetBinError(ibin);
    pt_spectra_single_positron->SetBinContent(ibin, content/binwidth);
    pt_spectra_single_positron->SetBinError(ibin, err/binwidth);

  }
  
  for(int icharge = 0; icharge < ncharge; icharge++){

    h2d_delta_phi_pair_mass[icharge] = (TH2D*)h2d_delta_phi_pair_mass_temp[icharge]->Clone();
    h2d_delta_phi_pair_mass[icharge]->Reset("ICESM");
    h2d_delta_phi_pair_mass[icharge]->SetName(Form("h2d_delta_phi_pair_mass_%s", charge_type[icharge]));
    
    int nbinsX = h2d_delta_phi_pair_mass_temp[icharge]->GetNbinsX();
    int nbinsY = h2d_delta_phi_pair_mass_temp[icharge]->GetNbinsY();
    
    for(int ibinY=1; ibinY < nbinsY+1; ibinY++){
    
      double binwidthY = h2d_delta_phi_pair_mass_temp[icharge]->GetYaxis()->GetBinWidth(ibinY);
    
      for(int ibinX=1; ibinX < nbinsX+1; ibinX++){
    
	double content = h2d_delta_phi_pair_mass_temp[icharge]->GetBinContent(ibinX, ibinY);
	double err = h2d_delta_phi_pair_mass_temp[icharge]->GetBinError(ibinX,ibinY);
	h2d_delta_phi_pair_mass[icharge]->SetBinContent(ibinX, ibinY, content/binwidthY);
	h2d_delta_phi_pair_mass[icharge]->SetBinError(ibinX, ibinY, err/binwidthY);
    
      }
    }

  }

  output->cd();
  
  pt_spectra_single_electron->Write();
  pt_spectra_single_positron->Write();
  pt_spectra->Write();
  dca_dist_single_electron->Write();
  
  for(int j=0; j<ncharge; j++){

    hist[j]->Write();
    hist_prof[j]->Write();

    h2d_ee_FG[j]->Write();
    h1d_ee_FG[j]->Write();
    h2d_delta_phi_pair_mass[j]->Write();

    h1d_delta_phi[j]->Write();


  }

  output->Close();

} 
