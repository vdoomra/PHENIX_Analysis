#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include "/direct/phenix+u/vdoomra/install/include/DileptonAnalysis/MyEvent.h"

using namespace std;
using namespace DileptonAnalysis;
const int centbin = 1;
const double step_size = 0.1;
const double pi = TMath::ACos(-1);

static const double Me = 0.000510998918;
static const double Me2 = Me*Me;
static const double c = 299792458;

const double DCENTERCUT = 0.05;
const double PC1_DPHI_CUT = 0.02, PC1_DZ_CUT = 0.5; // PC1 ghost cut

const int nbins_mass = 40;
const double mass_range[nbins_mass+1] = { 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75, 0.80, 0.87, 0.95, 1.05, 1.15, 1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.80, 2.95, 3.05, 3.10, 3.15, 3.20, 3.30, 3.45, 3.60, 3.70, 3.80, 4.00, 4.25, 4.5, 5.0, 5.5};

const int nbins_dca = 38;
const double dca_range[nbins_dca+1] = {0,40,80,120,160,200,240,280,320,360,400,440,480,520,560,600,640,680,720,760,800,840,880,920,960,1000,1040,1080,1120,1160,1200,1240,1280,1320,1360,1400,1440,1480,1520};

const int nbins_pt = 33;
const double pt_range[nbins_pt+1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.50, 6.00, 7.00, 8.00, 9.00, 10.00, 12.00, 14.00, 16.00, 18.00};

const int nbins_mass_prof = 10;
const double mass_range_prof[nbins_mass_prof+1] = {1.35, 1.65, 1.95, 2.25, 2.55, 2.85, 3.05, 3.25, 3.65, 3.85, 4.2 };

const int nlayers = 4;
const int ncharge = 2;
const int nsides = 2;
const int total_nlayers = 8;
const int counter = 3;

const char* charge_type[ncharge] = {"like","unlike"};

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

TH2D* h2d_delta_phi_pair_mass_temp[ncharge] = {NULL};
TH2D* h2d_delta_phi_pair_mass[ncharge] = {NULL};

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
  if(fabs(mytrk->GetEmcTOF()) > 5 ) return false;

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

void DCA_data_BG_phi_redef(const char* inFileA, const char* inFileB, const char* outFile, const int additional_rejection)
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

  TFile* input_A = new TFile(inFileA,"READ"); // ERT File
  if(!(input_A))
  {
    cout << "no input file" << endl;
    exit(1);
  }

  TFile* input_B = new TFile(inFileB,"READ"); // MB File
  if(!(input_B))
  {
    cout << "no input file" << endl;
    exit(1);
  }

  TTree* T_A = (TTree*)input_A->Get("T");
  TBranch* br_A = T_A->GetBranch("MyEvent");
  MyEvent* event_A = 0;
  br_A->SetAddress(&event_A);

  TTree* T_B = (TTree*)input_B->Get("T");
  TBranch* br_B = T_B->GetBranch("MyEvent");
  MyEvent* event_B = 0;
  br_B->SetAddress(&event_B);

  cout << "Trees read!" << endl; 

  TFile* output = new TFile(outFile,"RECREATE");
 
  int nevt_ert = T_A->GetEntries();
  int nevt_mb = T_B->GetEntries();

  for (int ievent_A = 0; ievent_A < nevt_ert; ievent_A++)
  {
    if (ievent_A%50000==0) cout << "Event: " << ievent_A << " / " << nevt_ert << endl;

    event_A->ClearEvent();
    br_A->GetEntry(ievent_A);
    
    int ntrack_A = event_A->GetNtrack();

    int nvtxhits_A = event_A->GetNVTXhit();
    double vtxx_A = event_A->GetPreciseX();
    double vtxy_A = event_A->GetPreciseY();
    double vtxz_A = event_A->GetPreciseZ();
    if(fabs(vtxz_A) > 10) continue;

    int ghost_A[40], real_A[40];
    int all_A[40];
    int nghost_A = 0, nreal_A = 0, nall_A = 0;


    for (int k1 = 0; k1 < ntrack_A; ++k1)
    {
      
      all_A[nall_A] = k1;
      ++nall_A;        
      int ghost_flag_A = 0;

      for (int k2 = 0; k2 < ntrack_A; ++k2)
      {
        if ( k1==k2 ) continue;
        if ( fabs((event_A->GetEntry(k1)).GetCrkphi() - (event_A->GetEntry(k2)).GetCrkphi() ) < DCENTERCUT  ) { ghost_flag_A = 1; break; }
      }
      if ( ghost_flag_A ) { ghost_A[nghost_A] = k1; ++nghost_A; }
      else { real_A[nreal_A] = k1; ++nreal_A; }
    }

    int nprimary_A = 0;
    int primary_id_A[10] = {-999};

    if(additional_rejection){

      for(int ireal_A=0; ireal_A<nreal_A; ireal_A++){

        MyTrack mytrk_A = event_A->GetEntry(real_A[ireal_A]);
        double px_A = mytrk_A.GetPtPrime()*TMath::Cos(mytrk_A.GetPhi0Prime());
        double py_A = mytrk_A.GetPtPrime()*TMath::Sin(mytrk_A.GetPhi0Prime());
        double pt_A = sqrt( px_A*px_A + py_A*py_A );

        int charge_A = mytrk_A.GetChargePrime();

        int charge_index_A = -999;
        if(charge_A<0) charge_index_A = 0;
        else charge_index_A = 1;

        int hit_index_A[total_nlayers] = {-999};

        hit_index_A[0] = mytrk_A.GetHitIndexL1();
        hit_index_A[1] = mytrk_A.GetHitIndexL2();
        hit_index_A[2] = mytrk_A.GetHitIndexL3();
        hit_index_A[3] = mytrk_A.GetHitIndexL4();
        hit_index_A[4] = mytrk_A.GetHitIndexL5();
        hit_index_A[5] = mytrk_A.GetHitIndexL6();
        hit_index_A[6] = mytrk_A.GetHitIndexL7();
        hit_index_A[7] = mytrk_A.GetHitIndexL8();

        bool veto_L1_A = false; bool veto_L2_A = false; bool veto_L3_A = false;

        for(int ilayer_A=1; ilayer_A<total_nlayers; ilayer_A++){

          if(hit_index_A[ilayer_A] == -999) continue;

          int layer_index_A = -999;
          if(ilayer_A == 1) layer_index_A = 0;
          else if(ilayer_A == 2 || ilayer_A == 3 || ilayer_A == 4) layer_index_A = 1;
          else layer_index_A = 2;

          MyVTXHit vtxhit_A = event_A->GetVTXHitEntry(hit_index_A[ilayer_A]);
          double zhit_A = vtxhit_A.GetZHit();
          double phi_hit_A = TMath::ATan2(vtxhit_A.GetYHit(), vtxhit_A.GetXHit());
          if(phi_hit_A < -pi/2) phi_hit_A += 2*pi;

          for(int ihit_A = 0; ihit_A < nvtxhits_A; ihit_A++){

            if(ihit_A == hit_index_A[ilayer_A]) continue;
            int dphi_index_A = -999;
            MyVTXHit vtx_near_hit_A = event_A->GetVTXHitEntry(ihit_A);
            int layer_A = vtx_near_hit_A.GetLayer();
            if( (ilayer_A==1) && ilayer_A != layer_A) continue;
            if( (ilayer_A==2 || ilayer_A==3 || ilayer_A==4) && layer_A !=2) continue;
            if( (ilayer_A==5 || ilayer_A==6 || ilayer_A==7) && layer_A !=3) continue;
            double phi_near_hit_A = TMath::ATan2(vtx_near_hit_A.GetYHit(), vtx_near_hit_A.GetXHit());
            if(phi_near_hit_A < -pi/2) phi_near_hit_A += 2*pi;
            double z_near_hit_A = vtx_near_hit_A.GetZHit();

            double dphi_A = phi_near_hit_A - phi_hit_A;
            double dz_A = z_near_hit_A - zhit_A;

            if(layer_A==1 && fabs(dz_A) > 0.2) continue; 
            else if((layer_A==2 || layer_A==3) && fabs(dz_A) > 0.4 ) continue;

            if(dphi_A<0) dphi_index_A = 0;
            else dphi_index_A = 1;

            double mean_A = veto_window_mean(layer_index_A, charge_index_A, dphi_index_A, pt_A);
            double sigma_A = veto_window_sigma(layer_index_A, charge_index_A, dphi_index_A, pt_A);

            if( fabs(dphi_A) > 0.001 && fabs(dphi_A) < fabs(mean_A)+15*sigma_A && layer_A==1) veto_L1_A = true;  
            if( fabs(dphi_A) > 0.001 && fabs(dphi_A) < fabs(mean_A)+15*sigma_A && layer_A==2) veto_L2_A = true; 
            if( fabs(dphi_A) > 0.001 && fabs(dphi_A) < fabs(mean_A)+15*sigma_A && layer_A==3) veto_L3_A = true;

          }
        }

        if( (veto_L1_A || veto_L2_A || veto_L3_A) ) continue;
        primary_id_A[nprimary_A] = real_A[ireal_A];
        nprimary_A++;

      }

    }

    else{

      nprimary_A = nreal_A;
      for(int ii=0; ii<nreal_A; ii++) { primary_id_A[ii] = real_A[ii];  }

    }

    for (int iprimary_A1 = 0; iprimary_A1 < nprimary_A; iprimary_A1++)
    {
      MyTrack mytrk_A1 = event_A->GetEntry(primary_id_A[iprimary_A1]);
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
      double layer2_A1 = -999;

      for(int ilayer=0; ilayer<total_nlayers; ilayer++){

        if(hit_index_A1[ilayer] == -999) continue;

        MyVTXHit vtxhit = event_A->GetVTXHitEntry(hit_index_A1[ilayer]);
        if(x1_A1 == -999){ x1_A1 = vtxhit.GetXHit();  y1_A1 = vtxhit.GetYHit(); continue; }
        if(x2_A1 == -999){ x2_A1 = vtxhit.GetXHit();  y2_A1 = vtxhit.GetYHit(); layer2_A1 = vtxhit.GetLayer(); }

      }

      double err_A1 = -999;
      if( layer2_A1 == 1) err_A1 = 0.0050;
      else err_A1 = 0.0080;

      double phi0_A1 = redefine_phi0(vtxx_A, vtxy_A, x1_A1, y1_A1, x2_A1, y2_A1, err_A1);
      double px_A1 = mytrk_A1.GetPtPrime()*TMath::Cos(phi0_A1);
      double py_A1 = mytrk_A1.GetPtPrime()*TMath::Sin(phi0_A1);
      double pz_A1 = mytrk_A1.GetPtPrime()*(TMath::Cos(mytrk_A1.GetThe0Prime()))/(TMath::Sin(mytrk_A1.GetThe0Prime()));
      double pt_A1 = sqrt( px_A1*px_A1 + py_A1*py_A1 );

      double phi_A1 = mytrk_A1.GetPhiDC();
      double zed_A1 = mytrk_A1.GetZDC();

      int charge_A1 = mytrk_A1.GetChargePrime();
      int isERT_A1 = mytrk_A1.GetisERT();
      double ecore_A1 = sqrt(pt_A1*pt_A1 + pz_A1*pz_A1 + Me2);
      bool isert_fired_A1 = false;

      double energy_threshold_A1 = -999;

      if(mytrk_A1.GetArm() == 1){  // West Arm

        if (mytrk_A1.GetSect() == 0) { energy_threshold_A1 = 1.652; } // 1.652
        else if (mytrk_A1.GetSect() == 1) { energy_threshold_A1 = 1.511; } // 1.511
        else if (mytrk_A1.GetSect() == 2) { energy_threshold_A1 = 1.541; } // 1.541
        else { energy_threshold_A1 = 1.497; } // 1.497

      }


      else {  // East Arm                                                                                                                               

        if (mytrk_A1.GetSect() == 0) { energy_threshold_A1 = 2.239; } // 2.239
        else if (mytrk_A1.GetSect() == 1) { energy_threshold_A1 = 2.561; } // 2.561
        else if(mytrk_A1.GetSect() == 2) { energy_threshold_A1 = 1.676; } // 1.676
        else { energy_threshold_A1 = 1.465; } // 1.465

      }

      isert_fired_A1 = isERT_A1 && ecore_A1 > energy_threshold_A1;

      if(!isert_fired_A1) continue;

      double crkz_A1    = mytrk_A1.GetCrkz();
      double crkphi_A1  = mytrk_A1.GetCrkphi();

      double R_A1 = pt_A1/(0.003*0.90); // Radius in cm

      TLorentzVector trk_A1;
      trk_A1.SetPxPyPzE(px_A1, py_A1, pz_A1, sqrt(px_A1*px_A1 + py_A1*py_A1 + pz_A1*pz_A1 + Me2));

      TVector3 trk3_A1;
      trk3_A1.SetXYZ(px_A1, py_A1, pz_A1);

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

      double L_A1 = sqrt((xc_A1-vtxx_A)*(xc_A1-vtxx_A) + (yc_A1-vtxy_A)*(yc_A1-vtxy_A));
      double dca_A1 = (L_A1-R_A1)*10000; // In Micro Meters

      for(int ievent_B = ievent_A+1; ievent_B < nevt_mb; ievent_B++){

        event_B->ClearEvent();
        br_B->GetEntry(ievent_B);
    
        int ntrack_B = event_B->GetNtrack();
        int nvtxhits_B = event_B->GetNVTXhit();
	      double vtxx_B = event_B->GetPreciseX();
	      double vtxy_B = event_B->GetPreciseY();
        double vtxz_B = event_B->GetPreciseZ();
        if(fabs(vtxz_B) > 10) continue;

        if(fabs(vtxz_A - vtxz_B) > 2.5) continue;  // Condition for mixing the events
    
        int ghost_B[40], real_B[40];
        int all_B[40];
        int nghost_B = 0, nreal_B = 0, nall_B = 0;
    
    
        for (int k1 = 0; k1 < ntrack_B; ++k1)
        {
          
          all_B[nall_B] = k1;
          ++nall_B;        
          int ghost_flag_B = 0;
    
          for (int k2 = 0; k2 < ntrack_B; ++k2)
          {
            if ( k1==k2 ) continue;
            if ( fabs((event_B->GetEntry(k1)).GetCrkphi() - (event_B->GetEntry(k2)).GetCrkphi() ) < DCENTERCUT ) { ghost_flag_B = 1; break; }
          }
          if ( ghost_flag_B ) { ghost_B[nghost_B] = k1; ++nghost_B; }
          else { real_B[nreal_B] = k1; ++nreal_B; }
        }
    
        int nprimary_B = 0;
        int primary_id_B[10] = {-999};

        if(additional_rejection){
    
          for(int ireal_B=0; ireal_B<nreal_B; ireal_B++){
        
            MyTrack mytrk_B = event_B->GetEntry(real_B[ireal_B]);
            double px_B = mytrk_B.GetPtPrime()*TMath::Cos(mytrk_B.GetPhi0Prime());
            double py_B = mytrk_B.GetPtPrime()*TMath::Sin(mytrk_B.GetPhi0Prime());
            double pt_B = sqrt( px_B*px_B + py_B*py_B );
    
            int charge_B = mytrk_B.GetChargePrime();
    
            int charge_index_B = -999;
            if(charge_B<0) charge_index_B = 0;
            else charge_index_B = 1;
    
            int hit_index_B[total_nlayers] = {-999};
    
            hit_index_B[0] = mytrk_B.GetHitIndexL1();
            hit_index_B[1] = mytrk_B.GetHitIndexL2();
            hit_index_B[2] = mytrk_B.GetHitIndexL3();
            hit_index_B[3] = mytrk_B.GetHitIndexL4();
            hit_index_B[4] = mytrk_B.GetHitIndexL5();
            hit_index_B[5] = mytrk_B.GetHitIndexL6();
            hit_index_B[6] = mytrk_B.GetHitIndexL7();
            hit_index_B[7] = mytrk_B.GetHitIndexL8();
    
            bool veto_L1_B = false; bool veto_L2_B = false; bool veto_L3_B = false;
    
            for(int ilayer_B=1; ilayer_B<total_nlayers; ilayer_B++){
          
              if(hit_index_B[ilayer_B] == -999) continue;
    
              int layer_index_B = -999;
              if(ilayer_B == 1) layer_index_B = 0;
              else if(ilayer_B == 2 || ilayer_B == 3 || ilayer_B == 4) layer_index_B = 1;
              else layer_index_B = 2;
    
              MyVTXHit vtxhit_B = event_B->GetVTXHitEntry(hit_index_B[ilayer_B]);
              double zhit_B = vtxhit_B.GetZHit();
              double phi_hit_B = TMath::ATan2(vtxhit_B.GetYHit(), vtxhit_B.GetXHit());
              if(phi_hit_B < -pi/2) phi_hit_B += 2*pi;
    
              for(int ihit_B = 0; ihit_B < nvtxhits_B; ihit_B++){
            
                if(ihit_B == hit_index_B[ilayer_B]) continue;
                int dphi_index_B = -999;
                MyVTXHit vtx_near_hit_B = event_B->GetVTXHitEntry(ihit_B);
                int layer_B = vtx_near_hit_B.GetLayer();
                if( (ilayer_B==1) && ilayer_B != layer_B) continue;
                if( (ilayer_B==2 || ilayer_B==3 || ilayer_B==4) && layer_B !=2) continue;
                if( (ilayer_B==5 || ilayer_B==6 || ilayer_B==7) && layer_B !=3) continue;
                double phi_near_hit_B = TMath::ATan2(vtx_near_hit_B.GetYHit(), vtx_near_hit_B.GetXHit());
                if(phi_near_hit_B < -pi/2) phi_near_hit_B += 2*pi;
                double z_near_hit_B = vtx_near_hit_B.GetZHit();
    
                double dphi_B = phi_near_hit_B - phi_hit_B;
                double dz_B = z_near_hit_B - zhit_B;
    
                if(layer_B==1 && fabs(dz_B) > 0.2) continue; 
                else if((layer_B==2 || layer_B==3) && fabs(dz_B) > 0.4 ) continue;
    
                if(dphi_B<0) dphi_index_B = 0;
                else dphi_index_B = 1;
    
                double mean_B = veto_window_mean(layer_index_B, charge_index_B, dphi_index_B, pt_B);
                double sigma_B = veto_window_sigma(layer_index_B, charge_index_B, dphi_index_B, pt_B);
    
                if( fabs(dphi_B) > 0.001 && fabs(dphi_B) < fabs(mean_B)+15*sigma_B && layer_B==1) veto_L1_B = true;  
                if( fabs(dphi_B) > 0.001 && fabs(dphi_B) < fabs(mean_B)+15*sigma_B && layer_B==2) veto_L2_B = true; 
                if( fabs(dphi_B) > 0.001 && fabs(dphi_B) < fabs(mean_B)+15*sigma_B && layer_B==3) veto_L3_B = true;
    
              }
            }
    
            if( (veto_L1_B || veto_L2_B || veto_L3_B) ) continue;
            primary_id_B[nprimary_B] = real_B[ireal_B];
            nprimary_B++;
    
          }
        }

        else{

          nprimary_B = nreal_B;
	        for (int i = 0; i<nreal_A; i++){ primary_id_B[i] = real_B[i] ;  }

        }
    
        for (int iprimary_B1 = 0; iprimary_B1 < nprimary_B; iprimary_B1++)
        {
          MyTrack mytrk_B1 = event_B->GetEntry(primary_id_B[iprimary_B1]);
          if(!single_cut(&mytrk_B1)) continue;
	  
          int hit_index_B1[total_nlayers] = {-999};

          hit_index_B1[0] = mytrk_B1.GetHitIndexL1();
          hit_index_B1[1] = mytrk_B1.GetHitIndexL2();
          hit_index_B1[2] = mytrk_B1.GetHitIndexL3();
          hit_index_B1[3] = mytrk_B1.GetHitIndexL4();
          hit_index_B1[4] = mytrk_B1.GetHitIndexL5();
          hit_index_B1[5] = mytrk_B1.GetHitIndexL6();
          hit_index_B1[6] = mytrk_B1.GetHitIndexL7();
          hit_index_B1[7] = mytrk_B1.GetHitIndexL8();

          double x1_B1 = -999; double y1_B1 = -999;
          double x2_B1 = -999; double y2_B1 = -999;
          int layer2_B1 = -999;

          for(int ilayer=0; ilayer<total_nlayers; ilayer++){

            if(hit_index_B1[ilayer] == -999) continue;

            MyVTXHit vtxhit = event_B->GetVTXHitEntry(hit_index_B1[ilayer]);
            if(x1_B1 == -999){ x1_B1 = vtxhit.GetXHit();  y1_B1 = vtxhit.GetYHit(); continue; }
            if(x2_B1 == -999){ x2_B1 = vtxhit.GetXHit();  y2_B1 = vtxhit.GetYHit(); layer2_B1 = vtxhit.GetLayer(); }

          }

          double err_B1 = -999;
          if( layer2_B1 == 1) err_B1 = 0.0050;
          else err_B1 = 0.0080;

          double phi0_B1 = redefine_phi0(vtxx_B, vtxy_B, x1_B1, y1_B1, x2_B1, y2_B1, err_B1);
	  float px_B1 = mytrk_B1.GetPtPrime()*TMath::Cos(phi0_B1);
          float py_B1 = mytrk_B1.GetPtPrime()*TMath::Sin(phi0_B1);
          float pz_B1 = mytrk_B1.GetPtPrime()*(TMath::Cos(mytrk_B1.GetThe0Prime()))/(TMath::Sin(mytrk_B1.GetThe0Prime()));
          float pt_B1 = sqrt( px_B1*px_B1 + py_B1*py_B1 );
    
          float phi_B1 = mytrk_B1.GetPhiDC();
          float zed_B1 = mytrk_B1.GetZDC();

          int charge_B1 = mytrk_B1.GetChargePrime();
          double crkz_B1    = mytrk_B1.GetCrkz();
          double crkphi_B1  = mytrk_B1.GetCrkphi();

          if (fabs(phi_A1-phi_B1)<PC1_DPHI_CUT && fabs(zed_A1-zed_B1)<PC1_DZ_CUT) continue;

          double R_B1 = pt_B1/(0.003*0.90); // Radius in cm

          TLorentzVector trk_B1;
          trk_B1.SetPxPyPzE(px_B1, py_B1, pz_B1, sqrt(px_B1*px_B1 + py_B1*py_B1 + pz_B1*pz_B1 + Me2));

          TVector3 trk3_B1;
          trk3_B1.SetXYZ(px_B1, py_B1, pz_B1);

          double a_B1 = ((x1_B1-x2_B1) -(y1_B1*y1_B1 - y2_B1*y2_B1)/(x1_B1-x2_B1))*0.5;
          double b_B1 = (y1_B1-y2_B1)/(x1_B1-x2_B1);

          double yc1_B1 = (y1_B1-a_B1*b_B1 + sqrt(R_B1*R_B1*(1+b_B1*b_B1) - (a_B1 + b_B1*y1_B1)*(a_B1 + b_B1*y1_B1)))/(1+b_B1*b_B1);
          double yc2_B1 = (y1_B1-a_B1*b_B1 - sqrt(R_B1*R_B1*(1+b_B1*b_B1) - (a_B1 + b_B1*y1_B1)*(a_B1 + b_B1*y1_B1)))/(1+b_B1*b_B1);

          double xc1_B1 = x1_B1 - a_B1 - b_B1*yc1_B1;
          double xc2_B1 = x1_B1 - a_B1 - b_B1*yc2_B1;

          double xc_B1 = -999; double yc_B1 = -999;

          if(phi0_B1 > -pi/2 && phi0_B1 < pi/2 && charge_B1 == -1){ 
	
	          if(yc1_B1<0){ yc_B1 = yc1_B1; xc_B1 = xc1_B1; }
	          else{ yc_B1 = yc2_B1; xc_B1 = xc2_B1; }

          }

          if(phi0_B1 > -pi/2 && phi0_B1 < pi/2 && charge_B1 == 1){

            if(yc1_B1>0){ yc_B1 = yc1_B1; xc_B1 = xc1_B1; }
	          else{ yc_B1 = yc2_B1; xc_B1 = xc2_B1; }

          }

          if(phi0_B1 > pi/2 && phi0_B1 < 3*pi/2 && charge_B1 == -1){

            if(yc1_B1>0){ yc_B1 = yc1_B1; xc_B1 = xc1_B1; }
	          else{ yc_B1 = yc2_B1; xc_B1 = xc2_B1; }

          }

          if(phi0_B1 > pi/2 && phi0_B1 < 3*pi/2 && charge_B1 == 1){

            if(yc1_B1<0){ yc_B1 = yc1_B1; xc_B1 = xc1_B1; }
            else{ yc_B1 = yc2_B1; xc_B1 = xc2_B1; }

          }

          double L_B1 = sqrt((xc_B1-vtxx_B)*(xc_B1-vtxx_B) + (yc_B1-vtxy_B)*(yc_B1-vtxy_B));
          double dca_B1 = (L_B1-R_B1)*10000; // In Micro Meters

          TLorentzVector pair = trk_A1 + trk_B1;
          double pair_mass = pair.M();
          double pair_pt = pair.Pt();

          double dca1 = sqrt( fabs(dca_A1*dca_A1 - dca_B1*dca_B1) );

          double angle = trk3_A1.Angle(trk3_B1);

          if(dca1 > 1520) continue;
          if(pair_pt > 5.0) continue;

          if(charge_A1 == charge_B1){

            h2d_ee_FG_temp[0]->Fill(pair_mass, pair_pt);
            hist_temp[0]->Fill(pair_mass, dca1);
            hist_prof[0]->Fill(pair_mass, dca1);	
            h1d_delta_phi_temp[0]->Fill(angle);
            h2d_delta_phi_pair_mass_temp[0]->Fill(pair_mass, angle);
  
          }
  
  
          else{ 
  
            h2d_ee_FG_temp[1]->Fill(pair_mass, pair_pt); 
            hist_temp[1]->Fill(pair_mass, dca1);
            hist_prof[1]->Fill(pair_mass, dca1);
            h2d_delta_phi_pair_mass_temp[1]->Fill(pair_mass, angle);
  
            if( (pair_mass > 1.15 && pair_mass < 2.95) || (pair_mass > 3.30 && pair_mass < 5.50) ) h1d_delta_phi_temp[1]->Fill(angle);
  
          }

        } // Track B1 Loop

      } // Event B Loop

    }  // Track A1 Loop

  } // Event A Loop

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
  
  for(int j=0; j<ncharge; j++){

    hist[j]->Write();
    hist_prof[j]->Write();

    h2d_ee_FG[j]->Write();
    h1d_ee_FG[j]->Write();
    h2d_delta_phi_pair_mass[j]->Write();

    h1d_delta_phi[j]->Write();


  }

  output->Close();

} // End make_TTree
