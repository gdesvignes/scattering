// Test
//
//
//
//Author: G. Desvignes - May 2016

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPaveStats.h"

#include "getopt.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"

#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"

#include "TMath.h"

#include <TApplication.h>
//#include <TGClient.h>
#include <TCanvas.h>
//#include <TRandom.h>
//#include <TGButton.h>
//#include <TRootEmbeddedCanvas.h>

using namespace std;
using namespace Pulsar;

// sigma, tau, DM, t0, Period,    A,b
int ipar1[9] = {0,1,2,3,4,5, 6,7,8};
int ipar2[9] = {0,1,2,3,4,5, 9,10,11};
int ipar3[9] = {0,1,2,3,4,5, 12,13,14};
int ipar4[9] = {0,1,2,3,4,5, 15,16,17};
int ipar5[9] = {0,1,2,3,4,5, 18,19,20};
int ipar6[9] = {0,1,2,3,4,5, 21,22,23};
int ipar7[9] = {0,1,2,3,4,5, 24,25,26};
int ipar8[9] = {0,1,2,3,4,5, 27,28,29};

struct GlobalChi2 { 
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,  ROOT::Math::IMultiGenFunction & f2, ROOT::Math::IMultiGenFunction & f3, ROOT::Math::IMultiGenFunction & f4, ROOT::Math::IMultiGenFunction & f5, ROOT::Math::IMultiGenFunction & f6, ROOT::Math::IMultiGenFunction & f7, ROOT::Math::IMultiGenFunction & f8) : 
      fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3), fChi2_4(&f4),fChi2_5(&f5), fChi2_6(&f6), fChi2_7(&f7), fChi2_8(&f8) {}

   // parameter vector is first background (in common 1 and 2) and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[9], p2[9], p3[9], p4[9]; 
      double p5[9], p6[9], p7[9], p8[9]; 

	//cout << par[0] << endl; fflush(stdout);

      for (int i = 0; i < 9; i++) {
	p1[i] = par[ipar1[i]];
        p2[i] = par[ipar2[i]];
        p3[i] = par[ipar3[i]];
        p4[i] = par[ipar4[i]];
	p5[i] = par[ipar5[i]];
        p6[i] = par[ipar6[i]];
        p7[i] = par[ipar7[i]];
        p8[i] = par[ipar8[i]];
      }

      return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3) + (*fChi2_4)(p4) + (*fChi2_5)(p5) + (*fChi2_6)(p6) + (*fChi2_7)(p7) + (*fChi2_8)(p8);
   } 

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;
   const  ROOT::Math::IMultiGenFunction * fChi2_4;
   const  ROOT::Math::IMultiGenFunction * fChi2_5;
   const  ROOT::Math::IMultiGenFunction * fChi2_6;
   const  ROOT::Math::IMultiGenFunction * fChi2_7;
   const  ROOT::Math::IMultiGenFunction * fChi2_8;
};

Double_t scatter(Double_t *x, Double_t *par) {
    Double_t sigma = par[0];
    Double_t tau = par[1];
    Double_t DM = par[2];
    Double_t t0_inf = par[3];
    Double_t Period = par[4];
    Double_t gamma = par[5];

    Double_t A = par[6];
    Double_t b = par[7];
    Double_t freq = par[8];
    Double_t delay = DM * 4.14879e3 * (1./(freq*freq)); // delay in s
    tau = tau * pow((freq*1e-3), gamma);
    //sigma = sigma / 360. * Period;

    Double_t T0 = t0_inf + delay/Period * 360.;
    while (T0 < 100) T0+=360.;
    while (T0 > 620.) T0-=360;
    //cout << x[0] << " " << T0 << endl;
    Double_t dt = (x[0] - T0) ;
    while (dt < 100) dt+=360.;
    while (dt > 180.) dt-=360.;

    //cout << x[0] << " "<< freq << " " << tau << " " << delay << " " <<  dt  << endl;

    Double_t t1 = tau * dt;
    Double_t t3 = sigma * sigma;
    Double_t t5 = tau * tau;
    Double_t t9 = TMath::Exp(-(0.4e1 * t1 - t3) / t5 / 0.4e1);
    Double_t t10 = 1.77245385090552;
    Double_t t19 = TMath::Erf((0.2e1 * t1 - t3) / sigma / tau / 0.2e1);
    Double_t result = (A*t9 * t10 * sigma * (t19 + 0.1e1) / 0.2e1+b);
    //cout << t1 << " " << t3 << " " << t5 << " " << t9 << endl;
    //cout <<"freq: "<< freq << " x: " << x[0] << " t1:  "<< t1 << " t9: " << t9 << " tau: " << tau << " delay: " << delay << " dt: " <<  dt  << endl;
    //cout <<"freq: "<< freq << " tau: " << tau << " dt: " <<  dt  << " x: " << x[0] << " result: " << result<< endl;
    return result;
}


int combinedFit(int argc, char **argv) {

  char filename[128]="test";
  int bscrunch = 1;
  double init_phase = 0.0;
  double win = 50;
  double init_dm = 0.0;
  double  win_lo1, win_hi1;
  double  win_lo2, win_hi2;
  double  win_lo3, win_hi3;
  double  win_lo4, win_hi4;
  double  win_lo5, win_hi5;
  double  win_lo6, win_hi6;
  double  win_lo7, win_hi7;
  double  win_lo8, win_hi8;
  double delay;
  double freq0, freq1, freq2, freq3;
  double freq4, freq5, freq6, freq7;
  bool fix_gamma = true;
  static struct option long_opts[] = {
      {"file",   0, NULL, 'f'},
      {0,0,0,0}
  };

  int option, opti;
  while ((option=getopt_long(argc,argv,"s:d:f:ghp:w:",long_opts,&opti))!=-1) {
      switch (option) {
          case 's':
              bscrunch = atof(optarg);
              break;
          case 'd':
              init_dm = atof(optarg);
              break;
          case 'f':
              strncpy(filename, optarg, 127);
              break;
          case 'g':
	      fix_gamma	= false;
              break;
          case 'p':
              init_phase = atof(optarg);
              break;
          case 'w':
              win = atof(optarg);
              break;
          case 'h':
          default:
              //usage();
              exit(0);
              break;
      }
  }

  Reference::To< Archive > archive = Archive::load(filename);
  if( !archive ) return -1;

  // Scrunch and Remove baseline
  archive->tscrunch();
  archive->fscrunch_to_nchan(8);
  if (bscrunch > 1) archive->bscrunch(bscrunch);
  int nchan = archive->get_nchan();
  int nbin = archive->get_nbin();
  double DM = archive->get_dispersion_measure();
  if (init_dm) DM = init_dm;


  if( archive->get_state() != Signal::Stokes) archive->convert_state(Signal::Stokes);
  archive->remove_baseline();

  // Get Data
  Pulsar::Integration* integration = archive->get_Integration(0);
  double period = integration->get_folding_period();

  // Get RMS
  vector< vector< double > > std_variance;
  integration->baseline_stats (0, &std_variance);

  double maxphase = 720.;
 
  // Create and fill histos
  TH1D * h1 = new TH1D("h1","histo 1", 2*nbin, 0, maxphase);
  TH1D * h2 = new TH1D("h2","histo 2", 2*nbin, 0, maxphase);
  TH1D * h3 = new TH1D("h3","histo 3", 2*nbin, 0, maxphase);
  TH1D * h4 = new TH1D("h4","histo 4", 2*nbin, 0, maxphase);
  TH1D * h5 = new TH1D("h5","histo 5", 2*nbin, 0, maxphase);
  TH1D * h6 = new TH1D("h6","histo 6", 2*nbin, 0, maxphase);
  TH1D * h7 = new TH1D("h7","histo 7", 2*nbin, 0, maxphase);
  TH1D * h8 = new TH1D("h8","histo 8", 2*nbin, 0, maxphase);

  char title[64];

  for (int ibin=0; ibin<2*nbin; ibin++) {
      h1->SetBinContent(ibin, integration->get_Profile(0,0)->get_amps()[ibin%nbin] / sqrt(std_variance[0][0]));
      h2->SetBinContent(ibin, integration->get_Profile(0,1)->get_amps()[ibin%nbin] / sqrt(std_variance[0][1]));
      h3->SetBinContent(ibin, integration->get_Profile(0,2)->get_amps()[ibin%nbin] / sqrt(std_variance[0][2]));
      h4->SetBinContent(ibin, integration->get_Profile(0,3)->get_amps()[ibin%nbin] / sqrt(std_variance[0][3]));
      h5->SetBinContent(ibin, integration->get_Profile(0,4)->get_amps()[ibin%nbin] / sqrt(std_variance[0][4]));
      h6->SetBinContent(ibin, integration->get_Profile(0,5)->get_amps()[ibin%nbin] / sqrt(std_variance[0][5]));
      h7->SetBinContent(ibin, integration->get_Profile(0,6)->get_amps()[ibin%nbin] / sqrt(std_variance[0][6]));
      h8->SetBinContent(ibin, integration->get_Profile(0,7)->get_amps()[ibin%nbin] / sqrt(std_variance[0][7]));
  }

  freq0 = integration->get_Profile(0,0)->get_centre_frequency();
  freq1 = integration->get_Profile(0,1)->get_centre_frequency();
  freq2 = integration->get_Profile(0,2)->get_centre_frequency();
  freq3 = integration->get_Profile(0,3)->get_centre_frequency();
  freq4 = integration->get_Profile(0,4)->get_centre_frequency();
  freq5 = integration->get_Profile(0,5)->get_centre_frequency();
  freq6 = integration->get_Profile(0,6)->get_centre_frequency();
  freq7 = integration->get_Profile(0,7)->get_centre_frequency();
  sprintf(title, "Freq %.2lf MHz", freq0);
  h1->SetTitle(title);
  sprintf(title, "Freq %.2lf MHz", freq1);
  h2->SetTitle(title);
  sprintf(title, "Freq %.2lf MHz", freq2);
  h3->SetTitle(title);
  sprintf(title, "Freq %.2lf MHz", freq2);
  h4->SetTitle(title);
  sprintf(title, "Freq %.2lf MHz", freq4);
  h5->SetTitle(title);
  sprintf(title, "Freq %.2lf MHz", freq5);
  h6->SetTitle(title);
  sprintf(title, "Freq %.2lf MHz", freq6);
  h7->SetTitle(title);
  sprintf(title, "Freq %.2lf MHz", freq7);
  h8->SetTitle(title);

  h8->SetXTitle("Rotational Phase (deg)");

  const int Npar = 6 + 3 * 8;

  double flo = integration->get_Profile(0,0)->get_centre_frequency();
  double fhi = integration->get_Profile(0,nchan-1)->get_centre_frequency();
  cout << "Frequencies: " << flo << " " << fhi << endl;

  double t0_inf = init_phase -  DM * 4.14879e3 * (1./(fhi*fhi)) / period * 360.;
  cout << "Input phase: " << init_phase << " Infinite-freq phase: " <<  t0_inf << endl;
  while (t0_inf < 0) t0_inf += 360.;
  //cout << "Input phase: " << init_phase << " Infinite-freq phase: " <<  t0_inf << endl;

  delay = DM * 4.14879e3 * (1./(freq0*freq0)) / period * 360.;
  win_lo1 = t0_inf + delay - win / 3.; win_hi1 = t0_inf + delay + win / 2.;
  while (win_lo1 > 640) win_lo1 -= 360;  while (win_hi1 > 640) win_hi1 -= 360;
  cout << win_lo1 << " " << win_hi1 << endl;
  TF1 * f1 = new TF1("f1",scatter, win_lo1, win_hi1, Npar);

  delay = DM * 4.14879e3 * (1./(freq1*freq1)) / period * 360.;
  win_lo2 = t0_inf + delay - win / 3.; win_hi2 = t0_inf + delay + win / 2.;
  while (win_lo2 > 640) win_lo2 -= 360;  while (win_hi2 > 640) win_hi2 -= 360;
  cout << win_lo2 << " " << win_hi2 << endl;
  TF1 * f2 = new TF1("f2",scatter, win_lo2, win_hi2, Npar);

  delay = DM * 4.14879e3 * (1./(freq2*freq2)) / period * 360.;
  win_lo3 = t0_inf + delay - win / 3.; win_hi3 = t0_inf + delay + win / 2.;
  while (win_lo3 > 640) win_lo3 -= 360;  while (win_hi3 > 640) win_hi3 -= 360;
  cout << win_lo3 << " " << win_hi3 << endl;
  TF1 * f3 = new TF1("f3",scatter, win_lo3, win_hi3, Npar);

  delay = DM * 4.14879e3 * (1./(freq3*freq3)) / period * 360.;
  win_lo4 = t0_inf + delay - win / 3.; win_hi4 = t0_inf + delay + win / 2.;
  while (win_lo4 > 640) win_lo4 -= 360;  while (win_hi4 > 640) win_hi4 -= 360;
  cout << win_lo4 << " " << win_hi4 << endl;
  TF1 * f4 = new TF1("f4",scatter, win_lo4, win_hi4, Npar);

  delay = DM * 4.14879e3 * (1./(freq4*freq4)) / period * 360.;
  win_lo5 = t0_inf + delay - win / 3.; win_hi5 = t0_inf + delay + win / 2.;
  while (win_lo5 > 640) win_lo5 -= 360;  while (win_hi5 > 640) win_hi5 -= 360;
  cout << win_lo5 << " " << win_hi5 << endl;
  TF1 * f5 = new TF1("f5",scatter, win_lo5, win_hi5, Npar);

  delay = DM * 4.14879e3 * (1./(freq5*freq5)) / period * 360.;
  win_lo6 = t0_inf + delay - win / 3.; win_hi6 = t0_inf + delay + win / 2.;
  while (win_lo6 > 640) win_lo6 -= 360;  while (win_hi6 > 640) win_hi6 -= 360;
  cout << win_lo6 << " " << win_hi6 << endl;
  TF1 * f6 = new TF1("f6",scatter, win_lo6, win_hi6, Npar);

  delay = DM * 4.14879e3 * (1./(freq6*freq6)) / period * 360.;
  win_lo7 = t0_inf + delay - win / 3.; win_hi7 = t0_inf + delay + win / 2.;
  while (win_lo7 > 640) win_lo7 -= 360;  while (win_hi7 > 640) win_hi7 -= 360;
  cout << win_lo7 << " " << win_hi7 << endl;
  TF1 * f7 = new TF1("f7",scatter, win_lo7, win_hi7, Npar);

  delay = DM * 4.14879e3 * (1./(freq7*freq7)) / period * 360.;
  win_lo8 = t0_inf + delay - win / 3.; win_hi8 = t0_inf + delay + win / 2.;
  while (win_lo8 > 640) win_lo8 -= 360;  while (win_hi8 > 640) win_hi8 -= 360;
  cout << win_lo8 << " " << win_hi8 << endl;
  TF1 * f8 = new TF1("f8",scatter, win_lo8, win_hi8, Npar);

  ROOT::Math::WrappedMultiTF1 wf1(*f1,1);
  ROOT::Math::WrappedMultiTF1 wf2(*f2,1);
  ROOT::Math::WrappedMultiTF1 wf3(*f3,1);
  ROOT::Math::WrappedMultiTF1 wf4(*f4,1);
  ROOT::Math::WrappedMultiTF1 wf5(*f5,1);
  ROOT::Math::WrappedMultiTF1 wf6(*f6,1);
  ROOT::Math::WrappedMultiTF1 wf7(*f7,1);
  ROOT::Math::WrappedMultiTF1 wf8(*f8,1);

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange range;
  range.SetRange(win_lo1, win_hi1);
  ROOT::Fit::BinData dataB1(opt, range); 
  ROOT::Fit::FillData(dataB1, h1);

  range.SetRange(win_lo2, win_hi2);
  ROOT::Fit::BinData dataB2(opt, range); 
  ROOT::Fit::FillData(dataB2, h2);

  range.SetRange(win_lo3, win_hi3);
  ROOT::Fit::BinData dataB3(opt, range); 
  ROOT::Fit::FillData(dataB3, h3);

  range.SetRange(win_lo4, win_hi4);
  ROOT::Fit::BinData dataB4(opt, range); 
  ROOT::Fit::FillData(dataB4, h4);

  range.SetRange(win_lo5, win_hi5);
  ROOT::Fit::BinData dataB5(opt, range); 
  ROOT::Fit::FillData(dataB5, h5);

  range.SetRange(win_lo6, win_hi6);
  ROOT::Fit::BinData dataB6(opt, range); 
  ROOT::Fit::FillData(dataB6, h6);

  range.SetRange(win_lo7, win_hi7);
  ROOT::Fit::BinData dataB7(opt, range); 
  ROOT::Fit::FillData(dataB7, h7);

  range.SetRange(win_lo8, win_hi8);
  ROOT::Fit::BinData dataB8(opt, range); 
  ROOT::Fit::FillData(dataB8, h8);

  ROOT::Fit::Chi2Function chi2_B1(dataB1, wf1);
  ROOT::Fit::Chi2Function chi2_B2(dataB2, wf2);
  ROOT::Fit::Chi2Function chi2_B3(dataB3, wf3);
  ROOT::Fit::Chi2Function chi2_B4(dataB4, wf4);
  ROOT::Fit::Chi2Function chi2_B5(dataB5, wf5);
  ROOT::Fit::Chi2Function chi2_B6(dataB6, wf6);
  ROOT::Fit::Chi2Function chi2_B7(dataB7, wf7);
  ROOT::Fit::Chi2Function chi2_B8(dataB8, wf8);

  GlobalChi2 globalChi2(chi2_B1, chi2_B2, chi2_B3, chi2_B4, chi2_B5, chi2_B6, chi2_B7, chi2_B8);

  ROOT::Fit::Fitter fitter;

  double par0[Npar] = { 15, 100, DM, t0_inf, period, -4, 20. ,0,1000., 20., 0,1000., 20., 0,1000., 20., 0,1000., 20., 0,1000., 20., 0,1000., 20., 0,1000., 20., 0,1000.};
  fitter.Config().SetParamsSettings(Npar,par0);
  fitter.Config().ParSettings(0).SetName("Sigma"); // fix 5-th parameter (Period)
  fitter.Config().ParSettings(1).SetName("Tau"); // fix 5-th parameter (Period)
  fitter.Config().ParSettings(2).SetName("DM"); // fix 5-th parameter (Period)
  fitter.Config().ParSettings(3).SetName("T0"); // fix 5-th parameter (Period)
  fitter.Config().ParSettings(4).Fix(); // fix 5-th parameter (Period)
  fitter.Config().ParSettings(4).SetName("Period"); // fix 5-th parameter (Period)
  fitter.Config().ParSettings(5).SetName("Gamma"); // fix 5-th parameter (Period)
  //fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");

  if (fix_gamma) fitter.Config().ParSettings(5).Fix();

  for (int i=0; i<nchan; i++) {
      fitter.Config().ParSettings(8+i*3).SetValue(integration->get_Profile(0,i)->get_centre_frequency());
      fitter.Config().ParSettings(8+i*3).Fix();
      sprintf(title, "Amp_%d", i);
      fitter.Config().ParSettings(6+i*3).SetName(title);
      sprintf(title, "Baseline_%d", i);
      fitter.Config().ParSettings(7+i*3).SetName(title);
      sprintf(title, "Freq_%d", i);
      fitter.Config().ParSettings(8+i*3).SetName(title);
  }

  // fit FCN function directly
  // (specify optionally data size and flag to indicate that is a chi2 fit)
  fitter.FitFCN(Npar, globalChi2,0, dataB1.Size()+dataB2.Size()+dataB3.Size()+dataB4.Size()  +  dataB5.Size()+dataB6.Size()+dataB7.Size()+dataB8.Size(),true);
  //gMinuit->SetErrorDef(4); //note 4 and not 2!
  fitter.CalculateMinosErrors();
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  cout  << integration->get_epoch().printdays(13) << " "
	<< result.Chi2() << " "
	<< result.Ndf() << " "
        << result.Value(0)*period/360.*1000. << " "
	<< result.Error(0)*period/360.*1000. << " "
	<< result.Value(1)*period/360. << " " 
	<< result.Error(1)*period/360. << " "
	<< result.Value(2) << " "
	<< result.Error(2) << endl;

  
  TCanvas * c1 = new TCanvas("Simfit","Scattering and DM fit",
                             10,10,700,700);
  gStyle->SetOptStat();
  gStyle->SetOptFit(100);
  c1->Divide(1, 9);
  c1->cd(1);
  f1->SetFitResult(result,ipar1);
  f1->SetLineColor(kRed);
  h1->GetListOfFunctions()->Add(f1);
  h1->Draw();

  
  c1->Update();
  TPaveStats *st = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats"); 
  
  sprintf(title, "#sigma (ms) = %.1f #pm %.1f", result.Value(0)*period/360.*1000., result.Error(0)*period/360.*1000.);
  st->AddText(title);
  sprintf(title, "#tau_{1GHz} (s) = %.2f #pm %.2f", result.Value(1)*period/360., result.Error(1)*period/360.);
  st->AddText(title);
  sprintf(title, "DM (pc cm^{-3}) = %.1f #pm %.1f", result.Value(2), result.Error(2));
  st->AddText(title);
  st->DrawClone();
  
  

  //h1->SetStats(0);
  //c1->Modified();

  //TPaveText *pt = new TPaveText(500, 10, 600,15);
  //pt->AddText("A TPaveText can contain severals line of text.");
  //pt->Draw();

  c1->cd(2);
  f2->SetFitResult(result,ipar2);
  f2->SetLineColor(kRed);
  h2->GetListOfFunctions()->Add(f2);
  h2->SetStats(0);
  h2->Draw();

  c1->cd(3);
  f3->SetFitResult(result,ipar3);
  f3->SetLineColor(kRed);
  h3->GetListOfFunctions()->Add(f3);
  h3->SetStats(0);
  h3->Draw();

  c1->cd(4);
  f4->SetFitResult(result,ipar4);
  f4->SetLineColor(kRed);
  h4->GetListOfFunctions()->Add(f4);
  h4->SetStats(0);
  h4->Draw();

  c1->cd(5);
  f5->SetFitResult(result,ipar5);
  f5->SetLineColor(kRed);
  h5->GetListOfFunctions()->Add(f5);
  h5->SetStats(0);
  h5->Draw();

  c1->cd(6);
  f6->SetFitResult(result,ipar6);
  f6->SetLineColor(kRed);
  h6->GetListOfFunctions()->Add(f6);
  h6->SetStats(0);
  h6->Draw();

  c1->cd(7);
  f7->SetFitResult(result,ipar7);
  f7->SetLineColor(kRed);
  h7->GetListOfFunctions()->Add(f7);
  h7->SetStats(0);
  h7->Draw();

  c1->cd(8);
  f8->SetFitResult(result,ipar8);
  f8->SetLineColor(kRed);
  h8->GetListOfFunctions()->Add(f8);
  h8->SetStats(0);
  h8->Draw();

  c1->Modified();
  c1->Update();

  return (0);
}

int main(int argc, char **argv) {
   TApplication theApp("App",&argc,argv);
   combinedFit(theApp.Argc(), theApp.Argv());
   theApp.Run();
   return 0;
}
