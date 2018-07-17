#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "FFTtools.h"
#include "FFTWComplex.h"
#include "TH1.h"

double findMaxInRangeAboutCenter(TGraph * g, int range)
{
  int mid = g->GetN()/2;
  TGraph * gcrop = FFTtools::cropWave(g, g->GetX()[mid-range], g->GetX()[mid+range]);
  double max = TMath::MaxElement(gcrop->GetN(), gcrop->GetY());
  double min = TMath::MinElement(gcrop->GetN(), gcrop->GetY());
  if(fabs(min)>max) max = fabs(min);
  return max;
}

void meanPeakCorr(const char * fn, double distance)
{
  TString fname = fn;
  TString fname1 = fname + "_Ch1.root";
  TString fname2 = fname + "_Ch2.root";

  double dt = 2e-6/9991;
  double causalTime = distance/3e8;
  int bins = 30;//1 + TMath::Ceil(causalTime/dt);
  printf("bins = %d\n", bins);

  TChain chain("traces");
  chain.Add(fname1.Data());
  TChain chain2("traces");
  chain2.Add(fname2.Data());

  TGraph * g1 = 0;
  TGraph * g2 = 0;
  TGraph * gCorr = 0;

  chain.SetBranchAddress("trace", &g1);
  chain2.SetBranchAddress("trace", &g2);

  double mPC = 0;

  TH1D * h1 = new TH1D("h1", "h1", 200,0,1);
  chain.GetEntry(0); //skip first because it was always bad

  for (int i = 1; i < chain.GetEntries(); i++)
  {
    bins = 30;
    chain.GetEntry(i);
    chain2.GetEntry(i);
    gCorr = FFTtools::getCorrelationGraph(g1,g2);
    double normalization = 1/(g1->GetN()/8192. * g1->GetRMS(2) * g2->GetRMS(2));
    mPC = normalization * findMaxInRangeAboutCenter(gCorr, bins);
    h1->Fill(mPC);
  }

  printf("distance = %g, mpc = %g, rms = %g\n", distance, h1->GetMean(), h1->GetRMS());
  delete h1;
}

void doAll()
{
  /*
     meanPeakCorr("/project/avieregg/gno_analysis/anechoic_test_analysis/root_runs_march/telewave_hba_space0", .42);
     meanPeakCorr("/project/avieregg/gno_analysis/anechoic_test_analysis/root_runs_march/telewave_hba_space1", .42);
     meanPeakCorr("/project/avieregg/gno_analysis/anechoic_test_analysis/root_runs_march/telewave_hba_space2", .42);
     meanPeakCorr("/project/avieregg/gno_analysis/anechoic_test_analysis/root_runs_march/telewave_hba_space3", .42);
     meanPeakCorr("/project/avieregg/gno_analysis/anechoic_test_analysis/root_runs_march/telewave_hba_space4", .42);
     meanPeakCorr("/project/avieregg/gno_analysis/anechoic_test_analysis/root_runs_march/telewave_hba_space5", .42);
     */

  meanPeakCorr("/project/avieregg/gno_analysis/anechoic_test_analysis/root_runs_march/ara_hpol_parallel_0", .42);
  meanPeakCorr("/project/avieregg/gno_analysis/anechoic_test_analysis/root_runs_march/ara_hpol_parallel_1", .42);

}

