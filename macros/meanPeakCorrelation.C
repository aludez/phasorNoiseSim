{
  //put in the two files you want to correlate 
  TString f1 = "root_runs/vpol_sim_wall_20k_phasor_space0_CH_1.root";
  TString f2 = "root_runs/vpol_sim_wall_20k_phasor_space0_CH_4.root";

  //this is how many bins you want to allow the correlation coefficient to float
  int bins = 30;

  TChain chain("traces");
  chain.Add(f1.Data());
  TChain chain2("traces");
  chain2.Add(f2.Data());
  
  TGraph * g1 = 0;
  chain.SetBranchAddress("trace", &g1);
  TGraph * g2 = 0;
  chain2.SetBranchAddress("trace", &g2);
  
  TH1F * h = new TH1F("", "", 200, 0, 1);
  
  TH1F * h3 = new TH1F("", "", 200, 0, 1);
  
  chain.GetEntry(0);

  //normalization is to correct for a 2^n thing that FFTtools does, this is for 10k point traces
  double normalization = g1->GetN()/8192.;
  if (g1->GetN() < 8192) normalization = g1->GetN()/1024.; // we also did some 2k point traces
  double best_corr = 0;
  double best_corr2 = 0;
  
  //loop and fill histograms with biggest correlation
  //skip the first entry because it was weird with fastframe mode acquisition

  for (int ii = 1; ii < chain.GetEntries(); ii++)
  {
    best_corr = 0;
    best_corr2 = 0;
    chain.GetEntry(ii);
    chain2.GetEntry(ii);
    
    double mean1 = g1->GetMean(2);
    double mean2 = g2->GetMean(2);
   
    //make sure it's all zero meaned
    for (int j=0; j < g1->GetN(); j++)
    {
      g1->GetY()[j] = g1->GetY()[j] - mean1;
      g2->GetY()[j] = g2->GetY()[j] - mean2;
    }
    
    int nPoints = g1->GetN();
    if(nPoints == 0) continue;
    
    TGraph * gc = FFTtools::getCorrelationGraph(g1,g2);
    best_corr = fabs(gc->Eval(0));
    
    int mid = gc->GetN()/2;
    //find largest correlation in specified range
    for (int iii = 1; iii< bins; iii++)
    {
      double s1 = gc->GetX()[mid + iii];
      double s2 = gc->GetX()[mid - iii];
      if(fabs(gc->Eval(s1)) > fabs(best_corr)) 
      {
        best_corr = fabs(gc->Eval(s1));
      }
      if(fabs(gc->Eval(s2)) > fabs(best_corr)) 
      {
        best_corr = fabs(gc->Eval(s2));
      }
    }
   //fill histogram with correctly normalized pearson correlation coefficient 
    h->Fill(best_corr/(g1->GetRMS(2) * g2->GetRMS(2) * normalization));
    delete gc;
  }

  //find mean peak correlation coefficient from histogram
  double mean_peak_corr = h->GetMean();

  h->Draw();
}
