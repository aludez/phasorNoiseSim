{
  //put in the two files you want to correlate 
  TString f1 = "root_runs/vpol_sim_wall_20k_phasor_space0_CH_1.root";
  TString f2 = "root_runs/vpol_sim_wall_20k_phasor_space0_CH_4.root";

  //pick the event you want to correlate
  int event = 1;

  TChain chain("traces");
  chain.Add(f1.Data());
  TChain chain2("traces");
  chain2.Add(f2.Data());
  
  TGraph * g1 = 0;
  chain.SetBranchAddress("trace", &g1);
  TGraph * g2 = 0;
  chain2.SetBranchAddress("trace", &g2);
  
  chain.GetEntry(0);

  //normalization is to correct for a 2^n thing that FFTtools does, this is for 10k point traces
  double normalization = g1->GetN()/8192.;
  if (g1->GetN() < 8192) normalization = g1->GetN()/1024.; // we also did some 2k point traces
  
  //skip the first entry because it was weird with fastframe mode acquisition
  chain.GetEntry(event);
  chain2.GetEntry(event);
  
  double mean1 = g1->GetMean(2);
  double mean2 = g2->GetMean(2);
  
  //make sure it's all zero meaned
  for (int j=0; j < g1->GetN(); j++)
  {
    g1->GetY()[j] = g1->GetY()[j] - mean1;
    g2->GetY()[j] = g2->GetY()[j] - mean2;
  }
  
  TGraph * gc = FFTtools::getCorrelationGraph(g1,g2);
  
  double g1RMS = g1->GetRMS(2);
  double g2RMS = g2->GetRMS(2);

  //loop thru graph and normalize correlation values correctly
  for (int j=0; j < gc->GetN(); j++)
  {
    gc->GetY()[j] = gc->GetY()[j]/(g1RMS * g2RMS * normalization);
  }

  gc->Draw("alp");
}
