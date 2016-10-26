#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TString.h"
#define FFTW_64_BIT

#include "/project/avieregg/gno_analysis/libRootFftwWrapper/include/FFTtools.h"
#include "/project/avieregg/gno_analysis/libRootFftwWrapper/include/FFTWComplex.h"
#include "TMath.h"

/*
phasorModel inputs:

dist_name is a string that gets appended to the output root files

ant_sep is the distance betweeen antenna feeds

phasor_num is the number of phasors used to generate each trace

geometry is a string that picks the geometry of the anechoic chamber. currently supported are "box", "shell". Where shell is a spherical shell and the other is a box

traces is the number of traces to be saved in the root file

antenna_type is either "vpol", "hba" or "lba" where lba and hba are the telewaves we used and vpol is the ara vpol.  This just changes what frequency response gets used

phasorModel outputs:
Two .root files are output with names corresponding to fname and fname2 in the code

At the bottom of the file is a doAll() function, which is just a placeholder to run this a bunch of times with different settings.
*/

void phasorModel(const char * dist_name, double ant_sep, int phasor_num=20000, const char * geometry = "box", int traces=500, const char * antenna_type = "vpol")
{
  //set up where to save the file, get TTrees and TGraphs ready  
  TString dn = dist_name;
  
  TString fname = "root_runs/vpol_sim_wall_20k_phasor_" + dn + "_CH_1.root"; 
  TFile file(fname.Data(), "RECREATE");
  TTree * tree = new TTree("traces", "traces");
  
  TString fname2 = "root_runs/vpol_sim_wall_20k_phasor_" + dn + "_CH_4.root"; 
  TFile file2(fname2.Data(), "RECREATE");
  TTree * tree2 = new TTree("traces", "traces");

  TGraph * g = 0;
  TGraph * g2 = 0;

  int num_traces = traces;

  //size of the shell and sphere
  double radius_s = 5.;
  double x_length = radius_s;
  double y_length = radius_s;
  double z_length = radius_s;

  //size of the wall
  double wall_size_x = 1.5;
  double wall_size_y = 1.5;
  double wall_size_z = 4.;

  //antenna positions
  double ant_x = 0;
  double ant_y = 0;
  double ant_zA = (ant_sep/2);
  double ant_zB = -1 * (ant_sep/2);

  //for generating phasors
  int num_phasors = phasor_num;
  TRandom3 tr1;
  tr1.SetSeed(0);
  double amp_mean = 9.25;
  TGraph * ant_response;
  //the numbers for shell and box are empirically determined to be the correct phasor amplitude for 100000 phasors.  the sqrt thing is scaling for different numbers of phasors.
  //also load up the frequency response because it is a convenient time to do so
  if (strcmp(antenna_type, "vpol") == 0)
  {
    if (strcmp(geometry, "shell") == 0) amp_mean = 9.25 * sqrt(100000/phasor_num);
    if (strcmp(geometry, "box") == 0) amp_mean = 4.1 * sqrt(100000/phasor_num);
    ant_response = new TGraph("L.csv");
  }
  
  if (strcmp(antenna_type, "hba") == 0)
  {
    if (strcmp(geometry, "shell") == 0) amp_mean = 13.5 * sqrt(100000/phasor_num);
    if (strcmp(geometry, "box") == 0) amp_mean = 6. * sqrt(100000/phasor_num);
    ant_response = new TGraph("H.csv");
  }
  
  if (strcmp(antenna_type, "lba") == 0)
  {
    if (strcmp(geometry, "shell") == 0) amp_mean = 14.5 * sqrt(100000/phasor_num);
    if (strcmp(geometry, "box") == 0) amp_mean = 6.4 * sqrt(100000/phasor_num);
    ant_response = new TGraph("L.csv");
  }
  
  double phasor_x = 0;
  double phasor_y = 0;
  double phasor_z = 0;
  double phasor_phase = 0;
  double amp = 0;


  double x_time[9991];

  //for the amp noise
  double amp_noise = 9/sqrt(5); //this is ~75k noise temp amp noise
  double amp_phase_A = 0;
  double amp_phase_B = 0;
  double amp_mag_A = 0;
  double amp_mag_B = 0;
  
  //arrays that will be used later for FFT components
  FFTWComplex fftw_A[4996];
  FFTWComplex fftw_B[4996];
  
  double * yAfilt;
  double * yBfilt;
  
  //time and frequency steps for TGraphs
  double dt = 2*pow(10,-6)/9991;
  double dF = 1/( 9991 * dt);

  //distances to antennae
  double dA = 0;
  double dB = 0;
  
  double phase = 0;
  double theta = 0;
  double rs = 0;

  FFTWComplex temp_A;
  FFTWComplex temp_B;
  
  //phase difference accumulated by waveform on the way to antenna A or B
  double dAphase = 0;
  double dBphase = 0;
  
  double filter_gain = 0;

  //set up what we want to put in the TTree
  tree->Branch("trace", &g);
  tree->Branch("FFTW", fftw_A);
  tree2->Branch("trace", &g2);
  tree2->Branch("FFTW", fftw_B);
  
  
  for (int t = 0; t < 9991; t++) x_time[t] = dt * t;
  //filter response from data sheet
  TGraph * hp = new TGraph("high_pass.txt");
  for (int j = 0; j < hp->GetN(); j++) hp->GetY()[j] = (pow(10, -hp->GetY()[j]/10));
  
  TGraph * lp = new TGraph("low_pass2.txt");
  for (int j = 0; j < lp->GetN(); j++) lp->GetY()[j] = (pow(10, -lp->GetY()[j]/10));
  
  //loop over desired number of traces 
  for (int l = 0; l < num_traces; l++)
    {
      if (l % 50 == 0) printf("%d\n",l);
     //clear fftw array to fill w phasors
      for (int i = 0; i < 4996; i++){
        fftw_A[i].re = 0;
        fftw_B[i].re = 0;
       
        fftw_A[i].im = 0;
        fftw_B[i].im = 0;
      }   
      for(int iii = 0; iii < num_phasors; iii++){

        //generate on a spherical shell
        if (strcmp(geometry, "shell") == 0){
          int np = 0;
          while (np < 1){
            double x1 = tr1.Uniform(-1, 1);
            double x2 = tr1.Uniform(-1, 1);
	    phasor_phase = tr1.Uniform(0, 2 * TMath::Pi());
	    amp = amp_mean * sqrt(-2 * log(tr1.Uniform(0,1)));
            double d = pow(x1,2) + pow(x2, 2);
            if( d < 1 ){
	      phasor_x = radius_s * 2 * x1 * sqrt(1 - pow(x1,2) - pow(x2,2));
	      phasor_y = radius_s * 2 * x2 * sqrt(1 - pow(x1,2) - pow(x2,2));
	      phasor_z = radius_s * (1 - 2 * d);
              dA = sqrt(pow(phasor_x - ant_x, 2) + pow(phasor_y - ant_y, 2) + pow(phasor_z - ant_zA, 2));
              dB = sqrt(pow(phasor_x - ant_x, 2) + pow(phasor_y - ant_y,2) + pow(phasor_z - ant_zB, 2));
              rs = sqrt(pow(phasor_x, 2) + pow(phasor_y, 2) + pow(phasor_z, 2));
              theta = TMath::ACos(phasor_z/rs);
              amp = amp * fabs(TMath::Sin(theta));
              np = 1;
            }
          }
        }

        //generate from a box 
        if (strcmp(geometry, "box") == 0){
          int np = 0;
          while (np < 1){
            int wall_num = tr1.Integer(6); //pick a wall to generate from
	    phasor_phase = tr1.Uniform(0, 2 * TMath::Pi());
	    amp = amp_mean * sqrt(-2 * log(tr1.Uniform(0,1)));
            phasor_x = tr1.Uniform(-wall_size_x, wall_size_x);
            phasor_y = tr1.Uniform(-wall_size_y, wall_size_y);
            phasor_z = tr1.Uniform(-wall_size_z, wall_size_z);
            if (wall_num == 0) phasor_y = -wall_size_y;
            if (wall_num == 1) phasor_y = wall_size_y;
            if (wall_num == 2) phasor_x = -wall_size_x;
            if (wall_num == 3) phasor_x = wall_size_x;
            if (wall_num == 4) phasor_z = -wall_size_z;
            if (wall_num == 5) phasor_z = wall_size_z;
            dA = sqrt(pow(phasor_x - ant_x, 2) + pow(phasor_y - ant_y, 2) + pow(phasor_z - ant_zA, 2));
            dB = sqrt(pow(phasor_x - ant_x, 2) + pow(phasor_y - ant_y,2) + pow(phasor_z - ant_zB, 2));
            np = 1;
          }
        }


        //calculate phase and magnitude of phasors at each antenna and add in fourier space
        phase = phasor_phase;
        int i = tr1.Integer(4995);
        dAphase = 2 * TMath::Pi() * i * dF * dA/(3 * pow(10,8));
        dBphase = 2 * TMath::Pi() * i * dF * dB/(3 * pow(10,8));
        filter_gain = hp->Eval(i * dF/pow(10,6)) * lp->Eval(i * dF/pow(10, 6));
        filter_gain = strcmp(antenna_type, "vpol") ? filter_gain*ant_response->Eval(i * dF) : filter_gain; //don't have the antenna response for vpol, so just use filter response
        
        temp_A.setMagPhase(filter_gain * amp/dA, dAphase + phase);
        temp_B.setMagPhase(filter_gain * amp/dB, dBphase + phase);

        fftw_A[i].re += temp_A.re;
        fftw_A[i].im += temp_A.im;
        fftw_B[i].re += temp_B.re;
        fftw_B[i].im += temp_B.im;
        
        
      }
      
      //generate 75k amp noise
      for(int i = 0; i<4996; i++){
        filter_gain = hp->Eval(i * dF/pow(10,6)) * lp->Eval(i * dF/pow(10,6));
        amp_phase_A = tr1.Uniform(0,2*TMath::Pi());
        amp_mag_A = amp_noise * sqrt(-2 * log(tr1.Uniform(0,1)));
        amp_phase_B = tr1.Uniform(0,2*TMath::Pi());
        amp_mag_B = amp_noise * sqrt(-2 * log(tr1.Uniform(0,1)));
        temp_A.setMagPhase(filter_gain * amp_mag_A, amp_phase_A);
        temp_B.setMagPhase(filter_gain * amp_mag_B, amp_phase_B);
        fftw_A[i].re += temp_A.re;
        fftw_A[i].im += temp_A.im;
        fftw_B[i].re += temp_B.re;
        fftw_B[i].im += temp_B.im;
      }

      //inverse fft back for original waveforms

      yAfilt = FFTtools::doInvFFT(9991, fftw_A);
      yBfilt = FFTtools::doInvFFT(9991, fftw_B);

      g = new TGraph(9991, x_time, yAfilt);
      g2 = new TGraph(9991, x_time, yBfilt);
      
      file.cd();
      tree->Fill();
      file2.cd();
      tree2->Fill();

      delete g;
      delete g2;
    }      
  file.cd();
  tree->Write();
  file2.cd();
  tree2->Write();
}

void doAll()
{
 //these spacings are what I used for making the figure we used for GNO paper
  phasorModel("0_cm", 0);
  phasorModel("270_cm", 2.7);
  phasorModel("150_cm", 1.5);
  phasorModel("65_cm", .65);
  phasorModel("50_cm", .5);
  phasorModel("25_cm", .25);
  phasorModel("10_cm", .1);
  phasorModel("5_cm", .05);
  phasorModel("1_cm", .01);
  phasorModel("35_cm", .35);
  phasorModel("15_cm", .15);
  phasorModel("20_cm", .20);
  phasorModel("45_cm", .45);
  phasorModel("40_cm", .40);
  phasorModel("55_cm", .55);
  
  phasorModel("space0", 0.73025);
  phasorModel("space1", 0.7874);
  phasorModel("space2", 1.1303);
  phasorModel("space4", 0.9398);
  phasorModel("space3", 1.7653);
  
}
