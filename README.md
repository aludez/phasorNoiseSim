Usage in command line:

root
.L toyModelImproved.C+

from here you can call the function phasorModel(*args).  phasorModel outputs two .root files containing TTrees filled with traces generated from the number of phasors and geometry specified.

The inputs to phasorModel are:
dist_name - a string that gets appeded to the output root files
  ant_sep - the separation between antenna feeds
save - whether you want to save the TTrees or not (1 yes 0 no)
  phasor_num - the number of phasors that generate each trace
  geometry - a string ("shell" or "box") that determines the layout where the phasors come from
  traces - number of waveforms that go into the TTree
  antenna_type - a string ("vpol", "lba" or "hba") that picks what antenna frequency response to use. "vpol" is ara vpol.
temp - noise temperature of the amplifier (300 noise temp outside is currently hardcoded in)


  The outputs of phasorModel are:
  two .root files containing TTrees filled with TGraphs of the traces generated by this simulation and a third .root file with histograms of the mean peak correlation coefficients

  At the bottom I put a function called doAll() with some phasorModel called a bunch of times with the settings I used to make the correlation vs spacing plot in the paper

  there are also 2 macros to draw graphs. mpcFig and mpcFig2Types, mpcFig draws one set of simulation points, the other overlays 2 

  The files high_pass.txt, low_pass2.txt, H.csv and L.csv are included because they are referenced in the toyPhasorModel.  H and L are the HBA and LBA frequency responses, and high_pass and low_pass2 are the high and low pass filters responses we used in the anechoic chamber tests.

  I've also added a compilable, executable version (named noiseModelApp) as long as you have libRootFftwWrapper installed.  Inputs are the same as above but you can compile it with make and run it from the command line.
# phasorNoiseSim
