#include <DFilter.h>
#include <iostream>
#include <vector>
#include <numeric>

// ROOT
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TRandom.h>

int main() {
  void datasource(std::vector<double> &t, std::vector<double> &d);
  void oscillationsource(std::vector<double> &d);
  TH1D* testbworth(std::vector<double>& data);  
  TH1D* testmav(std::vector<double>& data);
  TH1D* testLP(std::vector<double>& data);
  TH1D* testHP(std::vector<double>& data);
  TH1D* testmatch(std::vector<double>& data, std::vector<double>& st);
  TH1D* testQGshaper(std::vector<double>& data);

  // Start
  // output file creation
  TFile* ff = new TFile("histos.root","recreate");

  std::vector<double> stde;
  std::vector<double> data;
  std::vector<double> osc;
  datasource(stde, data);
  oscillationsource(osc);

  std::cout << "Made pulse, ndata = " << data.size() << std::endl;
  std::cout << "Made stde size = " << stde.size() << std::endl;
  std::cout << "Made oscillation size = " << osc.size() << std::endl;

  // test moving average filter
  TH1D* bwhist = testbworth(osc);

  // test moving average filter
  TH1D* mavhist = testmav(data);

  // test low pass filter
  TH1D* lrchist = testLP(data);

  // test high pass filter
  TH1D* hrchist = testHP(data);

  // test matched filter
  TH1D* mhist = testmatch(data,stde);

  // test matched filter
  TH1D* qghist = testQGshaper(data);

  // output
  int bin = 0;
  TH1D* dhist = new TH1D("data","data",data.size(),0,data.size()-1);
  for (double& entry : data) {
    dhist->SetBinContent(bin, entry);
    bin++;
  }
  bin = 0;
  TH1D* ohist = new TH1D("oscd","oscillation",osc.size(),0,osc.size()-1);
  for (double& entry : osc) {
    ohist->SetBinContent(bin, entry);
    bin++;
  }

  dhist->Write();
  ohist->Write();
  bwhist->Write();
  mavhist->Write();
  lrchist->Write();
  hrchist->Write();
  mhist->Write();
  qghist->Write();
  ff->Close();
  return 0;
}

TH1D* testbworth(std::vector<double>& data) {
  DFButterworth lrc;
  std::cout << "Instantiated DFButterworth." << std::endl;

  // single pass filter is default
  lrc.SetSamplingTimeBase(1.e-3); // 1 kHz; default, 1 ns
  lrc.SetLowFilterFreq(1.0e2);  // 100 Hz

  std::vector<double> newosc = lrc.Filter(data);  
  std::cout << "Butterworth executed" << std::endl;

  int nf = newosc.size();
  TH1D* fhist = new TH1D("bwf","filtered pulse",nf,0,nf-1);
  for (int i=0;i<nf;i++) 
    fhist->SetBinContent(i,newosc[i]);
  return fhist;
}


TH1D* testmav(std::vector<double>& data) {
  DFMovingAverage mav;
  std::cout << "Instantiated DFMovingAverage." << std::endl;

  mav.SetMovingAverageWidth(7);

  std::vector<double> newpulse = mav.Filter(data);  
  std::cout << "MAVFilter executed" << std::endl;

  int nf = newpulse.size();
  TH1D* fhist = new TH1D("mav","filtered pulse",nf,0,nf-1);
  for (int i=0;i<nf;i++) 
    fhist->SetBinContent(i,newpulse.at(i));
  return fhist;
}


TH1D* testLP(std::vector<double>& data) {
  DFLowRCfilter lrc;
  std::cout << "Instantiated DFLowRCfilter." << std::endl;

  // single pass filter is default
  lrc.SetSamplingTimeBase(1.0); // 1 ns
  lrc.SetLowRCfilterFreq(1.0e6); // 1 MHz

  std::vector<double> newpulse = lrc.Filter(data);  
  std::cout << "LPFilter executed" << std::endl;

  int nf = newpulse.size();
  TH1D* fhist = new TH1D("lrc","filtered pulse",nf,0,nf-1);
  for (int i=0;i<nf;i++) 
    fhist->SetBinContent(i,newpulse[i]);
  return fhist;
}


TH1D* testHP(std::vector<double>& data) {
  DFHighRCfilter hrc;
  std::cout << "Instantiated DFHighRCfilter." << std::endl;

  // single pass filter is default
  hrc.SetSamplingTimeBase(1.0); // 1 ns
  hrc.SetHighRCfilterFreq(1.0e6); // 1 MHz

  std::vector<double> newpulse = hrc.Filter(data);  
  std::cout << "HPFilter executed" << std::endl;

  int nf = newpulse.size();
  TH1D* fhist = new TH1D("hrc","filtered pulse",nf,0,nf);
  for (int i=0;i<nf;i++) 
    fhist->SetBinContent(i,newpulse[i]);
  return fhist;
}


TH1D* testQGshaper(std::vector<double>& data) {
  DFHighRCfilter hrc;
  std::cout << "Instantiated DFHighRCfilter." << std::endl;
  DFLowRCfilter lrc;
  std::cout << "Instantiated DFLowRCfilter." << std::endl;

  // Quasi Gaussian shaping from HP-4LP series filtering
  // single pass filter is default
  hrc.SetSamplingTimeBase(1.0); // 1 ns
  hrc.SetHighRCfilterFreq(2.0e6); // 2 MHz
  // 4-pass filter is for quasi-Gaussian shaping
  lrc.SetLowRCfilterNumberPasses(4); // 4 passes
  lrc.SetSamplingTimeBase(1.0); // 1 ns
  lrc.SetLowRCfilterFreq(1.0e7); // 10 MHz


  std::vector<double> stage1 = hrc.Filter(data);
  std::vector<double> newpulse = lrc.Filter(stage1);
  std::cout << "QGFilter executed" << std::endl;

  int nf = newpulse.size();
  TH1D* fhist = new TH1D("qgf","filtered pulse",nf,0,nf);
  for (int i=0;i<nf;i++) 
    fhist->SetBinContent(i,newpulse[i]);
  return fhist;
}


TH1D* testmatch(std::vector<double>& data, std::vector<double>& st) {
  DFMatched mat;
  std::cout << "Instantiated DFMatched." << std::endl;

  mat.SetStdEvent(st);

  std::vector<double> newpulse = mat.Filter(data);  
  std::cout << "MATCHFilter executed" << std::endl;

  int nf = newpulse.size();
  TH1D* fhist = new TH1D("matched","matched filtered pulse",nf,0,nf-1);
  for (int i=0;i<nf;i++) 
    fhist->SetBinContent(i,newpulse[i]);
  return fhist;
}


void oscillationsource(std::vector<double> &osc) {
  int samplelength = 10000; // 10 s at 1 kHz
  double f[4] = {25,110,150,250}; // frequencies [Hz]
  const double tpi = 2.0*std::acos(-1.0);
  for (int i=0;i<samplelength;++i)
    osc.push_back(std::sin(tpi*f[0]*i*10/samplelength)
		  +std::sin(tpi*f[1]*i*10/samplelength)
		  +std::sin(tpi*f[2]*i*10/samplelength)
		  +std::sin(tpi*f[3]*i*10/samplelength));
}


void datasource(std::vector<double> &stde, std::vector<double> &d) {
  double mypulse(double x, double *par);
  int pulselength = 1024;
  int stdlength  = 128;
  double par[4] = {1.0, 2.0, 100.0, 128.0}; // amplitude, risetime, decaytime, onset

  TRandom *rand = new TRandom();
  rand->SetSeed(0);

  // data with noise added, level 20 percent amplitude
  double nlevel = 0.2 * par[0];
  for (int i=0;i<pulselength;++i) d.push_back(mypulse( (double)i, par) + nlevel * (rand->Gaus(0,1))); // fill the data vector

  // standard event to look for, normalize for amplitude preservation with matched filter
  int onset = (int)par[3];
  for (int i=onset;i<stdlength+onset;++i) stde.push_back(mypulse( (double)i, par)); // fill the stdevent vector
  double cumsum = std::accumulate(stde.begin(), stde.end(), 0.0); // total sum
  for (unsigned int i=0;i<stde.size();++i) stde.at(i) /= cumsum; // normalize

}


double mypulse(double x, double *par)
{
  if (x < par[3]) return 0.0; // baseline fixed at 0.0
  double f = par[0]*(TMath::Exp(-(x - par[3])/par[1]) - TMath::Exp(-(x - par[3])/par[2]));
  return -f; // rising pulse, baseline at 0
}
