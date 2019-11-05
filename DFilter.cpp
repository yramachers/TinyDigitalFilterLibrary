#include <DFilter.h>
#include <iostream>
#include <cmath>


DFMatched::~DFMatched()
{
  fstde.clear();
}


DFMatched::DFMatched() : fstde(0)
{
}


void DFMatched::SetStdEvent(std::vector<double> &s)
{
  fstde.insert(fstde.begin(), s.begin(), s.end()); // straight copy
}


std::vector<double> DFMatched::Filter(std::vector<double> &record)
{
  std::vector<double> result;

  if (!fstde.empty()) {

    if (!record.empty())
      correlate(record, result); // work happens here
    else
      std::cout << "Error data not valid in filter" << std::endl;
  }
  else
    std::cout << "Error std event empty " << std::endl;
  
  return result; // empty return if clauses fail
}


void DFMatched::correlate(std::vector<double> &data, std::vector<double> &result)
{
  result.clear(); // clear in any case
  std::vector<double> padded = rightPadding(data, fstde.size()); // zero padding at the right for target overlap

  // cross correlation loops - apparently more efficient for small target on bigger data array
  for (unsigned int i=0;i<data.size();++i) {
    double sum = 0.0;
    for (unsigned int j=0;j<fstde.size();++j) sum += padded.at(i+j) * fstde.at(j);
    result.push_back(sum);
    sum = 0.0;
  }
}


std::vector<double> DFMatched::rightPadding(std::vector<double> &record, unsigned int width)
{
  std::vector<double> padded (record); // copy
  for (unsigned int i=0;i<width;++i) padded.push_back(0.0); // zero padding
  return padded;
}


// RC Low-pass filter
DFLowRCfilter::~DFLowRCfilter()
{
}


DFLowRCfilter::DFLowRCfilter() : ftimebase(1),
			     flowfreq(0), fNumberPasses(1)
{
}


void DFLowRCfilter::SetSamplingTimeBase(double ff)
{
  if (ff > 0.0) ftimebase = ff; // take only positive time base
  else ftimebase = 1.0; // unit nano seconds
}


void DFLowRCfilter::SetLowRCfilterFreq(double low)
{
  if (low >= 0.0) flowfreq = low; // take only positive frequencies
  else flowfreq = 0.0;
}


void DFLowRCfilter::SetLowRCfilterNumberPasses(int npass)
{
  if (npass > 1) fNumberPasses = npass;
  else fNumberPasses = 1;        // Default is single pass filter
}


std::vector<double> DFLowRCfilter::Filter(std::vector<double> &record)
{
  std::vector<double> result;

  double pi = std::acos(-1.0);
  double fc, nyfreq, unit;
  
  unit = 1.0e-9; // nano seconds
  nyfreq = 1.0 / (2.0 * ftimebase * unit);  // Nyquist frequency
  if (flowfreq < nyfreq) // must be smaller than Nyquist frequency
    fc = std::exp(- 2.0 * pi * flowfreq / nyfreq);
  else 
    fc = 0.0;           // otherwise no filtering on low frequencies

  // filter coefficients
  fresponse[0] = 1.0 - fc;
  fresponse[1] = fc;

  if (!record.empty()) {
    
    std::vector<double> copydata(record);  // dummy storage for the data
    result.resize(record.size());
    
    for (int j=1;j<=fNumberPasses;++j) { 
      result.at(0) = fresponse[0] * copydata.at(0); // fast recursive low pass
      for (unsigned int i=1; i<record.size(); ++i) 
	result.at(i) = fresponse[0] * copydata.at(i) + fresponse[1] * result.at(i-1);
      copydata = result; // for npass>1 copydata takes 
      // filtered data and is filtered again
    }
  }
  else
    std::cout << "Error data not valid in apply filter" << std::endl;
  
  return result;
}


// RC High-pass filter
DFHighRCfilter::~DFHighRCfilter()
{
}


DFHighRCfilter::DFHighRCfilter() : ftimebase(1),
			     fhighfreq(0), fNumberPasses(1)
{
}


void DFHighRCfilter::SetSamplingTimeBase(double ff)
{
  if (ff > 0.0) ftimebase = ff; // take only positive time base
  else ftimebase = 1.0; // unit nano seconds
}


void DFHighRCfilter::SetHighRCfilterFreq(double high)
{
  if (high >= 0.0) fhighfreq = high; // take only positive frequencies
  else fhighfreq = 0.0;
}


void DFHighRCfilter::SetHighRCfilterNumberPasses(int npass)
{
  if (npass > 1) fNumberPasses = npass;
  else fNumberPasses = 1;        // Default is single pass filter
}


std::vector<double> DFHighRCfilter::Filter(std::vector<double> &record)
{
  std::vector<double> result;
  
  double pi = std::acos(-1.0);
  double fc, nyfreq, unit;
  
  unit = 1.0e-9; // nano seconds
  nyfreq = 1.0 / (2.0 * ftimebase * unit);  // Nyquist frequency
  
  if (fhighfreq < nyfreq) // must be smaller than Nyquist frequency
    fc = std::exp(- 2.0 * pi * fhighfreq / nyfreq);
  else 
    fc = 0.0;           // otherwise no filtering on high frequencies
  
  // filter coefficients
  fresponse[0] = 0.5 * (1.0 + fc);
  fresponse[1] = - fresponse[0];
  fresponse[2] = fc;

  if (!record.empty()) {
    
    std::vector<double> copydata(record);  // dummy storage for the data
    result.resize(record.size());
    
    for (int j=1;j<=fNumberPasses;++j) { 
      result.at(0) = fresponse[0] * copydata.at(0); // fast recursive high pass
      for (unsigned int i=1; i<record.size(); ++i) 
	result.at(i) = fresponse[0] * copydata.at(i) + fresponse[1] * copydata.at(i-1)
	  + fresponse[2] * result.at(i-1);
      copydata = result; // for npass>1 copydata takes 
      // filtered data and is filtered again
    }
  }
  else
    std::cout << "Error data not valid in apply filter" << std::endl;
  
  return result;
}



DFMovingAverage::~DFMovingAverage()
{
}


DFMovingAverage::DFMovingAverage() : fMAwidth(5),
				     fresponse(0.0)
{
}


void DFMovingAverage::SetMovingAverageWidth(int width)
{
  if (width%2) fMAwidth = (width>3) ? width : 3 ;    // Minimum width here is 3
  else fMAwidth = (width>3) ? (width + 1) : 3 ; // take only odd integer widths
}


std::vector<double> DFMovingAverage::Filter(std::vector<double> &record)
{
  std::vector<double> result;
  
  if (fMAwidth>0.0)
    fresponse = 1.0 / fMAwidth;
  else
    return result; // empty return

  if (!record.empty()) {
    
    if (fMAwidth>=(0.5 * record.size())) { // window too big
      return result; // empty return
    }

    // Data has been set and the response calculated; so apply it
    
    double sum = 0.0;
    int start = (int)(std::floor(0.5 * fMAwidth));
    std::vector<double> padded = paddingData(record, start);
    
    for (unsigned int j=0; j<record.size(); ++j) { // all entries in record

      for (int i=-start; i<=start; ++i)	sum += padded.at(start + i + j);

      result.push_back(sum * fresponse);
      sum = 0.0;
    }
  }
    
  return result;
}

std::vector<double> DFMovingAverage::paddingData(std::vector<double> &record, int width)
{
  std::vector<double> padded;

  for (int i=width;i>0;--i) padded.push_back(record.at(i)); // flipped data padding left of start
  padded.insert(padded.end(), record.begin(), record.end()); // straight copy
  for (int i=1;i<=width;++i) padded.push_back(record.at(record.size() - i)); // flipped data padding right of end

  return padded;
}
