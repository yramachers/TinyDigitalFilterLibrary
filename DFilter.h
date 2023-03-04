#ifndef YR_DFilter
#define YR_DFilter
// Filter collection library
// Changelog, YR:
// March 2023: add Butterworth filter
// YR, University of Warwick 2019

#include <vector>
#include <array>

class DFButterworth {
  // Low-pass Butterworth filter, any order
  // Follows https://github.com/adis300/filter-c

 private:
  double ftimebase = 1.e-9; // unit nano seconds
  double flowfreq  = 0.0;
  int    fOrder    = 2;

 protected:

 public:
  DFButterworth() = default;
  virtual   ~DFButterworth() = default;
  std::vector<double>  Filter(std::vector<double> &record);

  // setter and getter
  void      SetSamplingTimeBase(double ff); // unit nano seconds

  void      SetLowFilterFreq(double low);
  double    GetLowFilterFreq() {return flowfreq;}

  void      SetFilterOrder(int o);
  int       GetFilterOrder() {return fOrder;}

};


class DFMatched {

 private:
  std::vector<double> fstde;

 protected:
  void correlate(std::vector<double> &data, std::vector<double> &result);
  std::vector<double> rightPadding(std::vector<double> &record, unsigned int width);

 public:
  DFMatched() = default;
  virtual  ~DFMatched(); // clears container
  std::vector<double>   Filter(std::vector<double> &record);

  // setter and getter
  void SetStdEvent(std::vector<double> &s);

};


class DFLowRCfilter {

 private:
  double ftimebase     = 1.e-9; // unit nano seconds
  double flowfreq      = 0.0;
  int    fNumberPasses = 1;
  std::array<double, 2> fresponse; // fixed size

 protected:

 public:
  DFLowRCfilter() = default;
  virtual   ~DFLowRCfilter() = default;
  std::vector<double>  Filter(std::vector<double> &record);

  // setter and getter
  void      SetSamplingTimeBase(double ff);

  void      SetLowRCfilterFreq(double low);
  double    GetLowRCfilterFreq() {return flowfreq;}

  void      SetLowRCfilterNumberPasses(int npass);
  int       GetLowRCfilterNumberPasses() {return fNumberPasses;}

};

class DFHighRCfilter {

 private:
  double ftimebase     = 1.e-9; // unit nano seconds
  double fhighfreq     = 0.0;
  int    fNumberPasses = 1;
  std::array<double, 3> fresponse; // fixed size

 protected:

 public:
  DFHighRCfilter() = default;
  virtual   ~DFHighRCfilter() = default;
  std::vector<double>  Filter(std::vector<double> &record);

  // setter and getter
  void      SetSamplingTimeBase(double ff); // unit nano seconds
  void      SetHighRCfilterFreq(double high);
  double    GetHighRCfilterFreq() {return fhighfreq;}
  void      SetHighRCfilterNumberPasses(int npass);
  int       GetHighRCfilterNumberPasses() {return fNumberPasses;}

};

class DFMovingAverage {

 private:
  int      fMAwidth = 5;        // Moving Average filter width, default=5
  double   fresponse;

 protected:
  std::vector<double> paddingData(std::vector<double> &record, int width);

 public:
  DFMovingAverage() = default;
  virtual  ~DFMovingAverage() = default;
  std::vector<double> Filter(std::vector<double> &record);

  // setter and getter
  void     SetMovingAverageWidth(int width);          // Manipulate the
  int      GetMovingAverageWidth(){return fMAwidth;}  // window width

};

#endif
