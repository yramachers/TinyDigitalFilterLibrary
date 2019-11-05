#ifndef YR_DFilter
#define YR_DFilter
// Filter collection library
// YR, University of Warwick 2019

#include <vector>
#include <array>

class DFMatched {

 private:
  std::vector<double> fstde;

 protected:
  void correlate(std::vector<double> &data, std::vector<double> &result);
  std::vector<double> rightPadding(std::vector<double> &record, unsigned int width);

 public:
  DFMatched();
  virtual  ~DFMatched();
  std::vector<double>   Filter(std::vector<double> &record);

  // setter and getter
  void SetStdEvent(std::vector<double> &s);

};


class DFLowRCfilter {

 private:
  double ftimebase;
  double flowfreq;
  int   fNumberPasses;
  std::array<double, 2> fresponse; // fixed size

 protected:

 public:
  DFLowRCfilter();
  virtual   ~DFLowRCfilter();
  std::vector<double>  Filter(std::vector<double> &record);

  // setter and getter
  void      SetSamplingTimeBase(double ff); // unit nano seconds

  void      SetLowRCfilterFreq(double low);
  double    GetLowRCfilterFreq() {return flowfreq;}

  void      SetLowRCfilterNumberPasses(int npass);
  int       GetLowRCfilterNumberPasses() {return fNumberPasses;}

};

class DFHighRCfilter {

 private:
  double ftimebase;
  double fhighfreq;
  int   fNumberPasses;
  std::array<double, 3> fresponse; // fixed size

 protected:

 public:
  DFHighRCfilter();
  virtual   ~DFHighRCfilter();
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
  int      fMAwidth;        // Moving Average filter width, default=5
  double   fresponse;

 protected:
  std::vector<double> paddingData(std::vector<double> &record, int width);

 public:
  DFMovingAverage();
  virtual  ~DFMovingAverage();
  std::vector<double> Filter(std::vector<double> &record);

  // setter and getter
  void     SetMovingAverageWidth(int width);          // Manipulate the
  int      GetMovingAverageWidth(){return fMAwidth;}  // window width

};

#endif
