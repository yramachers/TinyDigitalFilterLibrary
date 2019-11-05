#include "catch.hpp"
#include <DFilter.h>
#include <vector>
#include <iostream>


std::vector<double> squarePulse() 
{
  double d[] = {0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0}; // square pulse, 3 wide
  std::vector<double> data (d, d + sizeof(d) / sizeof(double));
  return data;
}

int check_pulseA()
{
  std::vector<double> target (3, 1.0/3.0); // 3 doubles, value 1.0, normed
  std::vector<double> data = squarePulse();

  DFMatched mf;
  mf.SetStdEvent(target);
  std::vector<double> filt = mf.Filter(data);

  std::cout << "Sq Pulse Response: " << std::endl;
  for (auto& value : filt) std::cout << value << ", ";
  std::cout << std::endl;

  return (int)filt.at(4); // middle should be 1
}


unsigned int check_pulseLength()
{
  std::vector<double> target (3, 1.0/3.0); // 3 doubles, value 1.0, normed
  std::vector<double> data = squarePulse();

  DFMatched mf;
  mf.SetStdEvent(target);
  std::vector<double> filt = mf.Filter(data);

  return (unsigned int)filt.size(); // middle should be 1
}


TEST_CASE( "Pulse A", "[pulsetest]" ) {
  REQUIRE( check_pulseA() == 1 );
}

TEST_CASE( "Pulse B", "[pulsetest2]" ) {
  REQUIRE( check_pulseLength() == 11 );
}

