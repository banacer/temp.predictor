// normal_distribution
#include <iostream>
#include <random>

using namespace std;
int main()
{
  const int nrolls=120;  // number of experiments
  const int nstars=100;    // maximum number of stars to distribute

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.05,0.02);

  int p[10]={};
  double start = 75.05;
  for (int i=0; i<nrolls; ++i) {
    double number = distribution(generator);
    start += number;
    cout << start << endl;
  }

  return 0;
}