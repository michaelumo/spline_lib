#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "matrix.h"
#include "spline.h"

int main(void){
  Matrix in(2,4);
  in(0,0) = 1;
  in(1,0) = 400.0;
  in(0,1) = 10.0;
  in(1,1) = 800;
  in(0,2) = 15.0;
  in(1,2) = 400.0;
  in(0,3) = 20.0;
  in(1,3) = 100.0;

  Spline sp;
  sp.cspline(in);
  std::ofstream my_file("output.csv");
  for(double i = in(0,0); i < in(0,in.getCols()-1); i+=in(0,in.getCols()-1)/1000.0){
    my_file<<i<<"  "<<sp.calc_y(i,in)<<"\n";
  }
  std::cout<<"DONE"<<std::endl;
  std::cout<<"Try: gnuplot> plot \"output.csv\" w l"<<std::endl;
}
