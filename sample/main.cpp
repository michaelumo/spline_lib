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
  double div = 1000.0;
  // cspline
  Spline sp;
  sp.cspline(in);
  std::ofstream my_file("cspline.csv");
  for(double i = in(0,0); i < in(0,in.getCols()-1); i+=in(0,in.getCols()-1)/div){
    my_file<<i<<"  "<<sp.calc_y(i,in)<<"\n";
  }
  my_file.close();

  //bspline
  my_file.open("bspline.csv");
  int k = 3;//degree
  sp.bspline(in, k);
  for(double i = 0; i <= 1; i+=1.0/div){
    Matrix vec(2,1);
    sp.calc_vec(i, in, vec);
    my_file<<vec(0,0)<<"  "<<vec(1,0)<<"\n";
  }
  my_file.close();

  std::cout<<"DONE"<<std::endl;
  std::cout<<"Try: gnuplot> plot \"cspline.csv\" w l, \"bspline.csv\" w l"<<std::endl;
}
