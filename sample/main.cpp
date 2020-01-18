/*! @file sample/main.cpp
 *  @version 1.0.0
 *  @date Jan 18 2020
 *
 *  @brief
 *  This example shows how to create a cubic spline on R^3 and bspline on R^2
 *
 *  @Copyright (c) 2020 Michael O. Umenyi
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */


#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "matrix.h"
#include "matrix_algorithm.h"
#include "spline.h"

int main(void){
  Matrix in(4,4);

  //t
  in(0,0) = 0;
  in(0,1) = 0.5;
  in(0,2) = 1.0;
  in(0,3) = 1.5;

  //x
  in(1,0) = sin(in(0,0)*M_PI);
  in(1,1) = sin(in(0,1)*M_PI);
  in(1,2) = sin(in(0,2)*M_PI);
  in(1,3) = sin(in(0,3)*M_PI);

  //y
  in(2,0) = cos(in(0,0)*M_PI);
  in(2,1) = cos(in(0,1)*M_PI);
  in(2,2) = cos(in(0,2)*M_PI);
  in(2,3) = cos(in(0,3)*M_PI);

  //z
  in(3,0) = exp(in(0,0)*M_PI);
  in(3,1) = exp(in(0,1)*M_PI);
  in(3,2) = exp(in(0,2)*M_PI);
  in(3,3) = exp(in(0,3)*M_PI);

  double div = 100.0;

  // cspline
  Spline sp;
  sp.cspline(in);
  std::ofstream my_file("cspline.csv");
  for(double i = in(0,0); i < in(0,in.getCols()-1); i+=in(0,in.getCols()-1)/div){
    sp.calc_point(i, in);
    my_file<< i<<"  "<<sp.point[0]<<"  "<<sp.point[1]<<"  "<<sp.point[2]<<"  "<<sp.velocity[0]<<"  "<<sp.velocity[1]<<"  "<<sp.velocity[2]<<"  "<<sp.accel[0]<<"  "<<sp.accel[1]<<"  "<<sp.accel[2]<<"\n";
  }
  my_file.close();

  // bspline
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
  std::cout<<"Try: gnuplot> splot \"cspline.csv\" u 2:3:4 w l, \"cspline.csv\" u 5:6:7 w l, \"cspline.csv\" u 8:9:10 w l"<<std::endl;
  std::cout<<"     gnuplot> plot \"bspline.csv\" w l"<<std::endl;
}
