#ifndef __SPLINE_H__
#define __SPLINE_H__

class Spline{
  private:
    Matrix f;
    std::vector<double> h;
  public:
    Spline();
    void cspline(Matrix &in);
    double calc_y(double i, Matrix &in);
};

Spline::Spline(){
  ;
}

void Spline::cspline(Matrix &in){//type of matrix is (x,y).t()
	std::vector<double> A;
  std::vector<double> B;
  std::vector<double> C;
  std::vector<double> D;
  h.clear();
  //calculate delta t
  for(int i = 0; i < in.getCols()-1; i++){
    h.push_back(in(0,i+1)-in(0,i));
  }
  A.push_back(0);
  B.push_back(0);
  C.push_back(0);
  D.push_back(0);
  A.push_back(0);
  for(int i = 1; i < h.size(); i++){
    B.push_back(2*(h[i-1]+h[i]));
    if(i != h.size()-1)C.push_back(h[i]);
    if(i>1)A.push_back(h[i-1]);
    D.push_back(6.0*((in(1,i+1)-in(1,i))/h[i])-6.0*((in(1,i)-in(1,i-1))/h[i-1]));
  }
  C.push_back(0);
  f = in.TDMA(A,B,C,D);
}

double Spline::calc_y(double i, Matrix &in){
  int j = 0;
  for(int k = 0; k < in.getCols()-1; k++){
    if(in(0,k) <= i && i<=in(0,k+1)){
      j = k;
      break;
    }
  }
  return f(j+1,0)/(6.0*h[j])*pow((i-in(0,j)),3.0)+f(j,0)/(6.0*h[j])*pow((in(0,j+1)-i),3.0)+(in(1,j+1)/h[j]-f(j+1,0)/6.0*h[j])*(i-in(0,j))+(in(1,j)/h[j]-f(j,0)/6.0*h[j])*(in(0,j+1)-i);
}
#endif
