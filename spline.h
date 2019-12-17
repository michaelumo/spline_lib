#ifndef __SPLINE_H__
#define __SPLINE_H__

class Spline{
  private:
    //cspline variables
    std::vector<Matrix> f;
    std::vector<double> h;
    // bspline variables
    std::vector<double> U;//Open Uniform Knot Vector
    int degree = 0;
    double B(int i, int k, double u);
  public:
    std::vector<float> point;
    Matrix pos, vel, accel;
    Spline();
    void cspline(Matrix &in);
    void bspline(Matrix &in, int k);
    void calc_point(double i, Matrix &in);
    void calc_vec(double u, Matrix &in, Matrix &vec);
};

Spline::Spline(){
  ;
}

void Spline::cspline(Matrix &in){//type of matrix is (t,x,y,z..).t()
  f.clear();
  int num = in.getRows();
  TDMA alg;
  point = std::vector<float>(num);
  for(int n = 1; n < num; n++){
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
      D.push_back(6.0*((in(n,i+1)-in(n,i))/h[i])-6.0*((in(n,i)-in(n,i-1))/h[i-1]));
    }
    C.push_back(0);
    f.push_back(alg.tdma(A,B,C,D));
  }
}

void Spline::calc_point(double i, Matrix &in){
  int j = 0;
  for(int k = 0; k < in.getCols()-1; k++){
    if(in(0,k) <= i && i<=in(0,k+1)){
      j = k;
      break;
    }
  }
  for(int n = 0; n < in.getRows()-1; n++){
    point.insert(point.begin()+n, f[n](j+1,0)/(6.0*h[j])*pow((i-in(0,j)),3.0)+f[n](j,0)/(6.0*h[j])*pow((in(0,j+1)-i),3.0)+(in(n+1,j+1)/h[j]-f[n](j+1,0)/6.0*h[j])*(i-in(0,j))+(in(n+1,j)/h[j]-f[n](j,0)/6.0*h[j])*(in(0,j+1)-i));
  }
}

void Spline::bspline(Matrix &in, int k){//k is degree, k < in.getCols()
  degree = k;
  if(k >= in.getCols()){std::cout<<"Error: degree is to high! Must be k < "<<in.getCols()<<std::endl; return;}
  int n = degree+in.getCols()+1;//number of knots
  for(int i = 0; i < degree; i++){
    U.push_back(0.0);
  }
  for(int i = 0; i <n-2*degree; i++){
    U.push_back(1.0/((double)n-2*degree-1)*i);
  }
  for(int i = 0; i < degree; i++){
    U.push_back(1.0);
  }
}

void Spline::calc_vec(double u, Matrix &in, Matrix &vec){
  for(int j = 0; j < in.getCols(); j++){
    vec(0,0)+=in(0,j)*B(j,degree,u);
    vec(1,0)+=in(1,j)*B(j,degree,u);
  }
}

double Spline::B(int i, int k, double u){// k: degree, i: index of control points, u: knot
  if(k == 0 && u >= U[i] && u < U[i+1])return 1.0;
  else if(k == 0)return 0.;
  double b1, b2;
  if((U[i+k]-U[i]) == 0)b1 = 0.;
  else b1 = (u-U[i])/(U[i+k]-U[i])*B(i,k-1,u);
  if((U[i+k+1]-U[i+1]) == 0)b2 = 0.;
  else b2 = (U[i+k+1]-u)/(U[i+k+1]-U[i+1])*B(i+1,k-1,u);
  return b1+b2;
}

#endif
