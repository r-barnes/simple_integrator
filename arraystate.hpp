#ifndef __array_state_
#define __array_state_

#include <array>
#include <cmath>
#include <vector>

template <class T, int N>
class arraystate : public std::array<T,N> {
  public:
    arraystate() {}
    arraystate(const std::vector<T> &init){
      for(int i=0;i<std::array<T,N>::size();++i)
        std::array<T,N>::operator[](i)=init[i];
    }

    arraystate<T,N>& operator+=(const arraystate<T,N> &other) {
      for(int i=0;i<std::array<T,N>::size();++i)
        std::array<T,N>::operator[](i)+=other[i];
    }
};

template <class T, int N>
arraystate<T,N> operator+( const arraystate<T,N> &a , double b ){
  arraystate<T,N> result(a);
  for(double& i: result)
    i+=b;
  return result;
}

template <class T, int N>
arraystate<T,N> operator+( double a, const arraystate<T,N> &b ){
  return b+a;
}

template <class T, int N>
arraystate<T,N> operator+( const arraystate<T,N> &a, const arraystate<T,N> &b ){
  arraystate<T,N> result(a);
  for(unsigned int i=0;i<result.size();++i)
    result[i]+=b[i];
  return result;
}

template <class T, int N>
arraystate<T,N> operator*( const arraystate<T,N> &a, double b ){
  arraystate<T,N> result(a);
  for(double& i: result)
    i*=b;
  return result;
}

template <class T, int N>
arraystate<T,N> operator*( double a, const arraystate<T,N> &b ){
  return b*a;
}

template <class T, int N>
arraystate<T,N> operator/( const arraystate<T,N> &a, double b){
  arraystate<T,N> result(a);
  for(double& i: result)
    i/=b;
  return result;
}

template <class T, int N>
double abs( const arraystate<T,N> &a ){
  double result=0;
  for(const double& i: a)
    result+=std::abs(i);
  return result;
}

template <class T, int N>
std::ostream& operator<<(std::ostream &out, const arraystate<T,N> &a){
  for(typename arraystate<T,N>::const_iterator i=a.begin();i!=a.end();++i)
    out<<*i<<" ";
  return out;
}


#endif
