
/*
--
-- Copyright (C) 2014  <rmurufas@hotmail.com>
--
-- This program is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You may read the full text of the GNU General Public License
-- at <http://www.gnu.org/licenses/>.
--
*/



// autodiff.h
//
// C++ templates for Automatic Differentiation
//

#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <math.h>
#include <iostream>
using namespace std;

#include "vmv.h" //VectoMatrixVar templates




// interface routines

template <uint dim>
double deriv(const var<dim> a, const uint i)
{ return a.grad(i); }

template <uint dim>
vector<double,dim> gradient(const var<dim> a){ return a.grad; };


// gradient functions


template <uint dim>
var<dim> sqrt(const var<dim> a)
{
 double m;
 var<dim> y;
 y.val=sqrt(a.val);
 m=2.0*y.val;
 for(uint i=0; i<dim; i++)
 {
  if( a.grad(i)== 0.0 ) y.grad(i)=0.0;
  else y.grad(i)=a.grad(i)/m;
 }
 return y;
}


template <uint dim>
var<dim> exp(const var<dim> a)
{
 var<dim> y;
 y.val=exp(a.val);
 for(uint i=0; i<dim; i++)
 {
  if(a.grad(i)==0.0) y.grad(i)=0.0;
  else y.grad(i)=y.val*a.grad(i);
 }
 return y;
}


template <uint dim>
var<dim> ln(const var<dim> a)
{
 var<dim> y;
 y.val=log(a.val);
 for(uint i=0; i<dim; i++)
 {
  if(a.grad(i)==0.0) y.grad(i)=0.0;
  else y.grad(i)=a.grad(i)/a.val;
 }
 return y;
}

 
template <uint dim>
var<dim> sin(const var<dim> a)
{
 double m;
 var<dim> y;
 y.val=sin(a.val);
 m=cos(a.val);
 for(uint i=0; i<dim; i++)
 {
  if(a.grad(i)==0.0) y.grad(i)=0.0;
  else y.grad(i)=m*a.grad(i);
 }
 return y;
}

 
template <uint dim>
var<dim> cos(const var<dim> a)
{
 double m;
 var<dim> y;
 y.val=cos(a.val);
 m=-sin(a.val);
 for(uint i=0; i<dim; i++)
 {
  if(a.grad(i)==0.0) y.grad(i)=0.0;
  else y.grad(i)=m*a.grad(i);
 }
 return y;
}

 
template <uint dim>
var<dim> atan(const var<dim> a)
{
 double m;
 var<dim> y;
 y.val=atan(a.val);
 m = 1.0 + a.val * a.val;
 for(uint i=0; i<dim; i++)
 {
  if(a.grad(i)==0.0) y.grad(i)=0.0;
  else y.grad(i)=a.grad(i)/m;
 }
 return y;
}





template <uint dim>
double value(const var<dim> a){ return a.val; }

template <uint dim>
vector<double,dim> value(const vector<var<dim>,dim> x)
{
 vector<double,dim> b;
 for(uint i=0; i<dim; i++) b(i)=x(i).val;
 return b;
}

template <uint dim>
void set_indep_var(vector<var<dim>,dim>& x, const vector<double,dim> r)
{
 for(uint i=0; i<dim; i++)
 {
  x(i).val = r(i);
  for(uint j=0; j<dim; j++) x(i).grad(j)=0.0;
  x(i).grad(i) = 1.0;
 }
}

template <uint dim>
matrix<double,dim> jacobian(const vector<var<dim>,dim> x)
{
 matrix<double,dim> jac;
 for(uint i=0; i<dim; i++)
 {
  for(uint j=0; j<dim; j++)
  {
   jac(i,j) = x(i).grad(j);
  }
 }
 return jac;
}


template <uint dim>
var<dim> makvar( const double r )
{
 var<dim> y;
 y.val=r;
 for(uint i=0; i<dim; i++) y.grad(i)=0.0;
 return y;
}


// linear solver routines

// LU decomposition
template <uint dim>
void dec(matrix<double,dim>& a, vector<uint,dim>& ip, bool& singular)
{
 uint i1,k1,m,n1,n;
 double t;

 singular=false;
 n=dim-1;
 i1=0;
 ip(n)=i1;
  
 n1=n-1;
 for(uint k=i1; k<=n1; k++)
 {
  k1=k+1;
  m=k;
  for(uint i=k1; i<=n; i++)
  {
   if( fabs(a(i,k)) > fabs(a(m,k)) ) m=i;
  } // end for i
  ip(k)=m;
  t=a(m,k);
  if( m != k )
  {
   a(m,k)=a(k,k);
   a(k,k)=t;
  }

  if( t==0 )
  {
   singular=true;
   return;
   //break;
  }
  t=1.0/t;
  for(uint i=k1; i<=n; i++) a(i,k) *= -t;
  for(uint j=k1; j<=n; j++)
  {
   t=a(m,j);
   a(m,j)=a(k,j);
   a(k,j)=t;
   if( t != 0 ) for(uint i=k1; i<=n; i++) a(i,j) += a(i,k)*t;
  } // end for j
 } // end for k

}; // end dec


// Solve, RHS is overwritten with solution
template <uint dim>
void sol( const matrix<double,dim>& a, vector<double,dim>& b, const vector<uint,dim>& ip)
{
 uint i1,n,j,k1,m,n1;
 double t;

 n=dim-1;
 i1=0;
 n1=n-1;
 
 for(uint k=i1; k<=n1; k++)
 {
  k1=k+1;
  m=ip(k);
  t=b(m);
  b(m)=b(k);
  b(k)=t;
  for(uint i=k1; i<=n; i++) b(i) += a(i,k)*t;
 
 } // end for k
 
 k1=n;
 for(uint k=i1; k<=n1; k++)
 {
  j=k1;
  k1 -= 1;
  b(j) /= a(j,j);
  t= -b(j);
  for(uint i=i1; i<=k1; i++) b(i) += a(i,j)*t;
 } // end for k
 b(i1) /= a(i1,i1);

}; // end sol




// vector binary operators begin
template <uint dim>
vector<double,dim> operator+(const vector<double,dim> u, const vector<double,dim> v)
{
 vector<double,dim> w = u.plus(v);
 return w;
}

template <uint dim>
vector<double,dim> operator-(const vector<double,dim> u)
{
 vector<double,dim> w = u.unaryminus(); // -u;
 return w;
}

template <uint dim>
vector<double,dim> operator-(const vector<double,dim> u, const vector<double,dim> v)
{
 vector<double,dim> w = u.minus(v); // u-v;
 return w;
}

template <uint dim>
vector<double,dim> operator*(const double u, const vector<double,dim> v)
{
 vector<double,dim> w = v.times(u); // v*u;
 return w;
}

template <uint dim>
vector<double,dim> operator*(const vector<double,dim> u, const double v)
{
 vector<double,dim> w = u.times(v); // u*v;
 return w;
}

template <uint dim>
vector<double,dim> operator/(const vector<double,dim> u, const double v)
{
 vector<double,dim> w = u.divide(v); // u/v;
 return w;
}
// vector binary operators end



// gradient binary operators begin
template <uint dim>
var<dim> operator+(const var<dim> a, const var<dim> b)
{
 var<dim> y;
 y.val=a.val+b.val;
 for(uint i=0; i<dim; i++)
 {
  y.grad(i) = a.grad(i) + b.grad(i);
 }
 return y;
}

template <uint dim>
var<dim> operator+(const double a, const var<dim> b)
{
 return makvar<dim>(a)+b;
}

template <uint dim>
var<dim> operator+(const var<dim> a, const double b)
{
 return a+makvar<dim>(b);
}

template <uint dim>
var<dim> operator-(const var<dim> a, const var<dim> b)
{
 var<dim> y;
 y.val=a.val-b.val;
 for(uint i=0; i<dim; i++)
 {
  y.grad(i) = a.grad(i) - b.grad(i);
 }
 return y;
}

template <uint dim>
var<dim> operator-(const double a, const var<dim> b)
{ return makvar<dim>(a)-b; }

template <uint dim>
var<dim> operator-(const var<dim> a, double b)
{ return a-makvar<dim>(b); }

// unary minus of a var:
template <uint dim>
var<dim> operator-(const var<dim> a)
{ return 0.0-a; }




template <uint dim>
var<dim> operator*(const var<dim> a, const var<dim> b)
{
 var<dim> y;
 y.val=a.val*b.val;
 for(uint i=0; i<dim; i++)
 {
  y.grad(i) = a.grad(i)*b.val + b.grad(i)*a.val;
 }
 return y;
}

template <uint dim>
var<dim> operator*(const double a, const var<dim> b)
{ return makvar<dim>(a)*b; }

template <uint dim>
var<dim> operator*(const var<dim> a, const double b)
{ return a*makvar<dim>(b); }

template <uint dim>
var<dim> operator/(const var<dim> a, const var<dim> b)
{
 var<dim> y;
 y.val = a.val/b.val;
 for(uint i=0; i<dim; i++)
 {
  if( a.grad(i)==0.0 ) y.grad(i)=0.0;
  else y.grad(i)=a.grad(i)/b.val;
  if( b.grad(i) != 0.0 ) y.grad(i) -= b.grad(i)*y.val/b.val;
 }
 return y;
}

template <uint dim>
var<dim> operator/(const double a, const var<dim> b)
{ return makvar<dim>(a)/b; }

template <uint dim>
var<dim> operator/(const var<dim> a, const double b)
{ return a/makvar<dim>(b); }


// TBD:  define a power operator too

// gradient binary operators end


// more utilities:

template <uint dim>
double l2norm( const vector<double,dim> r )
{
 double sumsq(0.0);
 for(uint i=0; i<dim; i++) sumsq += r(i)*r(i);
 return sqrt(sumsq);
}



template <uint dim>
double l2norm( const vector<var<dim>,dim> y )
{
 double sumsq(0.0);
 vector<double,dim> r=value<dim>(y);
 for(uint i=0; i<dim; i++) sumsq += r(i)*r(i);
 return sqrt(sumsq);
}




template <uint dim>
double maxnorm( const vector<double,dim> r )
{
 double max(0.0);
 for(uint i=0; i<dim; i++)
  if(fabs(r(i))>max) max=fabs(r(i));
 return max;
}



template <uint dim>
double maxnorm( const vector<var<dim>,dim> y )
{
 double max(0.0);
 vector<double,dim> r;
 r=value<dim>(y);
 for(uint i=0; i<dim; i++)
  if(fabs(r(i))>max) max=fabs(r(i));
 return max;
}


template <uint dim>
vector<var<dim>,dim> operator+( const vector<var<dim>,dim> y1, const vector<var<dim>,dim> y2 )
{
 vector<var<dim>,dim> z;
 for(uint i=0; i<dim; i++) z(i) = y1(i) + y2(i);
 return z;
}



template <uint dim>
vector<var<dim>,dim> operator-( const vector<var<dim>,dim> y1, const vector<var<dim>,dim> y2 )
{
 vector<var<dim>,dim> z;
 for(uint i=0; i<dim; i++) z(i) = y1(i) - y2(i);
 return z;
}


#endif

