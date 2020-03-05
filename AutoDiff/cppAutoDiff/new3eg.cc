
/*
--
-- Copyright (C) 2020  <rmurufas@hotmail.com>
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


// new3eg.cc
//
// demonstrates automatic differentiation to evaluate the
// Jacobian of a vector valued function of a vector variable.  
//
// This particular example shows the usage of damped Newton's 
// method for finding the root of a nonlinear function:
//
// find x,y such that
//
// F(x,y) = 0 where
//
// F(x,y) = atan2(y,x) - 1.0
//          x+y-1.0
//

#include <stdlib.h> // exit
#include <math.h>
#include <iostream>
//using namespace std;
using std::cout;
using std::endl;

#include "autodiff.h"


static const uint dim(2);

// for convenience, we declare these brief names:
typedef var<dim>              rvar;
typedef vector<double,dim>    rvector;
typedef vector<uint,dim>      ivector;
typedef vector<rvar,dim>      vect_var;
typedef matrix<double,dim>    rmatrix;


// define the system of equations here
vect_var f( const vect_var x )
{
 vect_var y;
 y(0) = atan2( x(1), x(0) ) - 1.0;
 y(1) =  x(0) + x(1) - 1.0;
 return y;
}



void output(const rmatrix j, const rvector x, const rvector y, const double res)
{
 cout.precision(12);
 cout<<"=================================="<<endl;
 cout<<"JAC:"<<endl;
 cout<<j(0,0)<<"  "<<j(0,1)<<endl;
 cout<<j(1,0)<<"  "<<j(1,1)<<endl;

 cout<<"-------------------------------------------------------------------------------"<<endl;
 cout<<"       x1            x2               y1                y2"<<endl;
 for(uint i=0; i<dim; i++)  cout<<x(i)<<"  ";
 for(uint i=0; i<dim; i++)  cout<<y(i)<<"  ";
 cout<<res<<endl;
 cout<<"-------------------------------------------------------------------------------"<<endl;
}

double uround;
void setmacheps(void) {

	double me = 1.0;
	while( 1.0+0.5*me > 1.0 ) {
		me=0.5*me;
	}
	uround=me*2.0;

}

int main(int argc, char** argv) {

	setmacheps();
	double epsilon=uround;

 //const double epsilon(1.0e-8); //(1.0e-7);
 const int maxdamp(25);
 rvector b,x0;
 rmatrix jac;
 vect_var x,y;
 ivector ip;
 double residual, old_res;
 bool singular;
 
 // set initial values:
 x0(0)=1.0;
 x0(1)=1.0;
 set_indep_var(x, x0 ); //initial guess

 y=f(x);
 old_res = l2norm(y);

 for(int iter=0; iter<55; iter++)
 {

  b=value(y); // value of vector function
  jac = jacobian(y);// analytic derivatives already here


  dec(jac,ip, singular); // LU factorization of jacobian
  if( singular ) 
  {
   cout<<"singular matrix in DEC"<<endl;
   exit(-1);
  }
  sol(jac,b,ip); // RHS=b is overwritten with solution


  set_indep_var(x, value(x)-b); //improve estimate of root
  y=f(x);
  residual=l2norm(y);


  // begin damping
  int ndamp=0;
  while( (residual>=old_res) && (residual>epsilon) )
  {
   ndamp++;
   if(ndamp>maxdamp) break;
   cout<<" damping...residual="<<residual<<endl;
   b = b/2.0;
   set_indep_var(x, value(x)+b);
   y=f(x);
   residual=l2norm(y);
  }; // end while (end damping)

  cout<<endl;
  cout<<"after damping, iter="<<iter<<endl;
  output( jac, value(x), value(y), residual );
  
  if(residual<epsilon) break; // good enough
  if((residual>old_res)&&(iter>10)) break; // diverging => as good as it gets

  old_res=residual;

 } // end for iter

 cout<<endl;
 cout<<"epsilon(double) = "<<epsilon<<endl;
 cout<<"...final values of (x,y) represent the root to the system"<<endl;
 cout<<" z1 = atan2( y, x ) - 1.0"<<endl;
 cout<<" z2 = x + y - 1.0"<<endl;
 cout<<endl;
 cout<<"representing that point on the line y = x + 1 that"<<endl;
 cout<<"is exactly 1 radian up from the positive x-axis"<<endl;
 cout<<endl;

 return 0;

}


