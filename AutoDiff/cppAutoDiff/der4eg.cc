
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


// der4eg.cc

#include <math.h>
#include <iostream>
using std::cout;
using std::endl;

#include "autodiff.h"


static const uint dim(3);

// for convenience, we declare these brief names:
typedef var<dim>              rvar;
typedef vector<double,dim>    rvector;
typedef vector<rvar,dim>      vect_var;


// define the system of equations here
rvar f( const vect_var x )
{
 rvar y = 
 	1.0 + x(0) + x(1) + x(2)
		+ x(0)*x(1) + x(1)*x(2) + x(0)*x(2)
		+ x(0)*x(1)*x(2) + exp( x(0)/x(1) + x(1)/x(2) );

 return y;
}


int main(int argc, char** argv) {


 rvector x0;
 vect_var x;
 rvar y;
 
 // set initial values:
 x0(0)=1.0;
 x0(1)=2.0;
 x0(2)=3.0;
 set_indep_var(x, x0 ); //initial guess

 y=f(x);


 cout<<endl;
 cout<<"    u = "<< value(y) <<endl;

 cout<<"du/dx = "<< deriv(y,0) <<endl;
 cout<<"du/dy = "<< deriv(y,1) <<endl;
 cout<<"du/dz = "<< deriv(y,2) <<endl;

 cout<<endl;

 return 0;

}


