


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




// vmv.h
//
// Vector-Matrix-Var template classes
// declarations and definitions
//

#ifndef VMV_H
#define VMV_H

#include <cassert>

typedef unsigned int uint;

template <typename T, uint dim>
class vector {
public:
 T * v;


//default constructor:
 vector<T,dim>() 
 {  
  v = new T[dim];
 }

// copy constructor:
 vector<T,dim>( const vector<T,dim>& rhs ) 
 {  
  v = new T[dim];
  for(uint i=0; i<dim; i++) v[i]=rhs(i);
 }


//destructor:
 ~vector(){ delete [] v; }
 
 
 // index operators:
 T & operator()(const uint i) 
 {  
  assert(i<dim);
  assert(i>=0);
  return v[i]; 
 }
 T   operator()(const uint i) const 
 {  
  assert(i<dim);
  assert(i>=0);
  return v[i]; 
 }
 

 // assignment
 void operator=( const vector<T,dim>& rhs )
 {
  for(uint i=0; i<dim; i++) v[i]=rhs(i);
 }

 // mathematical operator functions:

 vector<T,dim> plus(const vector<T,dim>& rhs) const
 {
  vector<T,dim> aux(*this);
  for(uint i=0; i<dim; i++) aux.v[i]+=rhs(i);
  return aux;
 }

 vector<T,dim> unaryminus() const
 {
  vector<T,dim> aux(*this);
  for(uint i=0; i<dim; i++) aux.v[i] *= -1.0;
  return aux;
 }

 vector<T,dim> minus(const vector<T,dim>& rhs) const
 {
  vector<T,dim> aux(*this);
  for(uint i=0; i<dim; i++) aux.v[i]-=rhs(i);
  return aux;
 }

 vector<T,dim> plus(const T rhs) const
 {
  vector<T,dim> aux(*this);
  for(uint i=0; i<dim; i++) aux.v[i]+=rhs;
  return aux;
 }

 vector<T,dim> minus(const T rhs) const
 {
  vector<T,dim> aux(*this);
  for(uint i=0; i<dim; i++) aux.v[i]-=rhs;
  return aux;
 }

 vector<T,dim>  times(const T rhs) const
 {
  vector<T,dim> aux(*this);
  for(uint i=0; i<dim; i++) aux.v[i]*=rhs;
  return aux;
 }

 vector<T,dim>  divide(const T rhs) const
 {
  vector<T,dim> aux(*this);
  for(uint i=0; i<dim; i++) aux.v[i]/=rhs;
  return aux;
 }

}; // end class vector



// define a square matrix:
template <typename T, uint dim>
class matrix {
private:
 vector<T,dim> * m;

public:

 // default constructor
 matrix<T,dim>()
 {
  m = new vector<T,dim> [dim];
 }

 // copy constructor
 matrix<T,dim>( const matrix<T,dim>& rhs )
 {
  m = new vector<T,dim> [dim];
  for(uint col=0; col<dim; col++) m[col]=rhs(col);
 }

 ~matrix(){ delete [] m;}


 // single index operators:
 vector<T,dim> & operator()(const uint col)
 {
  return m[col];
 }
 vector<T,dim>  operator()(const uint col) const
 {
  return m[col];
 }


 // double index operators:
 T & operator()(const uint row, const uint col)
 {
  return m[col](row);
 }
 T  operator()(const uint row, const uint col) const
 {
  return m[col](row);
 }


 // assignment
 void operator=( const matrix<T,dim>& rhs )
 {
  for(uint col=0; col<dim; col++) m[col]=rhs(col);
 }


}; // end class matrix


 
template <uint dim>
class var {
public:
 double val;
 vector<double,dim> grad;

 // default constructor:
 var<dim>()
 {
  val=0.0;
  for(uint i=0; i<dim; i++) grad(i)=0.0;
 };

 // copy constructor:
 var<dim>( const var<dim>& rhs ) 
   { val=rhs.val; grad=rhs.grad; };

 // unary minus
 var<dim> operator-() const
 {
  var<dim> aux(*this);
  aux.val *= -1.0;
  for(uint i=0; i<dim; i++) aux.grad(i) *= -1.0;
  return aux;
 }

 // assignment operator:
 void operator=(const var<dim>& rhs)
 {
  val = rhs.val;
  grad = rhs.grad;
 }

}; // end class var

 
#endif


