--
-- Copyright (C) 2020  <fastrgv@gmail.com>
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


-- demonstration of damped Newton's method for solving a 1-dimensional 
-- problem usind AD:
--
-- find the root of f(x) = x - cos x

with interfaces.c;
with text_io, mathtypes, aderiv;
use  text_io;

procedure new1eg is

	use mathtypes;
	use interfaces.c;

 package fio is new text_io.float_io(real);
 use  aderiv, fio ;

 pi      : constant real   := 3.141592654;
 maxiter : constant integer := 55;
 epsilon : real; --   := 1.0e-12;
 zero    : constant real   := 0.0;
 one     : constant real   := 1.0;

 y       :  var;
 x       :  indep_var;
 b, jac, residual, old_residual : real;

function f( x : indep_var ) return var is
  y : var;
begin
   y := cos(x) - x;
   return y;
end f;

begin

	epsilon:=uround;
	new_line;
	put_line("[real]Uround: "&real'image(uround));

-- set the initial value
  set_indep_var( x, 0.5 );

  y := f(x);
  old_residual := abs( value(y) );
  put("        x              y"); new_line;
  put(value(x),aft=>13); put("  "); put(value(y),aft=>13); new_line;

  for iter in 1..maxiter loop
    b     := value( y );
    jac   := deriv( y );
    b := b/jac;
    set_indep_var( x, value(x) - b );
    y := f( x );
    residual := abs( value(y) );
    while( (residual >= old_residual) and (residual > epsilon) ) loop
       put( " damping..." );  new_line;
       b := b / 2.0;
       set_indep_var( x, value(x) + b );
       y := f(x);
       residual := abs( value(y) );
    end loop;
    put(value(x),aft=>13); put("  ");put(value(y),aft=>13); new_line;
    exit when residual < epsilon;
    old_residual := residual;
  end loop;

  new_line;
  put("...final value of x is root of equation y = x-cos(x)");
  new_line;
  new_line;

end new1eg;
