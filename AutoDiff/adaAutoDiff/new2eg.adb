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



-- demonstration of damped Newton's method for solving a 2-dimensional 
-- problem using AD:
--
-- find the root (x,y) of the function
--
-- z = exp( -x + y ) - 0.1
-- w = exp( -x - y ) - 0.1
--


with math_lib, text_io, autodiff;
use text_io;

procedure new2eg is

   package gradpak is new autodiff( dimension => 2 );
   package fio is new text_io.float_io(float);

   use  gradpak, fio, math_lib ;

   pi      : constant float   := 3.141592654;
   maxiter : constant integer := 55;
   epsilon : constant float   := 1.0e-7;
   zero    : constant float   := 0.0;
   one     : constant float   := 1.0;

   b       :  rvector;
   jac     :  rmatrix;
   y       :  vect_var;
   x       :  indep_vect_var;
   ip      :  pervec;
   residual, old_residual : float;
   singular : boolean;

-- define the system of equations here
function f( x : indep_vect_var ) return vect_var is
  y : vect_var;
begin
   y(1) := exp( - x(1) + x(2) ) - 0.1;
   y(2) := exp( - x(1) - x(2) ) - 0.1;
   return y;
end f;


function norm( y : vect_var ) return float is
  max : float := 0.0;  r:rvector;
begin
  r := value(y);
  for i in index loop
    if( abs(r(i)) > max ) then
      max := abs(r(i));
    end if;
  end loop;
  return max;
end norm;


procedure output(c:rvector; res:float) is
begin
  put( "      x                      y                    residual");
  new_line;
  for i in index loop
    put( c(i), aft => 13 );
    put("  ");
  end loop;
  put("     ");
  put( res, exp => 3 );
  new_line;
--  new_line;
end output;


-------- begin main program -----------
begin

  -- set initial values
  set_indep_var( x, (4.3,2.0) );

  y := f(x);
  old_residual := norm( y );
  output( value(x), old_residual );

  for iter in 1..maxiter loop
    b     := value( y );
    jac   := jacobian( y );
    dec( jac, ip, singular );
    if( singular ) then
      put("singularity in matrix decomposition");
      exit;
    end if;
		put_line(" old b: " & float'image(b(1)) & " " & float'image(b(2)));
    sol( jac, b, ip );
		put_line(" new b: " & float'image(b(1)) & " " & float'image(b(2)));
    set_indep_var( x, value(x) - b );
    y := f( x );
    residual := norm( y );
    while( (residual >= old_residual) and (residual > epsilon) ) loop
       put( " damping..." );
       new_line;
       b := b / 2.0;
       set_indep_var( x, value(x) + b );
       y := f(x);
       residual := norm( y );
    end loop;
    output( value(x), residual );
    exit when residual < epsilon;
    old_residual := residual;
  end loop;


--
-- find the root (x,y) of the function
--
-- z = exp( -x + y ) - 0.1
-- w = exp( -x - y ) - 0.1
--

	new_line;
	put("...final values of (x,y) represent the root to the system");
	new_line;
	put(" z = exp( -x + y ) - 0.1");
	new_line;
	put(" w = exp( -x - y ) - 0.1");
	new_line;
	new_line;


end new2eg;

