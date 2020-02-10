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
-- TBD

with interfaces.c;
with mathtypes, text_io, autodiff;
use text_io;

procedure new3eg is

	use mathtypes;
	use interfaces.c;

   package gradpak is new autodiff( dimension => 2 );
   package fio is new text_io.float_io(real);

   use  gradpak, fio;

   pi      : constant real   := 3.141592654;
   maxiter : constant integer := 55;
   epsilon : real;
   zero    : constant real   := 0.0;
   one     : constant real   := 1.0;

   b       :  rvector;
   jac     :  rmatrix;
   y       :  vect_var;
   x       :  indep_vect_var;
   ip      :  pervec;
   residual, old_residual : real;
   singular : boolean;


uround3: real;
--value used to control error estimates:
procedure setmacheps is
	me: real := 1.0;
begin
	loop
		exit when 1.0+0.5*me = 1.0;
		me:=0.5*me;
	end loop;
	uround3:=me*2.0;
end;




--
-- define the system of equations here
function f( x : indep_vect_var ) return vect_var is
  y : vect_var;
begin
   y(1) := atan2( x(2), x(1) ) - 1.0;
   y(2) := x(1) + x(2)  - 1.0;
   return y;
end f;


function norm( y : vect_var ) return real is
  max : real := 0.0;  r:rvector;
begin
  r := value(y);
  for i in index loop
    if( abs(r(i)) > max ) then
      max := abs(r(i));
    end if;
  end loop;
  return max;
end norm;


procedure output(c:rvector; res:real) is
begin
  put( "      x                     y                   residual");
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

	setmacheps;
	epsilon:=uround3;

	new_line;
	put_line("[real]Uround: "&real'image(uround3));
	put_line("using epsilon: "&real'image(epsilon));


  -- set initial values
  set_indep_var( x, (1.0,1.0) );

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
		put_line(" old b: " & real'image(b(1)) & " " & real'image(b(2)));
    sol( jac, b, ip );
		put_line(" new b: " & real'image(b(1)) & " " & real'image(b(2)));
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
-- y1 = atan2( x2, x1 ) - 0.1
-- y2 = x1 + x2 - 1.0
--

	new_line;
	put("...final values of (x,y) represent the root to the system");
	new_line;
	put(" z1 = atan2( y, x ) - 1.0");
	new_line;
	put(" z2 = x + y  - 1.0");
	new_line;
	new_line;


end new3eg;

