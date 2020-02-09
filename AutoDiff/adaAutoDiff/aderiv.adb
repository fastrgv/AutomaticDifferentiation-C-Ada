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

with ada.numerics;
with mathtypes;  use mathtypes;
with interfaces.c; use interfaces.c;

-------------------------- body --------------------------

package body aderiv is

	use mathtypes.math_lib;
	use interfaces.c;

   zero : constant real := 0.0;
   one  : constant real := 1.0;
   two  : constant real := 2.0;
   onepi  : constant real := ada.numerics.pi;
   halfpi  : constant real := 0.5*ada.numerics.pi;



-- begin interface routines

procedure set_indep_var( x : in out indep_var;  r : real ) is
begin
  x.val := r;
  x.grad := one;
end set_indep_var;




function value(a:var) return real is
begin
  return a.val;
end value;

function value(a:indep_var) return real is
begin
  return a.val;
end value;

function deriv( a : var ) return real is
begin
  return a.grad;
end deriv;

-- end interface routines


-- begin conversion routines

function makvar( r:real ) return var is
 y:var;
begin
  y.val := r;
  y.grad := zero;
  return y;
end;

function makvar( i:integer ) return var is
begin
  return makvar( real(i) );
end;

-- end conversion routines


-- begin gradient operators

function "+"(a,b:var) return var is
 y:var;
begin
 y.val := a.val + b.val;
  y.grad:=a.grad+b.grad;
 return y;
end "+";

function "+"(a,b:indep_var) return var is
begin
 return var(a)+var(b);
end "+";

function "+"(a:indep_var; b:var) return var is
begin
 return var(a)+b;
end "+";

function "+"(a:var; b:indep_var) return var is
begin
 return a+var(b);
end "+";

function "+"(a:real; b:var) return var is
begin
 return makvar(a)+b;
end "+";

function "+"(a:real; b:indep_var) return var is
begin
 return makvar(a)+var(b);
end "+";

function "+"(a:integer; b:var) return var is
begin
 return makvar(a)+b;
end "+";
function "+"(a:integer; b:indep_var) return var is
begin
 return makvar(a)+var(b);
end "+";

function "+"(a:var; b:real) return var is
begin
 return a+makvar(b);
end "+";
function "+"(a:indep_var; b:real) return var is
begin
 return var(a)+makvar(b);
end "+";

function "+"(a:var; b:integer) return var is
begin
 return a+makvar(b);
end "+";
function "+"(a:indep_var; b:integer) return var is
begin
 return var(a)+makvar(b);
end "+";



function "-"(a,b:var) return var is
 y:var;
begin
 y.val:=a.val-b.val;
  y.grad:=a.grad-b.grad;
 return y;
end "-";
function "-"(a,b:indep_var) return var is
begin
 return var(a)-var(b);
end "-";
function "-"(a:indep_var; b:var) return var is
begin
 return var(a)-b;
end "-";
function "-"(a:var; b:indep_var) return var is
begin
 return a-var(b);
end "-";

function "-"(a:var; b:real) return var is
begin
 return a-makvar(b);
end "-";
function "-"(a:indep_var; b:real) return var is
begin
 return var(a)-makvar(b);
end "-";

function "-"(a:var; b:integer) return var is
begin
 return a-makvar(b);
end "-";
function "-"(a:indep_var; b:integer) return var is
begin
 return var(a)-makvar(b);
end "-";

function "-"(a:real; b:var) return var is
begin
 return makvar(a)-b;
end "-";
function "-"(a:real; b:indep_var) return var is
begin
 return makvar(a)-var(b);
end "-";

function "-"(a:integer; b:var) return var is
begin
 return makvar(a)-b;
end "-";
function "-"(a:integer; b:indep_var) return var is
begin
 return makvar(a)-var(b);
end "-";

-- unary minus
function "-"(a:var) return var is
begin
  return 0.0 - a;
end;
function "-"(a:indep_var) return var is
begin
  return 0.0 - var(a);
end;



function "*"(a,b:var) return var is
 y:var;
begin
 y.val:=a.val*b.val;
  y.grad:=a.val*b.grad+a.grad*b.val;
 return y;
end "*";

function "*"(a,b:indep_var) return var is
begin
  return var(a)*var(b);
end "*";

function "*"(a:indep_var; b:var) return var is
begin
  return var(a)*b;
end "*";

function "*"(a:var; b:indep_var) return var is
begin
  return a*var(b);
end "*";

function "*"(a:integer; b:var) return var is
begin
 return makvar(a)*b;
end "*";
function "*"(a:integer; b:indep_var) return var is
begin
 return makvar(a)*var(b);
end "*";

function "*"(a:var; b:integer) return var is
begin
 return a*makvar(b);
end "*";
function "*"(a:indep_var; b:integer) return var is
begin
 return var(a)*makvar(b);
end "*";

function "*"(a:real; b:var) return var is
begin
 return makvar(a)*b;
end "*";
function "*"(a:real; b:indep_var) return var is
begin
 return makvar(a)*var(b);
end "*";

function "*"(a:var; b:real) return var is
begin
 return a*makvar(b);
end "*";
function "*"(a:indep_var; b:real) return var is
begin
 return var(a)*makvar(b);
end "*";


function "/"(a,b:var) return var is
 y:var;
begin
 y.val:=a.val/b.val;
  if a.grad=zero
  then y.grad:=zero;
  else y.grad:=a.grad/b.val;
  end if;
  if b.grad /= zero
     then y.grad:=y.grad-b.grad*y.val/b.val;
  end if;
 return y;
end "/";

function "/"(a,b:indep_var) return var is
begin
  return var(a)/var(b);
end "/";
function "/"(a:indep_var; b:var) return var is
begin
  return var(a)/b;
end "/";
function "/"(a:var; b:indep_var) return var is
begin
  return a/var(b);
end "/";

function "/"(a:integer; b:var) return var is
begin
 return makvar(a)/b;
end "/";
function "/"(a:integer; b:indep_var) return var is
begin
 return makvar(a)/var(b);
end "/";

function "/"(a:var; b:integer) return var is
begin
 return a/makvar(b);
end "/";
function "/"(a:indep_var; b:integer) return var is
begin
 return var(a)/makvar(b);
end "/";

function "/"(a:real; b:var) return var is
begin
 return makvar(a)/b;
end "/";
function "/"(a:real; b:indep_var) return var is
begin
 return makvar(a)/var(b);
end "/";

function "/"(a:var; b:real) return var is
begin
 return a/makvar(b);
end "/";
function "/"(a:indep_var; b:real) return var is
begin
 return var(a)/makvar(b);
end "/";

function "**"(a,b:real) return real is
k:integer; u,v:real;
begin
  if(a=one) then
     u:=one;
  elsif(a=zero and b>zero) then
     u:=zero;
  else
     k:=integer(b);
     if(real(k)>b) then
        k:=k-1;
     end if;
     u:= a**k;
     v:=b-real(k);
     if(v /= zero) then
        u:=u*exp(v*log(a));
     end if;
  end if;
  return u;
end "**";

function "**"(a:real; b:var) return var is
m:real; res:var;
begin
  res.val := a**b.val;
    res.grad := zero;
  if(a>zero) then
    m:=res.val*log(a);
       if(b.grad /= zero) then
          res.grad := m*b.grad;
       end if;
  elsif(a<zero) then
      if(b.grad /= zero) then
-- attempting to raise a negative number to a nonconstant power
         raise numeric_error;
      end if;
  end if;
  return res;
end "**";
function "**"(a:real; b:indep_var) return var is
begin
  return a**var(b);
end "**";

function "**"(a:var; b:real) return var is
m:real; res:var;
begin
  res.val := a.val**b;
  if(b=one) then
    res.grad := a.grad;
  elsif( b=zero or ( a.val=zero and b>zero ) ) then
      res.grad := zero;
  else
    m := b * res.val / a.val;
      if(a.grad=zero) then
        res.grad := zero;
      else
        res.grad := m * a.grad;
      end if;
  end if;
  return res;
end "**";
function "**"(a:indep_var; b:real) return var is
begin
  return var(a)**b;
end "**";


function "**"(a,b:var) return var is
v,res:var;
begin
  res := a.val ** b;
  v   := a ** b.val;
    res.grad := res.grad + v.grad;
  return res;
end "**";
function "**"(a,b:indep_var) return var is
begin
  return var(a)**var(b);
end "**";


function "**"(a:integer; b:var) return var is
begin
  return real(a)**b;
end "**";
function "**"(a:integer; b:indep_var) return var is
begin
  return real(a)**var(b);
end "**";

function "**"(a:var; b:integer) return var is
begin
  return a**real(b);
end;
function "**"(a:indep_var; b:integer) return var is
begin
  return var(a)**real(b);
end;

function "<"(a,b: var) return boolean is
begin
	return a.val<b.val;
end;


-- end gradient operators


-- begin standard gradient functions

function recip(a:var) return var is
begin
	return 1.0/a;
end;

function sqr(a:var) return var is
begin
	return a*a;
end;

function sqrt(a:var) return var is
 y:var; m:real;
begin
 y.val:=sqrt(a.val);
 m:=two*y.val;
  if a.grad=zero then y.grad:=zero;
  else y.grad:=a.grad/m;
  end if;
 return y;
end sqrt;
function sqrt(a:indep_var) return var is
begin
  return sqrt(var(a));
end sqrt;

function exp(a:var) return var is
 y:var;
begin
 y.val:=exp(a.val);
  if a.grad=zero then y.grad:=zero;
  else y.grad:=y.val*a.grad;
  end if;
 return y;
end exp;
function exp(a:indep_var) return var is
begin
  return exp(var(a));
end exp;

function ln(a:var) return var is
 y:var;
begin
 y.val:=log(a.val);
  if a.grad=zero then y.grad:=zero;
  else y.grad:=a.grad/a.val;
  end if;
 return y;
end ln;
function ln(a:indep_var) return var is
begin
  return ln(var(a));
end ln;

function sin(a:var) return var is
 y:var; m:real;
begin
 y.val:=sin(a.val);
 m:=cos(a.val);
  if a.grad=zero then y.grad:=zero;
  else y.grad:=m*a.grad;
  end if;
 return y;
end sin;
function sin(a:indep_var) return var is
begin
  return sin(var(a));
end sin;

function cos(a:var) return var is
 y:var; m:real;
begin
 y.val:=cos(a.val);
 m:=-sin(a.val);
  if a.grad=zero then y.grad:=zero;
  else y.grad:=m*a.grad;
  end if;
 return y;
end cos;
function cos(a:indep_var) return var is
begin
  return cos(var(a));
end cos;


function tan(a:var) return var is
begin
	return sin(a)/cos(a);
end;

function asin(a:var) return var is
	a2: var := a*a;
begin
	return atan( a / sqrt(one-a2) );
end;


function atan(a:var) return var is
	y:var; 
	m: constant real := one + a.val*a.val;
begin
	y.val := arctan( a.val );
	y.grad:=a.grad/m;
	return y;
end atan;

function atan(a:indep_var) return var is
begin
  return atan(var(a));
end atan;

-- end standard gradient functions

-------------differential addendum-----------------------

--careful here...g represents known error in v
--...for purposes of tracking differentials.
--If g is unknown, use d2v().
function set_var(v, g : real) return var is
	z: var;
begin
	z.val := v;
	z.grad:= g;
	return z;
end;

--convert raw input data d into var;
--er represents error in its value,
--which we assume related to uround,
--in the absence of more specific knowledge:
function d2v( d: real ) return var is
	er: real;
begin
	if abs(d)>1.0 then
		er:=abs(d)*uround;
	else
		er:=uround;
	end if;
	return set_var(d,er);
end d2v;

function atan2(y,x: var) return var is
	z: var;
	m: constant real := y.val*y.val+x.val*x.val;
begin

	z.val := arctan(y.val,x.val); --intrinsic

	--z.grad := (x.val*y.grad-y.val*x.grad)/m;
	--quotient + chain rule

	--for the purposes of using differentials
	--to estimate errors I think I should do:
	z.grad := 
	(
		abs(x.val*y.grad) + abs(y.val*x.grad)
	) / m;

	return z;
end;



--value used to control error estimates:
procedure setmacheps is
	me: real := 1.0;
begin
	loop
		exit when 1.0+0.5*me = 1.0;
		me:=0.5*me;
	end loop;
	uround:=me*2.0;
end;


begin -- aderiv package body (initialization)

	setmacheps; -- prepares uround

end aderiv;

