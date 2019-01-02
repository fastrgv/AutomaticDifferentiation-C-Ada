--
-- Copyright (C) 2014  <fastrgv@gmail.com>
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

with math_lib;  use math_lib;

-------------------------- body --------------------------

package body aderiv is
   zero : constant float := 0.0;
   one  : constant float := 1.0;
   two  : constant float := 2.0;


-- begin interface routines

procedure set_indep_var( x : in out indep_var;  r : float ) is
begin
  x.val := r;
  x.grad := one;
end set_indep_var;

function value(a:var) return float is
begin
  return a.val;
end value;

function value(a:indep_var) return float is
begin
  return a.val;
end value;

function deriv( a : var ) return float is
begin
  return a.grad;
end deriv;

-- end interface routines


-- begin conversion routines

function makvar( r:float ) return var is
 y:var;
begin
  y.val := r;
  y.grad := zero;
  return y;
end;

function makvar( i:integer ) return var is
begin
  return makvar( float(i) );
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

function "+"(a:float; b:var) return var is
begin
 return makvar(a)+b;
end "+";

function "+"(a:float; b:indep_var) return var is
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

function "+"(a:var; b:float) return var is
begin
 return a+makvar(b);
end "+";
function "+"(a:indep_var; b:float) return var is
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

function "-"(a:var; b:float) return var is
begin
 return a-makvar(b);
end "-";
function "-"(a:indep_var; b:float) return var is
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

function "-"(a:float; b:var) return var is
begin
 return makvar(a)-b;
end "-";
function "-"(a:float; b:indep_var) return var is
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

function "*"(a:float; b:var) return var is
begin
 return makvar(a)*b;
end "*";
function "*"(a:float; b:indep_var) return var is
begin
 return makvar(a)*var(b);
end "*";

function "*"(a:var; b:float) return var is
begin
 return a*makvar(b);
end "*";
function "*"(a:indep_var; b:float) return var is
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

function "/"(a:float; b:var) return var is
begin
 return makvar(a)/b;
end "/";
function "/"(a:float; b:indep_var) return var is
begin
 return makvar(a)/var(b);
end "/";

function "/"(a:var; b:float) return var is
begin
 return a/makvar(b);
end "/";
function "/"(a:indep_var; b:float) return var is
begin
 return var(a)/makvar(b);
end "/";

function "**"(a,b:float) return float is
k:integer; u,v:float;
begin
  if(a=one) then
     u:=one;
  elsif(a=zero and b>zero) then
     u:=zero;
  else
     k:=integer(b);
     if(float(k)>b) then
        k:=k-1;
     end if;
     u:= a**k;
     v:=b-float(k);
     if(v /= zero) then
        u:=u*exp(v*log(a));
     end if;
  end if;
  return u;
end "**";

function "**"(a:float; b:var) return var is
m:float; res:var;
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
function "**"(a:float; b:indep_var) return var is
begin
  return a**var(b);
end "**";

function "**"(a:var; b:float) return var is
m:float; res:var;
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
function "**"(a:indep_var; b:float) return var is
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
  return float(a)**b;
end "**";
function "**"(a:integer; b:indep_var) return var is
begin
  return float(a)**var(b);
end "**";

function "**"(a:var; b:integer) return var is
begin
  return a**float(b);
end;
function "**"(a:indep_var; b:integer) return var is
begin
  return var(a)**float(b);
end;


-- end gradient operators


-- begin standard gradient functions

function sqrt(a:var) return var is
 y:var; m:float;
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
 y:var; m:float;
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
 y:var; m:float;
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

function atan(a:var) return var is
 y:var; m:float;
begin
 y.val := arctan( a.val );
 m := one + a.val * a.val;
  if a.grad=zero then y.grad:=zero;
  else y.grad:=a.grad/m;
  end if;
 return y;
end atan;
function atan(a:indep_var) return var is
begin
  return atan(var(a));
end atan;

-- end standard gradient functions


end aderiv;

