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

with mathtypes; 
with ada.numerics;

with interfaces.c; 

package body autodiff is

use mathtypes;
use mathtypes.math_lib;
use interfaces.c;

   zero : constant real := 0.0;
   one  : constant real := 1.0;
   two  : constant real := 2.0;
	onepi: constant real := ada.numerics.pi;
	halfpi: constant real := 0.5*onepi;

-- begin linear solver routines

procedure dec(a:in out rmatrix; ip:out pervec; singular:out boolean) is
 i1,k1, m, n1, n, rows, cols : index;
 t : real;
begin
 singular := false;
 n := index'last;
 i1 := index'first;
 rows := n;
 cols := n;
 ip(n):=i1;
 n1:= index'pred(n) ;
 for k in i1..n1 loop
  k1 := index'succ(k);
  m:=k;
  for i in k1..n loop
   if (abs(a(i,k))>abs(a(m,k))) then
     m:=i;
   end if;
  end loop;
  ip(k):=m;
  t:=a(m,k);
  if (m /= k) then
   a(m,k):=a(k,k);
   a(k,k):=t;
  end if;
  if (t=zero) then
   singular := true;
   return;
  end if;
  t:=one/t;
  for i in k1..n loop
   a(i,k) := -a(i,k)*t;
  end loop;
  for j in k1..n loop
   t:=a(m,j);
   a(m,j):=a(k,j);
   a(k,j):=t;
   if (t /= zero) then
     for i in k1..n loop
      a(i,j):=a(i,j)+a(i,k)*t;
     end loop;
   end if;
  end loop;
 end loop;
end dec;

procedure sol(a:rmatrix; b:in out rvector; ip:pervec) is
 i1, n, j, k1, m, n1, rows, cols : index;
 t : real;
begin
  n := index'last;
  i1 := index'first;
  rows := n;
  cols := n;
  n1 := index'pred(n);
  for k in i1..n1 loop
   k1:= index'succ(k);
   m:=ip(k);
   t:=b(m);
   b(m):=b(k);
   b(k):=t;
   for i in k1..n loop
    b(i):=b(i)+a(i,k)*t;
   end loop;
  end loop;
  k1 := n;
  for k in i1..n1 loop
   j := k1;
   k1 := index'pred(k1);
   b(j):=b(j)/a(j,j);
   t:=-b(j);
   for i in i1..k1 loop
    b(i):=b(i)+a(i,j)*t;
   end loop;
  end loop;
  b(i1):=b(i1)/a(i1,i1);
end sol;

-- end linear solver routines


-- begin interface routines

function value(a:var) return real is
begin
  return a.val;
end value;

function value(a:indep_var) return real is
begin
  return a.val;
end value;

function deriv(a:var; i:index) return real is
begin
  return a.grad(i);
end deriv;


procedure set_indep_var( x : in out indep_vect_var;  r : rvector ) is
begin
  for i in index loop
    x(i).val := r(i);
    for j in index loop
      x(i).grad(j) := zero;
    end loop;
    x(i).grad(i) := one;
  end loop;
end set_indep_var;


--procedure set_var( x : in out var;  r: real; g: rvector ) is
--begin
--	x.val:=r;
--	x.grad:=g;
--end set_var;




function value(x:vect_var) return rvector is
  b:rvector;
begin
  for i in index loop
    b(i) := x(i).val;
  end loop;
  return b;
end value;

function value(x:indep_vect_var) return rvector is
  b:rvector;
begin
  for i in index loop
    b(i) := x(i).val;
  end loop;
  return b;
end value;


function gradient(a:var) return rvector is
begin
  return a.grad;
end gradient;

function jacobian(x:vect_var) return rmatrix is
  jac:rmatrix;
begin
  for i in index loop
    for j in index loop
      jac(i,j) := x(i).grad(j);
    end loop;
  end loop;
  return jac;
end jacobian;

function "+"(u,v:rvector) return rvector is
  w:rvector;
begin
  for i in index loop
    w(i) := u(i) + v(i);
  end loop;
  return w;
end "+";

function "-"(u,v:rvector) return rvector is
  w:rvector;
begin
  for i in index loop
    w(i) := u(i) - v(i);
  end loop;
  return w;
end "-";

function "*"(u:real; v:rvector) return rvector is
  w:rvector;
begin
  for i in index loop
    w(i) := u * v(i);
  end loop;
  return w;
end "*";

function "*"(v:rvector; u:real) return rvector is
  w:rvector;
begin
  for i in index loop
    w(i) := u * v(i);
  end loop;
  return w;
end "*";

function "/"(v:rvector; u:real) return rvector is
  w:rvector;
begin
  for i in index loop
    w(i) := v(i) / u;
  end loop;
  return w;
end "/";

-- end interface routines

-- begin conversion routines

function makvar( r:real ) return var is
 y:var;
begin
  y.val := r;
  for i in index loop
    y.grad(i) := zero;
  end loop;
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
 for i in index loop
  y.grad(i):=a.grad(i)+b.grad(i);
 end loop;
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
 for i in index loop
  y.grad(i):=a.grad(i)-b.grad(i);
 end loop;
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
 for i in index loop
  y.grad(i):=a.val*b.grad(i)+a.grad(i)*b.val;
 end loop;
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
 for i in index loop
  if a.grad(i)=zero
  then y.grad(i):=zero;
  else y.grad(i):=a.grad(i)/b.val;
  end if;
  if b.grad(i) /= zero
     then y.grad(i):=y.grad(i)-b.grad(i)*y.val/b.val;
  end if;
 end loop;
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
  for i in index loop
    res.grad(i) := zero;
  end loop;
  if(a>zero) then
    m:=res.val*log(a);
    for i in index loop
       if(b.grad(i) /= zero) then
          res.grad(i) := m*b.grad(i);
       end if;
    end loop;
  elsif(a<zero) then
    for i in index loop
      if(b.grad(i) /= zero) then
-- attempting to raise a negative number to a nonconstant power
         raise numeric_error;
      end if;
    end loop;
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
    for i in index loop
      res.grad(i) := zero;
    end loop;
  else
    m := b * res.val / a.val;
    for i in index loop
      if(a.grad(i)=zero) then
        res.grad(i) := zero;
      else
        res.grad(i) := m * a.grad(i);
      end if;
    end loop;
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
  for i in index loop
    res.grad(i) := res.grad(i) + v.grad(i);
  end loop;
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


-- end gradient operators


-- begin standard gradient functions

function sqrt(a:var) return var is
 y:var; m:real;
begin
 y.val:=sqrt(a.val);
 m:=two*y.val;
 for i in index loop
  if a.grad(i)=zero then y.grad(i):=zero;
  else y.grad(i):=a.grad(i)/m;
  end if;
 end loop;
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
 for i in index loop
  if a.grad(i)=zero then y.grad(i):=zero;
  else y.grad(i):=y.val*a.grad(i);
  end if;
 end loop;
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
 for i in index loop
  if a.grad(i)=zero then y.grad(i):=zero;
  else y.grad(i):=a.grad(i)/a.val;
  end if;
 end loop;
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
 for i in index loop
  if a.grad(i)=zero then y.grad(i):=zero;
  else y.grad(i):=m*a.grad(i);
  end if;
 end loop;
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
 for i in index loop
  if a.grad(i)=zero then y.grad(i):=zero;
  else y.grad(i):=m*a.grad(i);
  end if;
 end loop;
 return y;
end cos;
function cos(a:indep_var) return var is
begin
  return cos(var(a));
end cos;





function atan(a:var) return var is
 y:var; m:real;
begin
 y.val := arctan( a.val );
 m := one + a.val * a.val;
 for i in index loop
  y.grad(i):=a.grad(i)/m;
 end loop;
 return y;
end atan;
function atan(a:indep_var) return var is
begin
  return atan(var(a));
end atan;





-- end standard gradient functions


------ addendum begin ---------------------------------

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


begin -- autodiff package body [initialization]

	setmacheps; -- prepares uround

end autodiff;

