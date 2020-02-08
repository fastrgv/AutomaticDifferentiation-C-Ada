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

with math_lib;   use math_lib;

generic
   dimension : in positive;

package autodiff is

   type var is private;
   type indep_var is limited private;

   subtype index is integer range 1..dimension;
   type rvector is array( index ) of float;
   type rmatrix is array(index,index) of float;
   type pervec  is array(index) of index;

   type vect_var is array(index) of var;
   type indep_vect_var is array(index) of indep_var;



------ linear solver routines -----------

procedure dec(a:in out rmatrix; ip:out pervec; singular:out boolean) ;

procedure sol(a:rmatrix; b:in out rvector; ip:pervec) ;



------ interface routines -----------------------------------------

function value(a:var) return float;
function value(a:indep_var) return float;

function deriv( a:var; i:index ) return float;

procedure set_indep_var(x : in out indep_vect_var; r:rvector);

function value(x:vect_var) return rvector;

function value(x:indep_vect_var) return rvector;

function gradient(a:var) return rvector;

function jacobian(x:vect_var) return rmatrix;


------- auxillary routines ---------------------------------

function "+"(u,v:rvector) return rvector;

function "-"(u,v:rvector) return rvector;

function "*"(u:float; v:rvector) return rvector;

function "*"(v:rvector; u:float) return rvector;

function "/"(v:rvector; u:float) return rvector;


--  gradient operators

function "+"(a,b:var) return var ;
function "+"(a,b:indep_var) return var ;
function "+"(a:indep_var; b:var) return var ;
function "+"(a:var; b:indep_var) return var ;

function "+"(a:float; b:var) return var ;
function "+"(a:float; b:indep_var) return var ;

function "+"(a:integer; b:var) return var ;
function "+"(a:integer; b:indep_var) return var ;

function "+"(a:var; b:float) return var ;
function "+"(a:indep_var; b:float) return var ;

function "+"(a:var; b:integer) return var ;
function "+"(a:indep_var; b:integer) return var ;


function "-"(a,b:var) return var ;
function "-"(a,b:indep_var) return var ;
function "-"(a:indep_var; b:var) return var ;
function "-"(a:var; b:indep_var) return var ;

function "-"(a:var; b:float) return var ;
function "-"(a:indep_var; b:float) return var ;

function "-"(a:var; b:integer) return var ;
function "-"(a:indep_var; b:integer) return var ;

function "-"(a:float; b:var) return var ;
function "-"(a:float; b:indep_var) return var ;

function "-"(a:integer; b:var) return var ;
function "-"(a:integer; b:indep_var) return var ;

-- unary minus
function "-"(a:var) return var;
function "-"(a:indep_var) return var;


function "*"(a,b:var) return var ;
function "*"(a,b:indep_var) return var ;
function "*"(a:indep_var; b:var) return var ;
function "*"(a:var; b:indep_var) return var ;

function "*"(a:integer; b:var) return var ;
function "*"(a:integer; b:indep_var) return var ;

function "*"(a:var; b:integer) return var ;
function "*"(a:indep_var; b:integer) return var ;

function "*"(a:float; b:var) return var ;
function "*"(a:float; b:indep_var) return var ;

function "*"(a:var; b:float) return var ;
function "*"(a:indep_var; b:float) return var ;


function "/"(a,b:var) return var ;
function "/"(a,b:indep_var) return var;
function "/"(a:var; b:indep_var) return var;
function "/"(a:indep_var; b:var) return var;

function "/"(a:integer; b:var) return var ;
function "/"(a:integer; b:indep_var) return var ;

function "/"(a:var; b:integer) return var ;
function "/"(a:indep_var; b:integer) return var ;

function "/"(a:float; b:var) return var ;
function "/"(a:float; b:indep_var) return var ;

function "/"(a:var; b:float) return var ;
function "/"(a:indep_var; b:float) return var ;


function "**"(a,b:float) return float ;
function "**"(a:float; b:var) return var ;
function "**"(a:var; b:float) return var ;
function "**"(a,b:var) return var ;
function "**"(a:integer; b:var) return var;
function "**"(a:var; b:integer) return var;

function "**"(a:float; b:indep_var) return var ;
function "**"(a:indep_var; b:float) return var ;
function "**"(a,b:indep_var) return var ;
function "**"(a:integer; b:indep_var) return var;
function "**"(a:indep_var; b:integer) return var;



--  standard gradient functions

function sqrt(a:var) return var ;
function sqrt(a:indep_var) return var;

function exp(a:var) return var ;
function exp(a:indep_var) return var;

function ln(a:var) return var ;
function ln(a:indep_var) return var;

function sin(a:var) return var ;
function sin(a:indep_var) return var;

function cos(a:var) return var ;
function cos(a:indep_var) return var;

function atan(a:var) return var ;
function atan(a:indep_var) return var;


private

   type      var is record
                      val:float;
                      grad:rvector;
                    end record;

   type indep_var is new var;

end autodiff;
