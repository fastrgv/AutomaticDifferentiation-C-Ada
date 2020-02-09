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

with mathtypes; use mathtypes;

package aderiv is

	pragma Elaborate_Body;

   type var is private;
   type indep_var is limited private;

	uround: real;
	procedure setmacheps;



------ interface routines -----------------------------------------

procedure set_indep_var( x : in out indep_var;  r : real );

function value(a:var) return real;
function value(a:indep_var) return real;

function deriv( a:var ) return real;


--  gradient operators

function "+"(a,b:var) return var ;
function "+"(a,b:indep_var) return var ;
function "+"(a:indep_var; b:var) return var ;
function "+"(a:var; b:indep_var) return var ;

function "+"(a:real; b:var) return var ;
function "+"(a:real; b:indep_var) return var ;

function "+"(a:integer; b:var) return var ;
function "+"(a:integer; b:indep_var) return var ;

function "+"(a:var; b:real) return var ;
function "+"(a:indep_var; b:real) return var ;

function "+"(a:var; b:integer) return var ;
function "+"(a:indep_var; b:integer) return var ;


function "-"(a,b:var) return var ;
function "-"(a,b:indep_var) return var ;
function "-"(a:indep_var; b:var) return var ;
function "-"(a:var; b:indep_var) return var ;

function "-"(a:var; b:real) return var ;
function "-"(a:indep_var; b:real) return var ;

function "-"(a:var; b:integer) return var ;
function "-"(a:indep_var; b:integer) return var ;

function "-"(a:real; b:var) return var ;
function "-"(a:real; b:indep_var) return var ;

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

function "*"(a:real; b:var) return var ;
function "*"(a:real; b:indep_var) return var ;

function "*"(a:var; b:real) return var ;
function "*"(a:indep_var; b:real) return var ;


function "/"(a,b:var) return var ;
function "/"(a,b:indep_var) return var;
function "/"(a:var; b:indep_var) return var;
function "/"(a:indep_var; b:var) return var;

function "/"(a:integer; b:var) return var ;
function "/"(a:integer; b:indep_var) return var ;

function "/"(a:var; b:integer) return var ;
function "/"(a:indep_var; b:integer) return var ;

function "/"(a:real; b:var) return var ;
function "/"(a:real; b:indep_var) return var ;

function "/"(a:var; b:real) return var ;
function "/"(a:indep_var; b:real) return var ;


function "**"(a,b:real) return real ;
function "**"(a:real; b:var) return var ;
function "**"(a:var; b:real) return var ;
function "**"(a,b:var) return var ;
function "**"(a:integer; b:var) return var;
function "**"(a:var; b:integer) return var;

function "**"(a:real; b:indep_var) return var ;
function "**"(a:indep_var; b:real) return var ;
function "**"(a,b:indep_var) return var ;
function "**"(a:integer; b:indep_var) return var;
function "**"(a:indep_var; b:integer) return var;

function "<"(a,b:var) return boolean;

--  standard gradient functions

function recip(a:var) return var;
function sqr(a:var) return var;

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

function tan(a:var) return var;
function asin(a:var) return var;

function atan(a:var) return var ;
function atan(a:indep_var) return var;

-------------- addendum ----------------------------------

--careful here...g represents known error in v
--for purposes if calculating differentials:
function set_var(v,g: real) return var;

--convert raw input data d into var;
--er represents error in its value,
--which we assume related to uround,
--in the absence of more specific knowledge:
function d2v( d: real ) return var;

function atan2(y,x:var) return var ;



private

   type      var is record
                      val:real;
                      grad:real;
                    end record;

   type indep_var is new var;

end aderiv;
