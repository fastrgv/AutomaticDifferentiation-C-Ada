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


package aderiv is

	pragma Elaborate_Body;

   type var is private;
   type indep_var is limited private;

	uround: float;
	procedure setmacheps;


	-- carefully generate a var from a float:
	function d2v( d: float ) return var;


------ interface routines -----------------------------------------

--careful here...g represents error in v:
function set_var(v,g: float) return var;

procedure set_indep_var( x : in out indep_var;  r : float );

function value(a:var) return float;
function value(a:indep_var) return float;

function deriv( a:var ) return float;


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

function atan2(dy,dx:var) return var ;

private

   type      var is record
                      val:float;
                      grad:float;
                    end record;

   type indep_var is new var;

end aderiv;
