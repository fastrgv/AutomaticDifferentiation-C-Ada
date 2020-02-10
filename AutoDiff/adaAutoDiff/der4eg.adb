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


with interfaces.c;
with mathtypes, text_io, autodiff;
use text_io;

procedure der4eg is

	use mathtypes;
	use interfaces.c;

	package fio is new text_io.float_io(real);
	use fio;

   package gradpak is new autodiff( dimension => 3 );
   use  gradpak;


   y       :  var;
   x       :  indep_vect_var;


	function f( x : indep_vect_var ) return var is
	  y : var;
	begin
		y := 1.0 + x(1) + x(2) + x(3)
			+ x(1)*x(2) + x(2)*x(3) + x(1)*x(3)
			+ x(1)*x(2)*x(3) + exp( x(1)/x(2) + x(2)/x(3) );
		return y;
	end f;


-------- begin main program -----------
begin

	-- set initial values
	set_indep_var( x, (1.0,2.0,3.0) );

	y := f(x);


	new_line;

	put("     u = ");
	fio.put( value(y), 5, 5, 0);
	new_line;

	put(" du/dx = ");
	fio.put( deriv(y,1), 5, 5, 0);
	new_line;

	put(" du/dy = ");
	fio.put( deriv(y,2), 5, 5, 0);
	new_line;

	put(" du/dz = ");
	fio.put( deriv(y,3), 5, 5, 0);
	new_line;


end der4eg;

