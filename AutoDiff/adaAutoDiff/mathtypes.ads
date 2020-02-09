
with interfaces.c;
with ada.numerics.generic_elementary_functions;


package mathtypes is

	-- use this for serious work:
	subtype real is interfaces.c.double;

	-- use this to scrutinize error estimates:
	--subtype real is interfaces.c.c_float;

	package math_lib is new 
		ada.numerics.generic_elementary_functions(real);

end mathtypes;

