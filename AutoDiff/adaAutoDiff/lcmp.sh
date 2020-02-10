
export PATH=$HOME/opt/GNAT/2019/bin:$PATH

# use this to ensure a complete recompilation:
if [ -d ./obj/ ]; then
	rm ./obj/*
else
	mkdir obj
fi


gnatmake new1eg -O3 -gnat12 -D ./obj
gnatmake new2eg -O3 -gnat12 -D ./obj

gnatmake new3eg -O3 -gnat12 -D ./obj

gnatmake der4eg -O3 -gnat12 -D ./obj

