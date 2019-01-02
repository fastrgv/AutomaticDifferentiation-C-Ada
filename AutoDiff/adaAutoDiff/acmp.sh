
# expose AdaCore 2018 compiler's default location:
export PATH=$HOME/opt/GNAT/2018/bin:$PATH

gnatmake $1 -O3 -gnat12  --subdirs=./obj

cp ./obj/$1 .

