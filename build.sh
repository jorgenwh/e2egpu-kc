if [ ! -d temp ]; then
  mkdir temp
fi

cd temp 
cmake ..
make
cd ..
./temp/e2egpu_kc
