#!/bin/bash

declare VALGRIND=0

for arg in $@; do
  if [ $arg = "v" ]; then
    VALGRIND=1
  fi
done

if [ ! -d temp ]; then
  mkdir temp
fi

cd temp 
cmake ..
make
cd ..

if [ $VALGRIND = 1 ]; then
  valgrind \
    --leak-check=full \
    --show-leak-kinds=all \
    --track-origins=yes \
    --verbose \
    --log-file=valgrind_temp.txt \
    ./temp/e2egpu_kc
  cat valgrind_temp.txt
  rm valgrind_temp.txt
else
  ./temp/e2egpu_kc
fi
