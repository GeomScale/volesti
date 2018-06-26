#!/bin/bash

cmake .
make

if [[ $1 == "-exp" ]]; then
./vol -f1 ./test_data/birk3.ine -mag 1
./vol -f1 ./test_data/birk4.ine -mag 4
./vol -f1 ./test_data/birk5.ine -mag 7
./vol -f1 ./test_data/birk6.ine -mag 13
./vol -f1 ./test_data/birk7.ine -mag 21
#----------------------------------#
./vol -f1 ./test_data/cross_2.ine -mag 0
./vol -f1 ./test_data/cross_3.ine -mag 0
./vol -f1 ./test_data/cross_4.ine -mag 1
./vol -f1 ./test_data/cross_5.ine -mag 1
#--------------------------------#
./vol -f1 ./test_data/cube10.ine -mag 0
./vol -f1 ./test_data/cube11.ine -mag 0
./vol -f1 ./test_data/cube12.ine -mag 0
./vol -f1 ./test_data/cube13.ine -mag 0 
./vol -f1 ./test_data/cube14.ine -mag 0
./vol -f1 ./test_data/cube15.ine -mag 0
./vol -f1 ./test_data/cube16.ine -mag 0
./vol -f1 ./test_data/cube17.ine -mag 0
./vol -f1 ./test_data/cube18.ine -mag 0
./vol -f1 ./test_data/cube19.ine -mag 0
./vol -f1 ./test_data/cube20.ine -mag 0
./vol -f1 ./test_data/cube30.ine -mag 0
./vol -f1 ./test_data/cube40.ine -mag 0
#--------------------------------#
./vol -f1 ./test_data/prod_simplex_5_5.ine -mag 5
./vol -f1 ./test_data/prod_simplex_10_10.ine -mag 14
./vol -f1 ./test_data/prod_simplex_15_15.ine -mag 25
./vol -f1 ./test_data/prod_simplex_20_20.ine -mag 37
#--------------------------------------------#
./vol -f1 ./test_data/simplex10.ine -mag 7
./vol -f1 ./test_data/simplex20.ine -mag 19
./vol -f1 ./test_data/simplex30.ine -mag 33
./vol -f1 ./test_data/simplex40.ine -mag 49
./vol -f1 ./test_data/simplex50.ine -mag 65
#----------------------------------#
./vol -f1 ./test_data/skinny_cube2.ine -mag 0
./vol -f1 ./test_data/skinny_cube3.ine -mag 0
./vol -f1 ./test_data/skinny_cube4.ine -mag 0
./vol -f1 ./test_data/skinny_cube5.ine -mag 0
fi
