#/usr/bin/sh
prog=clique

for j in  "c1.txt" "c2.txt" "c3.txt" "c4.txt" "c5.txt" "c6.txt" "c7.txt" "c8.txt" "c9.txt" "c10.txt"
do
 date > results_$j

 for i in  1 2 3 4 5 6 7 8 9 10
 do
     echo Executing run $i of $prog with configuration file $j 
     ./$prog $j >> results_$j
     echo "************************************" >> results_$j
 done
done