#!/bin/bash

tim_file=$1
output_file="mc_output.txt"
var1=$2
var2=$3
var3=$4
var4=$5

chis=$(tempo2 -f temp.par $tim_file | grep "Chisqr/nfree" | awk '{print $9}') #the 9th entry on that line is the value of the chisqr/nfree
echo $var1 $var2 $var3 $var4 $chis >> $output_file

rm temp.par #Remove the temporary par file now that you're done using it
