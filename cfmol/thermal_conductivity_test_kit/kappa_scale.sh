#!/bin/bash
#This script scales thermal conductivity with respect to a bulk reference

echo "enter kappa: "
read -e kappa
echo "enter kappaerror: "
read -e kappaerror

echo "Is scaling value 196.8 +/- 20 W/m-K? (1-y/0-n)"
read -e answer

if [ "$answer" = "0" ]; then
    echo "enter scaling kappa: "
    read -e kappascaling
    echo "enter uncertainty: "
    read -e kappascaleerror
else
    kappascaling=196.8
    kappascaleerror=20
fi

kappascaled=$(bc <<< "scale = 10; ($kappa / $kappascaling)")

kappascalederror=$(bc <<< "scale = 10; $kappascaled*(sqrt((($kappaerror/$kappa)^2)+(($kappascaleerror/$kappascaling)^2)))")

echo $kappascaled "+/-" $kappascalederror

