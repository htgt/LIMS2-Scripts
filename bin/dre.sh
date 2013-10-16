#! /bin/bash

# apply dre recombinase appropriately for plates <= CEP 26 and > 26

printf "Modifying Dre recombinase for CEP plates\n\n"

printf "Adding Dre recombinase to CEP plates <= 26...\n"

add_dre_plate_names=(
"CEP00001"
"CEP00002"
"CEP00003"
"CEP00005"
"CEP00006"
"CEP00007"
"CEP00008"
"CEP00009"
"CEP00010"
"CEP00011"
"CEP00012"
"CEP00013"
"CEP00014"
"CEP00015"
"CEP00016"
"CEP00017"
"CEP00018"
"CEP00019"
"CEP00020"
"CEP00021"
"CEP00022"
"CEP00023"
"CEP00024"
"CEP00025"
"CEP00026"
)

plate_counter=0
for plate in ${add_dre_plate_names[*]}
do
    printf "\n==> Plate No. %d: %s\n\n" `expr $plate_counter + 1` $plate
    perl add_dre.pl $plate
    let plate_counter++
done

printf "Removing Dre recombinase from CEP plates > 26...\n"

remove_dre_plate_names=(
"CEP00027"
"CEP00028"
"CEP00029"
"CEP00030"
"CEP00031"
"CEP00032"
"CEP00033"
"CEP00034"
"CEP00035"
"CEP00036"
"CEP00037"
"CEP00038"
"CEP00039"
"CEP00040"
"CEP00041"
"CEP00042"
"CEP00043"
"CEP00044"
"CEP00045"
"CEP00046"
"CEP00047"
"CEP00048"
"CEP00049"
"CEP00050"
"CEP00053"
"CEP00054"
"CEP00055"
"CEP00056"
"CEP00058"
"CEP00059"
)

plate_counter=0
for plate in ${remove_dre_plate_names[*]}
do
    printf "\n==> Plate No. %d: %s\n\n" `expr $plate_counter + 1` $plate
    perl remove_dre.pl $plate
    let plate_counter++
done


