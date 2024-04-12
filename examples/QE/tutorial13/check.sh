#!/bin/bash

read -p "Enter a number: " num

if [[ $num =~ \. ]]; then
    num=1
    echo "Converted to integer: $num"
fi


