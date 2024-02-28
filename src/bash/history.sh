#!/bin/bash
echo "-----------------------------------------------------"
echo "Showing latest 10 mainprogram command"
echo "-----------------------------------------------------"

grep "mainprogram " ~/.bash_history | sed '/history/d' | sed '/vi/d'| tail -n 10
