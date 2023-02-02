#!/bin/bash

#shopt -s lastpipe # enable lastpipe

  for parN in 1.0 1.2 1.4 1.5 1.6 1.8 2.0
  do
  for (( iIter=0; iIter<$nInterations; iIter++ ))
    do
     echo par = ${parN} iteration = ${iIter}
   done
done
