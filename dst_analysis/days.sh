#!/bin/bash
#. ~/wait_jobs.sh
for (( day_number=68; day_number<89;  day_number++ ))
do
     . ./sendScript_SL.sh 0$day_number
done

. ~/wait_jobs.sh
