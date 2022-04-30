#!/bin/bash

script=sendScript_SL.sh


for (( day_number=67; day_number<89;  day_number++ ))
     do
	 . ${script} $day_number accepted ag158ag_3200A
	 . ${script} $day_number rejected ag158ag_3200A
done


for (( day_number=88; day_number<91;  day_number++ ))
     do
	 . ${script} $day_number accepted ag123ag_2500A
	 . ${script} $day_number rejected ag123ag_2500A
done

for (( day_number=90; day_number<92;  day_number++ ))
     do
	 . ${script} $day_number accepted ag123ag_200A
	 . ${script} $day_number rejected ag123ag_200A
done
