#!/bin/csh

@ day = 244

#cp input.nml inputDOY244.nml

while ( $day <= 273 ) 

 sed s/244/${day}/ inputDOY244.nml >! input.nml
 convert_aura

 @ day = $day + 1

end
