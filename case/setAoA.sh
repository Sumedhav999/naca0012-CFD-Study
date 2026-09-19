#!/bin/bash
A=$1
UX=$(python3 -c "import math;print(f'{30*math.cos(math.radians($A)):.6f}')")
UY=$(python3 -c "import math;print(f'{30*math.sin(math.radians($A)):.6f}')")
DX=$(python3 -c "import math;print(f'{math.cos(math.radians($A)):.6f}')")
DY=$(python3 -c "import math;print(f'{math.sin(math.radians($A)):.6f}')")
sed -i "s|^internalField uniform .*|internalField uniform ($UX $UY 0);|" 0/U
sed -i "s|^\( *dragDir *\).*|\1( $DX  $DY 0);|" system/forceCoeffs
sed -i "s|^\( *liftDir *\).*|\1(-$DY  $DX 0);|" system/forceCoeffs
echo "AoA=$A  U=($UX $UY 0)  drag=($DX $DY 0)  lift=(-$DY $DX 0)"
