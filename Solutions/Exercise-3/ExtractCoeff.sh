COEFFS_OUTPUT=10
COEFF=$1
cat coeffs.data | awk -v Coeff=${COEFF} '$1==Coeff {print $2,$3}' | awk -v CO=${COEFFS_OUTPUT} '{printf "%5d   %s   %s\n",(NR-1)*CO,$1,$2}'
