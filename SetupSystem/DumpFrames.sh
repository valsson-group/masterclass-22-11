GMX_BIN=gmx_mpi

for i in `seq 1 10`
do 
  time=`echo "${i}*200" | bc -l`
  echo "0" | ${GMX_BIN} trjconv -f NaCl_NPT.trr -s NaCl_NPT.tpr -o NaCl_StartingStructure-${i}.gro -dump ${time} 
done 











