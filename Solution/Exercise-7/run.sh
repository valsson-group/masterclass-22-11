GMX_BIN=gmx_mpi
NumProcs=2
RunFilename=NaCl_VES_NPT-300K

mpirun -np ${NumProcs} ${GMX_BIN}  mdrun -deffnm ${RunFilename} -v -plumed plumed.dat
echo "0" | ${GMX_BIN} trjconv -f ${RunFilename}.xtc -s ${RunFilename}.tpr -pbc whole -o ${RunFilename}.pbc-whole.xtc


