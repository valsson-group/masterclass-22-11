GMX_BIN=gmx_mpi
StartingGeometry=NaCl_StartingStructure-1.gro
RunFilename=NaCl_VES_NPT-300K
${GMX_BIN}  grompp -f MD-NPT.mdp -c ${StartingGeometry} -p NaCl.top -o ${RunFilename}.tpr  -maxwarn 1
