START_POP=11
ENS_NUM=1000



DYN_POP="dyn_pop"$START_POP"_"
END_POP=$((START_POP+1))
FORCE_NAME="forces_population"$END_POP"_{}.dat"
STRESS_NAME="pressures_population"$END_POP"_{}.dat"
ENER_NAME="energies_supercell_population"$END_POP".dat"
DYN_POP="dyn_pop"$START_POP"_"
OUT_POP="OUT"$END_POP".dat"
cat 6_ParseOutput.py | sed -e "s/FORCENAME/$FORCE_NAME/g" \
                  | sed -e "s/STRESSNAME/$STRESS_NAME/g" \
                  | sed -e "s/ENERNAME/$ENER_NAME/g" \
                  > R6_ParseOutput.py

cat 0_Relax.py | sed -e "s/DYNPOP/$DYN_POP/g" \
	       | sed -e "s/ENDPOP/$END_POP/g" \
	       | sed -e "s/ENSNUM/$ENS_NUM/g" \
               > Relax.py
cat 7_SubRelax.sh | sed -e "s/OUTPOP/$OUT_POP/g" \
	          > R7_SubRelax.sh
python R6_ParseOutput.py
#
sleep 5
#
sbatch R7_SubRelax.sh
