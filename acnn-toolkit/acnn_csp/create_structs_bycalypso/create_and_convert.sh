rm -r results
rm caly.log POSCAR* *py step
calypso.x > caly.log 2>&1


a=$(basename $(pwd))

wait

for i in {1..5000}; do
echo $i
sed '6 i Li Pb H' POSCAR_$i > POSCAR-$i
cabal poscar res < POSCAR-$i > ../../Base/LiPbH-${i}-${a}-N64.res
rres ../../Base/LiPbH-${i}-${a}-N64.res -s LiPbH-${i}-${a}-N64-st
cabal res res 0 < ../../Base/LiPbH-${i}-${a}-N64.res > ../../Base/LiPbH-${i}-${a}-N64-st.res
rm ../../Base/LiPbH-${i}-${a}-N64.res

done
