dyn_name=$1
echo "a" > POSCAR
echo "1.0" >> POSCAR
sed -n '5,7p' $dyn_name >> POSCAR
atom1=`sed -n '8p' $dyn_name | awk '{gsub("\047","",$2); print $2}' `
atom2=`sed -n '9p' $dyn_name | awk '{gsub("\047","",$2); print $2}' `
atom3=`sed -n '10p' $dyn_name | awk '{gsub("\047","",$2); print $2}' `
echo $atom1 $atom2 $atom3 >> POSCAR
echo "1 2 24" >> POSCAR
echo "C" >> POSCAR
sed -n '11,37p' $dyn_name | awk {'print $3 "   " $4 "   " $5'}>> POSCAR


