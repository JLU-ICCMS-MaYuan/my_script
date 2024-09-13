rm V3Hessian.dyn0
echo "2 2 2" >> V3_Hessian.dyn0 
echo "4" >> V3_Hessian.dyn0 
for i in `seq 1 4`
do
grep "q =" V3_Hessian.dyn$i | head -1 | awk '{print $4i"  "$5"  "$6}' >> V3_Hessian.dyn0
done
source /work/env/intel2018
~/soft/q-e-qe-7.0/bin/q2r.x < q2r.in > q2r.out
~/soft/q-e-qe-7.0/bin/matdyn.x < matdyn.in > matdyn.out

