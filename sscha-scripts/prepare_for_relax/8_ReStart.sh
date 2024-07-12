START_POP=11 #计算第0代的动力学矩阵
ENS_NUM=1000



DYN_POP="dyn_pop"$START_POP"_"
END_POP=$((START_POP+1))
SCF_NAME="scf_population"$END_POP"_"
cat 1_CreatEns.py | sed -e "s/DYNPOP/$DYN_POP/g" \
	          | sed -e "s/ENSNUM/$ENS_NUM/g" \
		  | sed -e "s/POPID/$END_POP/g" \
		  > R1_CreatEns.py

cat 2_CreatQEInput.py | sed -e "s/SCFNAME/$SCF_NAME/g" \
	              > R2_CreatQEInput.py

#conda activate sscha

python R1_CreatEns.py

sleep 5

python R2_CreatQEInput.py
sleep 5
python 3_Creat_Sub.py
sleep 5
. 4_SubAllJobs.sh
#
#
#
#
