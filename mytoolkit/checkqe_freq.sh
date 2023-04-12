#!/bin/bash
echo "你需要指定两个参数，比如你有10个目录，名字分别为1，2，3，4，5，6，7，8，9，10"
echo "那么你就这样使用该脚本：checkqe_frep.sh 1 10"

begin=$1
end=$2
for ((i=$begin; i<=$end; i++)); do
#for i in {$beign..$end}; do

cd $i
echo $i
grep freq *dyn | head 
cd ..

done
