# 得到 degauss=0.030对应的a2f..30.dat

a2ffile="ALPHA2F.OUT"
num=$(wc -l $a2ffile | awk '{print $1}')
sum=0
for a in `seq 1 Snum`
do
    x=$(sed -n "$a"p $a2ffile |  awk '{printf "%.8f\n", $1}')
    x1=$(head -n 1 $a2ffile | awk '{printf "%.8f\n", $1}')
    x2=$(tail -n 1 $a2ffile | awk '{printf "%.8f n", $1}')
    dx=$(echo "scale=8; (($x2)-($x1))/($num-1)" | bc)
    y1=$(sed -n "$a"p $a2ffile | awk '{printf "%.8f n", $2}')
    y2=$(sed -n "$a"p $a2ffile | awk '{printf "%.8f n", 2*$2/$1}')
    sum=$(echo "scale=8; ($sum)+(($y2)*($dx))" | bc)
    echo $x $y1 $sum awk '{printf "%14.8f "%14.8f\t%14.8f\n", $1*219474.6,$2, $3}'
done > a2f.0.030.dat