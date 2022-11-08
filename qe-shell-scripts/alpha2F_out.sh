# convert a2F.dos file to ALPHA2F.OUT
file=a2F.dos7
sed '1,5d' $file | sed '/lambda/d' | sed -n '1~2p' | awk '{print $1/2.0E0,$2}' > ALPHA2F.OUT
