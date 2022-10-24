# calculate minimum distance from COLVAR_SB file
awk '{if(NR==1)print "#! FIELDS time s1 s2 s3 s4 s5"; else{for(i=2;i<=NF;i+=2){s[i/2]=$i; j=i+1; if($j<$i) s[i/2]=$j}; if(NF==11) print $1,s[1],s[2],s[3],s[4],s[5]; if(NF==7) print $1,s[1],100.0,100.0,s[2],s[3]}}' COLVAR_SB > SALT_BRIDGE
# calculate salt bridge formation frequency with cutoff = 0.32 nm
for i in `seq 2 6`; do awk '{if(NR>1 && $i<0.32)tot++}END{print tot/(NR-1)}' i=$i SALT_BRIDGE; done > tmp
# prepare first column
cat << EOF > legenda
D30-K417K
E35-Q493K
D38-Q493K
K31-E35
D38-K353
EOF
# paste them together into SALT_BRIDGE.BAR
paste legenda tmp | awk '{printf "%9s %6.3lf\n",$1,$2}'> SALT_BRIDGE.BAR 
# clean up
rm SALT_BRIDGE tmp legenda
