repdirs=$(ls -1d rep??)
nrep=$(echo $repdirs | awk '{ print NF }')
dir=$(ls -1d rep?? | head -n 1)

for file in $dir/run.log ; do

grep "Repl pr" $file | tail -n +50 |
awk -v file="$file" -v nrep="$nrep" '{if(NF%2==1)j=1; else j=0; for(i=0;i<=NF-3;i++) {
  count[2*i+j]++;
  prob[2*i+j]+=$(i+3);
  }
  } END{
    print("```")
    print(file)
    print("exchanges:",count[0])
    for(i=0;i<nrep-1;i++) print(i,"x",i+1,prob[i]/count[i])
    print("```")
  }'

done
