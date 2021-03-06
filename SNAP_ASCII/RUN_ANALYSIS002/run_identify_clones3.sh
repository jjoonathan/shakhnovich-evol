#!/bin/bash
#BSUB -J identifyclones
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -q normal_serial

for SNAP in `seq --format %010g 100501 500 150001 |awk '{print "RUN1.gen" $0 ".snap.002.a"}'`; do  
                                      
  RUN=`ls ${SNAP}| awk 'BEGIN{FS="."}{print $2}'`
  echo ${RUN}

  rm ${RUN}.clones
  for i in `seq 0 9`; do
    echo ${i}
    grep -e "G ${i} " ${SNAP} | awk '{print $9}' | sort | uniq > temp3.uniq
 
    rm temp3.clones
    for j in `cat temp3.uniq`;do 
      FREQ=`grep -c ${j} ${SNAP}`;
      echo ${FREQ} ${j} >> temp3.clones
    done
 
    #dominant clones
    sort -n -r temp3.clones | awk '{printf "%d %d %s\n",'${i}', $1, $2}' >> ${RUN}.clones 
  done
  echo "..."

done
