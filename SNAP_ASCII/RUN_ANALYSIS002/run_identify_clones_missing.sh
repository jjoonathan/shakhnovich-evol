#!/bin/bash
#BSUB -J missing_rerun
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -q normal_serial


grep -e clones missing.e | awk '{print "RUN1."$2}' | sed 's/clones:/snap.002.a/' > list.dat
                                      
for SNAP in `cat list.dat`; do
  RUN=`ls ${SNAP}| awk 'BEGIN{FS="."}{print $2}'`
  echo ${RUN}

  rm ${RUN}.clones
  for i in `seq 0 9`; do
    echo ${i}
    grep -e "G ${i} " ${SNAP} | awk '{print $9}' | sort | uniq > temp1.uniq
 
    rm temp1.clones
    for j in `cat temp1.uniq`;do 
      FREQ=`grep -c ${j} ${SNAP}`;
      echo ${FREQ} ${j} >> temp1.clones
    done
 
    #dominant clones
    sort -n -r temp1.clones | awk '{printf "%d %d %s\n",'${i}', $1, $2}' >> ${RUN}.clones 
  done
  echo "..."

done
