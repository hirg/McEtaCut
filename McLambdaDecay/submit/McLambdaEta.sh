#!/bin/bash
date

# . ./McLambdaEta.sh

if [ $# -eq 0 ]
then
  PID="Lambda"
  Energy="200GeV"
  Name="Mc"$PID$Energy
  for((counter=0;counter<100;counter=counter+1))
  do
    cp ./run.csh ./run$Name${counter}.csh

    echo "mkdir $Name$counter" >> ./run$Name${counter}.csh
    echo "cp McLambdaEta.C $Name$counter" >> ./run$Name${counter}.csh
    echo "cd $Name$counter" >> ./run$Name${counter}.csh

    echo "root4star -b -q -x McLambdaEta.C++'('6','0','0','5000000','$counter')'" >> ./run$Name${counter}.csh
    echo "cd .." >> ./run$Name${counter}.csh
    echo "rm -rf $Name$counter" >> ./run$Name${counter}.csh
    qsub -hard -l scratchfree=500,h_cpu=2:00:00,h_vmem=1.8G,projectio=1 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/job$Name$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/job$Name$counter.err ./run$Name$counter.csh
    mv ./run$Name${counter}.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Script/
  done
else
  echo "Wrong number of parameters!!"
fi

