#!/bin/tcsh 
# Note, all code in this scrpt is executed under tcsh
if (  "$SLURM_JOB_PARTITION" =~ *"chos"  ) then  
    echo  task-in-chos
    chosenv
    ls -l /proc/chos/link
else
    echo  task-in-shifter
    echo inShifter:`env|grep  SHIFTER_RUNTIME`
    cat /etc/*release
    #
    # - - - -  D O   N O T  T O U C H  T H I S   S E C T I O N- - - - 
    #
    whoami    
    echo  load STAR enviroment in shifter
    set NCHOS = sl64
    set SCHOS = 64
    set DECHO = 1
    set SCRATCH = $WRK_DIR/out-star1
    setenv GROUP_DIR /common/star/star${SCHOS}/group/
    source $GROUP_DIR/star_cshrc.csh    
    #
    # - - - -   Y O U   C A N   C H A N G E   B E L O W  - - - -
    #    
endif  

    cd ${WRK_DIR} # important
     
    echo  starting new-McLambdaDecay with execName=$EXEC_NAME, NUM_EVE=$NUM_EVE, workerName=`hostname -f`, startDate=`date`

    echo testing STAR setup $STAR_VER in `pwd`
    starver $STAR_VER 
    env |grep STAR

    echo 'my new STAR ver='$STAR'  test root '
    root -b -q 
    if ( $? != 0) then
	echo STAR environment is corrupted, aborting job
	echo $STAR
	which root
	exit
    endif
 
    #echo EEEEE ;   exit

    echo `date`" Fire: $EXEC_NAME for numEve=$NUM_EVE  [wiat]"
    cd ${WRK_DIR}/McEtaCut/McLambdaDecay
    echo `pwd`
    /usr/bin/time -v  $EXEC_NAME -b -q McLambdaEta.C+\($Energy,$PID,$Cent,$NUM_EVE,$SLURM_ARRAY_TASK_ID\) >& $LOG_PATH/Lambda_${SLURM_ARRAY_TASK_ID}.log

