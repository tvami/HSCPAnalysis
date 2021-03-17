input=$1
joblist=$2

sed "s;JOBLIST;$joblist;;s;SANDBOX;$input;;" condor_task_template_EU.cfg > tmp.cfg

condor_submit tmp.cfg
rm tmp.cfg
