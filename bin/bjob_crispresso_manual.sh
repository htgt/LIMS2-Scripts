#!/usr/local/bin/bash
PYTHONPATH=/nfs/team87/farm3_lims2_vms/software/python_local/lib/python/
PATH=/nfs/team87/farm3_lims2_vms/software/crispresso_dependencies/bin:$PATH
JOB="$(head -$LSB_JOBINDEX crispresso.sh | tail -1)"
echo $LSB_JOBINDEX
echo $JOB
$JOB
