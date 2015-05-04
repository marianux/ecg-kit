#!/bin/sh
#
# Filename: matlab.sh
#
# ./matlab.sh "simplesleep(10,2,'output.txt')"
#

on_die()
{
	#Some error happened.
	echo "Some error happened in $HOSTNAME. Error log should be above this." >/dev/fd/2
	exit 1	
}


# Execute function on_die() receiving TERM signal
#
trap 'on_die' TERM INT

# run matlab in text mode 

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'/home/bio/mllamedo/wfdb/lib64'
export PATH=$PATH:'/home/bio/mllamedo/wfdb/bin:/home/bio/mllamedo/ari/bin:/home/bio/mllamedo/ecg-kit/common'
export MATLABPATH='/home/bio/mllamedo/ecg-kit/:/home/bio/mllamedo/ecg-kit/common/:/home/bio/mllamedo/ecg-kit/examples/'
#flag that the program does not ended correctly up to the moment.
export A2HBC_ES='1'
exec /extra/software/bin/matlab  -nodesktop -nodisplay -nosplash -r "$*"

if [ "$A2HBC_ES" -eq "1" ]
then
	#Some error happened.
	echo "Some error happened in $HOSTNAME. Error log should be above this." >/dev/fd/2
fi

#Program exit status (ES)
exit $A2HBC_ES
