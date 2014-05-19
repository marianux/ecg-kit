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

export MATLABPATH='/home/bio/mllamedo/ecg_classification/a2hbc/'
#flag that the program does not ended correctly up to the moment.
export A2HBC_ES='1'
exec matlab  -nojvm -nodisplay -nosplash -r "$*"

if [ "$A2HBC_ES" -eq "1" ]
then
	#Some error happened.
	echo "Some error happened in $HOSTNAME. Error log should be above this." >/dev/fd/2
fi

#Program exit status (ES)
exit $A2HBC_ES
