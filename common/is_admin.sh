cp -f $1 $2
bOk=$?
echo $bOk
rm -f $2
if [ $bOk -eq 0 -a $? -eq 0 ]; then
    exit 0
else
    exit 1
fi
