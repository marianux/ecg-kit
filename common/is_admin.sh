if [ "$(id -u)" -eq "0" ]; then
	echo "root user"
	exit 0
fi
echo "You must be a root user"
exit 1

