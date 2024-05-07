#sage $1
echo "Running sage $1.sage...\n\n"
sage $1.sage

rm -rf $1.sage.py