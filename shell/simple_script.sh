#!/bin/ksh

if [[ $# -eq 0 ]];then
   print "No Arguments"
fi

print $*
print $1

value=7

if [[ $value -eq 7 ]];then
   print "$value is 7"
fi

count=10

while [[ $count -gt 0 ]];do
   print "\$count is $count"
   (( count -= 1 ))
done

typeset -Z4 X
X=1
print "X is $X"

EXP=`cat new_file`
print "EXP is $EXP"

#until [[ $answer = "yes" ]];do
#   print -n "Please enter \"yes\": "
#   read answer
#   print ""
#done


for foo in $(ls);do
   if [[ -d $foo ]];then
      print "$foo is a directory"
   else
      print "$foo is not a directory"
   fi
done

cat > new_file << EOF
#!/bin/ksh
#This is an empty script created
#by simple_script.sh, coooooool...
EOF
