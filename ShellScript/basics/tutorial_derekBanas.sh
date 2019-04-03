#!/bin/bash

# Comment
# echo "Hello World"

# myName = "Vini"
# declare -r NUM1=5
# num2=4
# num3=$((NUM1+num2))
# num4=$((NUM1-num2))
# num5=$((NUM1*num2))
# num6=$((NUM1/num2))

# echo "5 + 4 = $num3"
# echo "5 - 4 = $num4"
# echo "5 * 4 = $num5"
# echo "5 / 4 = $num6"

# echo $((5**2))
# echo $(( 5%4 ))

# Assignment operators allow for shorthand arithmetic 
# +=, -=, *=, /=
# rand=5
# let rand+=4
# echo "$rand"

# Shorthand increment and decrement
# echo "rand++ = $(( rand++ ))"
# echo "++rand = $(( ++rand ))"
# echo "rand-- = $(( rand-- ))"
# echo "--rand = $(( --rand ))"

# Use Python to add floats
# num7=1.2
# num8=3.4
# num9=$(python -c "print($num7+$num8)")
# echo $num9

# You can print over multiple lines with a Here Script
# cat prints a file or any string past to it
# cat << END
# This text
# prints on
# many lines
# END

######################
# # Define function
######################

# getDate() {
    
#     # Get current date and time
#     date
    
#     # Return returns an exit status number between 0 - 255
#     return
# }

# getDate

# # This is a global variable
# name="Derek"

# # Local variable values aren't available outside of the function
# demLocal() {
#     local name="Paul"
#     return
# }

# demLocal

# echo "$name"

# # A function that receives 2 values and prints a sum
# getSum() {

#     # Attributes are retrieved by referring to $1, $2, etc.
#     local num3=$1
#     local num4=$2
    
#     # Sum values
#     local sum=$((num3+num4))
    
#     # Pass values back with echo
#     echo $sum
# }

# num1=5
# num2=6

# # You pass atributes by separating them with a space
# # Surround function call with $() to get the return value
# sum=$(getSum num1 num2)
# echo "The sum is $sum"

######################
# Conditionals / Input
######################

# You can use read to receive input which is stored in name
# The p option says that we want to prompt with a string
# read -p "What is your name? " name
# echo "Hello $name"

# read -p "How old are you? " age
    
# You place your condition with in []
# Include a space after [ and before ]
# Integer Comparisons: eq, ne, le, lt, ge, gt
# if [ $age -ge 16 ]
# then
#     echo "You can drive"

# # Check another condition
# elif [ $age -eq 15 ]
# then
#     echo "You can drive next year"
    
# # Executed by default
# else
#     echo "You can't drive"
    
# # Closes the if statement
# fi

# Extended integer test
# read -p "Enter a number : " num

# if ((num == 10)); then
#     echo "Your number equals 10"
# fi

# if ((num > 10)); then
#     echo "It is greater then 10"
# else
#     echo "It is less then 10"
# fi

# if (( ((num % 2)) == 0 )); then
#     echo " It is even"
# fi

# You can use logical operators like &&, || and !
# if (( ((num > 0)) && ((num < 11)) )); then
#   echo "$num is between 1 and 10"
# fi

# Create a file and then if that worked open it in Vim
#touch samp_file && vim samp_file

# If samp_dir doesn't exist make it
#[ -d samp_dir ] || mkdir samp_dir

# Delete file rm samp_file
# Delete directory rmdir samp_dir

# Testing strings
# str1=""
# str2="Sad"
# str3="Happy"

# # Test if a string is null
# if [ "$str1" ]; then
#     echo "$str1 is not null"
# fi

# if [ -z "$str1" ]; then
#     echo "str1 has no value"
# fi

# # Check for equality
# if [ "$str2" == "$str3" ]; then
#     echo "$str2 equals $str3"
# elif [ "$str2" != "$str3" ]; then
#     echo "$str2 is not equal to $str3"
# fi

# if [ "$str2" > "$str3" ]; then
#     echo "$str2 is greater then $str3"
# elif [ "$str2" < "$str3" ]; then
#     echo "$str2 is less then $str3"
# fi

# # Check the file test_file1 and test_file2
# file1="./test_file1"
# file2="./test_file2"

# if [ -e "$file1" ]; then
#     echo "$file1 exists"
    
#     if [ -f "$file1" ]; then
#         echo "$file1 is a normal file"
#     fi
    
#     if [ -r "$file1" ]; then
#         echo "$file1 is readable"
#     fi
    
#     if [ -w "$file1" ]; then
#         echo "$file1 is writable"
#     fi
    
#     if [ -x "$file1" ]; then
#         echo "$file1 is executable"
#     fi
    
#     if [ -d "$file1" ]; then
#         echo "$file1 is a directory"
#     fi
    
#     if [ -L "$file1" ]; then
#         echo "$file1 is a symbolic link"
#     fi
    
#     if [ -p "$file1" ]; then
#         echo "$file1 is a named pipe"
#     fi
    
#     if [ -S "$file1" ]; then
#         echo "$file1 is a network socket"
#     fi
    
#     if [ -G "$file1" ]; then
#         echo "$file1 is owned by the group"
#     fi
    
#     if [ -O "$file1" ]; then
#         echo "$file1 is owned by the userid"
#     fi
    
# fi

# With extended test [[ ]] you can use Regular Expressions
	
# read -p "Validate Date : " date

# pat="^[0-9]{8}$"

# if [[ $date =~ $pat ]]; then
#     echo "$date is valid"
# else
#     echo "$date is not valid"
# fi

# Read multiple values
	
# read -p "Enter 2 Numbers to Sum : " num1 num2

# sum=$((num1+num2))

# echo "$num1 + $num2 = $sum"

# # Hide the input with the s code
# read -sp "Enter the Secret Code" secret

# if [ "$secret" == "password" ]; then
#     echo "Enter"
# else
#     echo "Wrong Password"
# fi

# Use case to when it makes more sense then if

# read -p "How old are you : " age

# # Check the value of age
# case $age in

# # Match numbers 0 - 4
# [0-4]) 
#     echo "To young for school"
#     ;; # Stop checking further
    
# # Match only 5
# 5)
#     echo "Go to kindergarten"
#     ;;
    
# # Check 6 - 18
# [6-9]|1[0-8])
#     grade=$((age-5))
#     echo "Go to grade $grade"
#     ;;
    
# # Default action
# *)
#     echo "You are to old for school"
#     ;;
# esac # End case

# 8. Ternary Operator performs different actions based on a condition
# #!/bin/bash
# can_vote=0
# age=18

# ((age>=18?(can_vote=1):(can_vote=0)))
# echo "Can Vote : $can_vote"

######################
# Parameter Expansions and Strings
######################

# rand_str="A random string"
	
# # Get string length
# echo "String Length : ${#rand_str}"

# # Get string slice starting at index (0 index)
# echo "${rand_str:2}"

# # Get string with starting and ending index
# echo "${rand_str:2:7}"

# # Return whats left after A
# echo "${rand_str#*A }"
 
######################
# Looping
######################

# While
# num=1
	
# while [ $num -le 10 ]; do
#     echo $num
#     num=$((num + 1))
# done

# Continue and Break
# num=1

# while [ $num -le 20 ]; do

#     # Don't print evens
#     if (( ((num % 2)) == 0 )); then
#         num=$((num + 1))
#         continue
#     fi
    
#     # Jump out of the loop with break
#     if ((num >= 15)); then
#         break
#     fi
    
#     echo $num
#     num=$((num + 1))
# done

# Until loops until the loop is true
	
# num=1

# until [ $num -gt 10 ]; do
#     echo $num
#     num=$((num + 1))
# done

# Use read and a loop to output file info
	
# while read avg rbis hrs; do

#     # printf allows you to use \n
#     printf "Avg: ${avg}\nRBIs: ${rbis}\nHRs: ${hrs}\n"
    
# # Pipe data into the while loop
# done < barry_bonds.txt

# There are many for loop options. Here is the C form.
  	
# for (( i=0; i <= 10; i=i+1 )); do
#     echo $i
# done

# We can cycle through ranges
# for i in {A..Z}; do
#     echo $i
# done

######################
# Arrays
######################

# Messing with arrays
	
# Create an array
# fav_nums=(3.14 2.718 .57721 4.6692)

# echo "Pi : ${fav_nums[0]}"

# # Add value to array
# fav_nums[4]=1.618

# echo "GR : ${fav_nums[4]}"

# # Add group of values to array
# fav_nums+=(1 7)

# # Output all array values
# for i in ${fav_nums[*]}; do
#     echo $i;
# done

# # Output indexes
# for i in ${!fav_nums[@]}; do
#     echo $i;
# done

# # Get number of items in array
# echo "Array Length : ${#fav_nums[@]}"

# # Get length of array element
# echo "Index 3 length : ${#fav_nums[3]}"

# # Sort an array
# sorted_nums=($(for i in "${fav_nums[@]}"; do
#     echo $i;
# done | sort))

# for i in ${sorted_nums[*]}; do
#     echo $i;
# done

# # Delete array element
# unset 'sorted_nums[1]'

# # Delete Array
# unset sorted_nums

######################
# Positional Parameters
######################

# Print the first argument
echo "1st Argument : $1"

sum=0

# $# tells you the number of arguments
while [[ $# -gt 0 ]]; do

    # Get the first argument
    num=$1
    sum=$((sum + num))
    
    # shift moves the value of $2 into $1 until none are left
    # The value of $# decrements as well
    shift
done

echo "Sum : $sum"