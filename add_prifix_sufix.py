import sys

list_of_lists = []

for line in sys.stdin:
    new_list = [elem for elem in line.split()]
    list_of_lists.append(new_list)

input_list=list_of_lists[0]

output_list=['done('+a+")" for a in input_list]
output_string=" && ".join(output_list)

print output_string

