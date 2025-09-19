with open("Input.txt", "r") as input:
    lines = input.read().splitlines()
#print(lines)

length = len(lines)
#print(str(length))

i = 0
while i < (len(lines)):
    if (("#" in lines[i]) | ("" in lines[i])):
        lines.pop(i)
    i += 1

print(str(len(lines)))
print(lines)

#print(str(float(lines[2])))