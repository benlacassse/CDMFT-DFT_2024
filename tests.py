def read_lines_from_file(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            return lines
    except FileNotFoundError:
        print(f"The file {filename} does not exist.")
        return []

# Example usage
lines = read_lines_from_file('test.txt')
# for line in lines:
#     print(line)
_iclus = 0
p_str = []
nsym = True
print(len(lines))
for i in range(24):
    if lines[i] == '\n':
        current_line = i+1
        break
    elif "super-cluster" in lines[i]:
        _iclus +=1
    else:
        # the replacement here is necessary to differentiate
        # bath operators of each cluster within PyQCM
        if nsym:
            s = 1
            if "5_1" in lines[i]:
                p_str.append(lines[i][:4]+str(_iclus)+'='+str(s)+'*'+lines[i][:2]+'3_'+str(_iclus)+'\n')
            elif "6_1" in lines[i]:
                p_str.append(lines[i][:4]+str(_iclus)+'='+str(s)+'*'+lines[i][:2]+'4_'+str(_iclus)+'\n')
            else:
                p_str.append(lines[i].replace("_1", "_"+str(_iclus)))
        else:
            p_str.append(lines[i].replace("_1", "_"+str(_iclus)))

print(p_str)