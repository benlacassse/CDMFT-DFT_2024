import re
import copy
import numpy as np
import sys

class DMFTFileModifier:
    def __init__(self, file_name, layers=1, nodes=4):
        self.indmfl_name = file_name + '.indmfl'
        self.indmfi_name = file_name + '.indmfi'
        self.layers = layers
        self.nodes = nodes
        self.lines = []
        self.content = ""
        self.modified_content = ""
        self.old_dimensions = 0
        self.ind_components = 0
        self.cluster_matrix = []

    def read_file(self, file_name):
        try:
            with open(file_name, 'r') as file:
                self.lines = file.readlines()
            with open(file_name, 'r') as file:
                self.content = file.read()
        except FileNotFoundError:
            raise FileNotFoundError(f"The file '{file_name}' was not found.")
        except Exception as e:
            raise Exception(f"An error occurred while reading the file: {e}")

    def validate_layers(self):
        if self.layers == 1:
            raise ValueError('No changes due to no additional layers.')
        if self.layers > 5:
            raise ValueError('Number of layers too high for the code at the moment.')

    def replace_cix(self):
        cix_count = 0
        modified_lines = copy.deepcopy(self.lines)  
        for i, line in enumerate(self.lines):
            if 'cix' in line:
                cix_count += 1
                count = 0
                matches = re.findall(r'\d+', line)
                for match in re.finditer(r'\d+', line):
                    count += 1
                    if count == 3:
                        position = match.start()
                        break
                if cix_count <= self.layers*self.nodes:
                    modified_lines[i] = line[:position] + ' '*(len(matches[2])-1) + '1' + line[position + len(matches[2]):]
                elif cix_count <= self.layers*self.nodes*2:
                    modified_lines[i] = line[:position] + ' '*(len(matches[2])-1) + str(int(matches[2])-self.layers+1) + line[position + len(matches[2]):]
                else:
                    break
        self.lines = modified_lines
        self.modified_content = "".join(modified_lines)
        
    def replace_parameters(self):
        search_line = False
        modify_next_line = False
        modified_lines = copy.deepcopy(self.lines)
        for i, line in enumerate(self.lines):
            if modify_next_line:
                matches = re.findall(r'\d+', line)
                if matches:
                    new_max_dimensions = str(self.old_dimensions * self.layers)
                    line = self.replace_number(line, str(self.old_dimensions), new_max_dimensions)
                    new_max_ind_components = str(self.ind_components * (self.generate_diagonal_pattern_matrix()[1]+self.layers//2 + self.layers%2))
                    line = self.replace_number(line, str(self.ind_components), new_max_ind_components)
                modified_lines[i] = line
                modify_next_line = False

            if search_line:
                matches = re.findall(r'\d+', line)
                if matches:
                    kcix_blocks = matches[0]
                    if int(kcix_blocks) != (self.nodes + 1) * self.layers:
                        raise ValueError('Input of layers and nodes do not match the number of independent kcix blocks.')
                    else:
                        new_kcix_blocks = str(int(kcix_blocks) - self.layers + 1)
                        if len(new_kcix_blocks) < len(kcix_blocks):
                            new_kcix_blocks = new_kcix_blocks.ljust(len(kcix_blocks))
                        line = self.replace_number(line, kcix_blocks, new_kcix_blocks)
                        
                        self.old_dimensions = int(matches[1])
                        new_max_dimensions = str(self.old_dimensions * self.layers)
                        line = self.replace_number(line, str(self.old_dimensions), new_max_dimensions)

                        self.ind_components = int(matches[2])
                        new_max_ind_components = str(self.ind_components * (self.generate_diagonal_pattern_matrix()[1]+self.layers//2 + self.layers%2))
                        line = self.replace_number(line, str(self.ind_components), new_max_ind_components)
                    modified_lines[i] = line

                search_line = False
                modify_next_line = True

            if line.strip().startswith('#='):
                search_line = True
        self.lines = modified_lines
        self.modified_content = "".join(modified_lines)

    def replace_number(self, line, old, new):
        return re.sub(r'\b{}\b'.format(re.escape(old)), new, line, 1)

    def independent_components(self):
        modified_lines = copy.deepcopy(self.lines)
        line_count = 0
        found_line = False
        ind_comp_list = []
        for line in self.lines:
            if line.strip().startswith('#---------------- # Independent'):
                line_count += 1
                found_line = True
                continue
            if found_line:
                if line_count > self.layers:
                    break
                if line.rstrip() not in ind_comp_list:
                    ind_comp_list.append(line.rstrip())
                found_line = False
        new_ind_comp = ' '.join(ind_comp_list)
        num_comp = self.ind_components * self.generate_diagonal_pattern_matrix()[1]
        found_line = False
        for i in range(num_comp):
            x = i+self.ind_components*len(ind_comp_list)+1
            new_ind_comp += ' ' + f"'x2y2-{x}'"
        for i, line in enumerate(self.lines):
            if line.strip().startswith('#---------------- # Independent'):
                found_line = True
                continue
            if found_line:
                modified_lines[i] = new_ind_comp + '\n'
                break
            
        self.lines = modified_lines
        self.modified_content = "".join(modified_lines)

    def modify_line(self, line):
        parts = line.split()
        
        modified_line = []
        add_value = self.ind_components * self.generate_diagonal_pattern_matrix()[1]

        for part in parts:
            if part.isdigit():
                number = int(part)
                if number != 0:
                    modified_number = number + add_value
                    modified_line.append(str(modified_number))
                else:
                    modified_line.append(part)
            else:
                raise ValueError("The line contains invalid characters.")

        return ' '.join(modified_line)

    def create_sigind_matrix(self):
        matrix_data = []
        sigind_found = False
        j = 0
        for line in self.lines:
            if 'Sigind ' in line:
                sigind_found = True
                continue
            if sigind_found == True : 
                if j < self.old_dimensions:
                    row = list(map(int, line.split()))
                    matrix_data.append(row)
                    j += 1

        matrix_array = np.array(matrix_data)
        maximum = 0
        for row in matrix_array:
            if max(row) > maximum:
                maximum = max(row)

        matrix_array_trimmed = matrix_array[~np.all(matrix_array == 0, axis=1)]
        matrix_array_trimmed = matrix_array_trimmed[:, ~np.all(matrix_array_trimmed == 0, axis=0)]
        matrix_array_trimmed = matrix_array_trimmed.reshape(1, 1, 4, 4)

        coefficients_matrix = self.generate_diagonal_pattern_matrix()[0]
        coefficients_matrix = (maximum * coefficients_matrix).reshape(self.layers, self.layers, 1, 1)

        matrix_4d = matrix_array_trimmed + coefficients_matrix
        return matrix_4d

    def replace_sigind_matrix(self):
        modified_lines = copy.deepcopy(self.lines)
        sigind_count = 0
        
        for i, line in enumerate(self.lines):
            this_line = line.rstrip()
            if 'Sigind ' in line:
                sigind_count += 1
                j=0
                continue
            if sigind_count == 1 : 
                if j < self.old_dimensions:
                    if line == '0 '*(self.old_dimensions-1)+'0\n':
                        for k in range(self.layers):
                            new_line = line.rstrip()+' '
                            if k == 0:
                                modified_lines[i] = new_line*self.layers+'\n'
                            else:
                                modified_lines.insert(i + self.old_dimensions,  new_line*self.layers+'\n')
                            
                    else:
                        for k in range(self.layers):
                            new_line = ''
                            for l in range(self.layers):
                                for m in range(self.nodes):
                                    new_line = new_line + str(self.create_sigind_matrix()[k][l][j][m]) + ' '
                            
                            if k == 0:
                                modified_lines[i] = new_line + '0 '*(self.old_dimensions-self.nodes)*self.layers+'\n'
                            else:
                                modified_lines.insert(i + self.old_dimensions+(k-1)*(j+1), new_line + '0 '*(self.old_dimensions-self.nodes)*self.layers+'\n')
                    j += 1

            elif sigind_count > 1:
                if j < self.old_dimensions/self.nodes:
                    this_line = line.rstrip()
                    modified_lines[i + self.old_dimensions*(self.layers-1)] = self.modify_line(this_line) +'\n'
                    j += 1


        self.modified_content = "".join(modified_lines)
        self.lines = modified_lines

    def generate_diagonal_pattern_matrix(self):
        n = self.layers

        diagonal_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i == j:
                    diagonal_matrix[i, j] = 1

        return diagonal_matrix, n

    def update_file(self, file_name):
        with open(file_name, 'w') as file:
            file.write(self.modified_content)
            file.write('\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script_name.py <file_name> <layers> <nodes>")
        sys.exit(1)
    
    file_name = sys.argv[1]
    layers = int(sys.argv[2])
    nodes = int(sys.argv[3])

    modifier = DMFTFileModifier(file_name, layers, nodes)
    try:
        modifier.read_file(modifier.indmfl_name)
        modifier.read_file(modifier.indmfi_name)
        modifier.validate_layers()
        modifier.replace_cix()
        modifier.replace_parameters()
        modifier.independent_components()
        modifier.replace_sigind_matrix()
        modifier.update_file(modifier.indmfl_name)
        modifier.update_file(modifier.indmfi_name)
        print("Files updated successfully.")
    except Exception as e:
        print(f"An error occurred: {e}")
