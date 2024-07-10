import re
import copy
import numpy as np

class DMFTFileModifier:
    def __init__(self, file_name, layers=1, nodes=4):
        self.file_name = file_name
        self.layers = layers
        self.nodes = nodes
        self.lines = []
        self.content = ""
        self.modified_content = ""
        self.old_dimensions = 0
        self.ind_components = 0

    def read_file(self):
        try:
            with open(self.file_name, 'r') as file:
                self.lines = file.readlines()
            with open(self.file_name, 'r') as file:
                self.content = file.read()
        except FileNotFoundError:
            raise FileNotFoundError(f"The file '{self.file_name}' was not found.")
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
                    modified_lines[i] = line[:position] + ' '*(len(matches[2])-1) + str(int(matches[2])-layers+1) + line[position + len(matches[2]):]
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
                    line = self.replace_number(line, str(self.old_dimensions), new_max_dimensions)
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
                        new_max_ind_components = str(self.ind_components * self.layers)
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
        search_line = False
        modified_lines = copy.deepcopy(self.lines)
        for i, line in enumerate(self.lines):
            if search_line:
                line = (line.rstrip() + ' ')*self.layers + '\n'
                modified_lines[i] = line
                break
            if line.strip().startswith('#---------------- # Independent'):
                search_line = True 
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

    def block_finder(self):
        search_line = False
        block_ind = 0
        block_list = [] 
        for line in self.lines:
            if 'Sigind ' in line:
                search_line = True
                continue
            if search_line == True :
                if block_ind <  self.nodes:
                    j = 0
                    for char in line.strip():
                        j +=1
                        if char.isdigit():
                            if char != '0':
                                line_block = line[:j]
                    block_list.append(line_block)
                    block_ind += 1
        return block_list

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
                    # PROBLEM!!!!


        self.modified_content = "".join(modified_lines)
        self.lines = modified_lines

    def generate_diagonal_pattern_matrix(self):
        n = self.layers
        if n == 2:
            matrix = np.array([[1, 2],[2, 1]])
            jump = 1
        elif n == 3:
            matrix = np.array([[1, 4, 3],[4, 1, 3], [3, 3, 2]])
            jump = 2
        elif n == 4:
            matrix = np.array([[1, 6, 4, 5],[6, 1, 5, 4], [4, 5, 2, 3], [5, 4, 3, 2]])
            jump = 4
        elif n == 5:
            matrix = np.array([[1, 9, 7, 8, 6],[9, 1, 8, 7, 6], [7, 8, 2, 5, 4], [8, 7, 5, 2, 4], [6, 6, 4, 4, 3]])
            jump = 6

        for i, row in enumerate(matrix):
            for j, coefficient in enumerate(row):
                matrix[i][j] -= 1

        return matrix, jump

    def process_lines(self, lines):
        value_lines = []
        zero_lines = []
        
        for line in lines:
            stripped_line = line.replace(' ', '')
            if all(ch == '0' for ch in stripped_line):
                zero_lines.append(line)
            else:
                value_lines.append(line)
        
        return value_lines + zero_lines

    def organize_matrix(self):
        line_found = False
        modified_lines = copy.deepcopy(self.lines)
        matrix_lines = []  
        j = 0
        for i, line in enumerate(self.lines):
            this_line = line.rstrip()
            if 'Sigind ' in line:
                line_found = True
                continue
            if line_found == True: 
                matrix_lines.append(this_line)
                if len(matrix_lines) == self.old_dimensions * self.layers:
                    line_found = False
                    break

        organized_matrix_lines = self.process_lines(matrix_lines)

        for i, line in enumerate(self.lines):
            this_line = line.rstrip()
            if 'Sigind ' in line:
                line_found = True
                continue
            if line_found == True: 
                modified_lines[i] = organized_matrix_lines[j] + '\n'
                j += 1
                if j == self.old_dimensions * self.layers:
                    break
        self.modified_content = "".join(modified_lines)
        self.lines = modified_lines

    def transformation_matrix(self):
        if self.layers == 2:
            repetition = [2]
        elif self.layers == 3:
            repetition = [2, 1]
        elif self.layers == 4:
            repetition = [2, 2]
        elif self.layers == 5:
            repetition = [2, 2, 1]

        modified_lines = []
        start_index = self.lines.index('#---------------- # Transformation matrix follows -----------\n') +1
        modified_lines = self.lines[:start_index]
        end = self.lines[start_index+self.old_dimensions:]
        for i in range(start_index, start_index+ self.old_dimensions, self.nodes):
            block = []
            for line in self.lines[i:i+self.nodes]:
                if '\n' not in line:
                    line += '\n'
                block.append(line)
            new_block = self.adjust_lines(block, repetition)
            for row in new_block:
                modified_lines.append(row)
        for row in end:
            modified_lines.append(row)
        self.modified_content = "".join(modified_lines)
        self.lines = modified_lines 
                    
    def adjust_lines(self, block, repetition):
        new_block = []
        before = 0
        after = self.layers
        for i, value in enumerate(repetition):
            after -= value
            for a in range(value):
                k = 0
                for line in block:
                    new_line =  '0.00000000  0.00000000    '*int(self.old_dimensions/self.nodes)*((value-1)*k+a) + line[1:-1] +'    0.00000000  0.00000000'*int(self.old_dimensions/self.nodes)*((value-1)*(self.nodes-k)-a)
                    new_line =' ' + '0.00000000  0.00000000    '*self.old_dimensions*before + new_line + '    0.00000000  0.00000000'*self.old_dimensions*after + '\n'
                    new_block.append(new_line)
                    k += 1
            before += value
        return new_block

    def delete_matrices(self):
        cix_num_count = 0
        modified_lines = []
        
        i = 0
        while i < len(self.lines):
            line = self.lines[i]
            if 'cix-num' in line:
                cix_num_count += 1
                if cix_num_count > 1 and cix_num_count <= self.layers:
                    i = i + 5 + 2*self.old_dimensions 
                    continue
                elif cix_num_count > self.layers:
                    line = self.replace_number(line, str(cix_num_count), str(cix_num_count + 1 - self.layers))
                    
            modified_lines.append(line)
            i += 1
        self.modified_content = "".join(modified_lines)
        self.lines = modified_lines

    def write_modified_content(self, new_file_name):
        with open(new_file_name, 'w') as file:
            file.write(self.modified_content)

    def process_file(self):
        new_file_name = self.file_name[:-7] + '_modified.indmfl'
        self.validate_layers()
        self.read_file()
        self.replace_cix()
        self.replace_parameters()
        self.independent_components()
        self.delete_matrices()
        self.create_sigind_matrix()
        self.replace_sigind_matrix()
        self.organize_matrix()
        self.transformation_matrix()
        self.write_modified_content(new_file_name)
        print(f"The file '{self.file_name}' has been modified and saved as '{new_file_name}'.")

if __name__ == "__main__":
    file_name = 'dmft_U12_bz2.indmfl'
    # file_name = 'dmft_U12_bz2_3couches.indmfl'
    # file_name = 'dmft_U12_bz2_4couches.indmfl'
    # file_name = 'dmft_U12_bz2_5couches.indmfl'
    layers = 2
    nodes = 4


    modifier = DMFTFileModifier(file_name, layers, nodes)
    try:
        modifier.process_file()
    except Exception as e:
        print("Error:", e)
