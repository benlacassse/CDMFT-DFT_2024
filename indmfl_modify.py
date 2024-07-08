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
        # Use regular expressions to ensure exact match
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
            
    def modify_line(self, line, index):
        modified_line = []
        for char in line.strip():
            if char.isdigit():
                if char != '0':
                    modified_line.append(str(int(char) + self.ind_components*index))
                else:
                    modified_line.append(char)
            elif char == ' ':
                modified_line.append(char)
            else:
                raise ValueError("The line contains invalid characters.")

        return ''.join(modified_line)

    def block_finder(self):
        search_line = False
        block_ind = 0
        block_list = [] 
        for i, line in enumerate(self.lines):
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

    def sigind_matrix(self):
        sigind_count = 0
        modified_lines = copy.deepcopy(self.lines)  
        liste = list(np.hstack((np.arange(self.layers)[::-1], np.arange(1, self.layers))))
        for i, line in enumerate(self.lines):
            this_line = line.rstrip()
            if 'Sigind ' in line:
                sigind_count += 1
                j=1
                continue
            if sigind_count == 1 : 
                if j <= self.old_dimensions:
                    for l in range(self.layers):
                        new_line = ''
                        for k in range(self.layers):
                            new_line = new_line + self.modify_line(this_line, liste[self.layers-1+k-l]) + ' '
                            
                        if l == 0:
                            modified_lines[i] = new_line + '\n'
                        else:
                            modified_lines.insert(i + self.old_dimensions+(l-1)*j, new_line + '\n')
                    j += 1
            elif sigind_count > 1:
                if j <= self.old_dimensions/self.nodes:
                    this_line = line.rstrip()
                    modified_lines[i + self.old_dimensions*(self.layers-1)] = self.modify_line(this_line, self.layers-1) +'\n'
                    j += 1


        self.modified_content = "".join(modified_lines)
        self.lines = modified_lines

    def process_line(self, line):
        # Remove spaces and process the digits
        digits = line.replace(' ', '')
        non_zero_chars = [ch for ch in digits if ch != '0']
        zero_count = digits.count('0')
        
        # Append zeros to the end and re-add spaces between digits
        processed_digits = non_zero_chars + ['0'] * zero_count
        return ' '.join(processed_digits)

    def process_lines(self, lines):
        processed_lines = []
        zero_only_lines = []
        
        for line in lines:
            stripped_line = line.replace(' ', '')
            if all(ch == '0' for ch in stripped_line):
                zero_only_lines.append(line)
            else:
                processed_lines.append(self.process_line(line))
        
        # Append lines with only zeros at the end
        return processed_lines + zero_only_lines

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
        modified_lines = []
        start_index = self.lines.index('#---------------- # Transformation matrix follows -----------\n') +1
        modified_lines = self.lines[:start_index]
        end = self.lines[start_index+self.old_dimensions:]
        for i in range(start_index, start_index+ self.old_dimensions, self.nodes):
            block = []
            for line in self.lines[i:i+self.nodes]:
                if '\n' not in line:
                    line = line + '\n'
                block.append(line)
            for j in range(self.layers):
                new_block = self.adjust_lines(block, j)
                for row in new_block:
                    modified_lines.append(row)
        for row in end:
            modified_lines.append(row)
        self.modified_content = "".join(modified_lines)
        self.lines = modified_lines 
                    
    def adjust_lines(self, block, a):
        new_block = []
        k = 0
        for line in block:
            new_line = ' '+'0.00000000  0.00000000    '*int(self.old_dimensions/self.nodes)*((self.layers-1)*k+a) + line[1:-1] +'    0.00000000  0.00000000'*int(self.old_dimensions/self.nodes)*((self.layers-1)*(self.nodes-k)-a)+ '\n'
            new_block.append(new_line)
            k += 1
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

    def process_file(self, new_file_name):
        self.validate_layers()
        self.read_file()
        self.replace_cix()
        self.replace_parameters()
        self.independent_components()
        self.delete_matrices()
        self.sigind_matrix()
        self.organize_matrix()
        self.transformation_matrix()
        self.write_modified_content(new_file_name)
        print(f"The file '{self.file_name}' has been modified and saved as '{new_file_name}'.")


if __name__ == "__main__":
    # Define the file name and parameters
    file_name = 'dmft_U12_bz2.indmfl'
    # file_name = 'dmft_U12_bz2 - Copy.indmfl'
    new_file_name = 'dmft_U12_bz2_modified.indmfl'
    layers = 2
    nodes = 4


    # Create an instance of the DMFTFileModifier class and process the file
    modifier = DMFTFileModifier(file_name, layers, nodes)
    try:
        modifier.process_file(new_file_name)
    except Exception as e:
        print("Error:", e)
