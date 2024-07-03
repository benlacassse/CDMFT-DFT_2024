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
        # Check if the line is only made of zeros or spaces
        if set(line.strip()) <= {'0', ' '}:
            return line
        
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

    def transformation_matrix(self):
        search_line = False
        modified_lines = copy.deepcopy(self.lines)  
        j=1
        for i, line in enumerate(self.lines):
            if search_line:
                if j <= self.old_dimensions:
                    this_line = line.rstrip()[1:]
                    if line.rstrip()[0] == ' ':
                        new_line = ' '
                    else:
                        new_line = '-'
                    for l in range(self.layers):
                        for k in range(self.layers):
                            if l == k:
                                new_line = new_line + this_line + '    '
                            else:
                                new_line = new_line + '0.00000000  0.00000000    '*self.old_dimensions
                            
                        if l == 0:
                            modified_lines[i] = new_line + '\n'
                        else:
                            modified_lines.insert(i + self.old_dimensions+(l-1)*j, new_line + '\n')
                        new_line = ' '
                    j += 1
                else:
                    search_line = False

            if line.strip().startswith('#---------------- # T'):
                search_line = True

        self.modified_content = "".join(modified_lines)
        self.lines = modified_lines
    
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
        self.replace_parameters()
        self.independent_components()
        self.delete_matrices()
        self.sigind_matrix()
        self.transformation_matrix()
        self.write_modified_content(new_file_name)
        print(f"The file '{self.file_name}' has been modified and saved as '{new_file_name}'.")


if __name__ == "__main__":
    # Define the file name and parameters
    file_name = 'dmft_U12_bz2.indmfl'
    new_file_name = 'dmft_U12_bz2_modified.indmfl'
    layers = 2
    nodes = 4


    # Create an instance of the DMFTFileModifier class and process the file
    modifier = DMFTFileModifier(file_name, layers, nodes)
    try:
        modifier.process_file(new_file_name)
    except Exception as e:
        print("Error:", e)
