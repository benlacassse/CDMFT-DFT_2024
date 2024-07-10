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
                if cix_count <= self.layers * self.nodes:
                    modified_lines[i] = line[:position] + ' ' * (len(matches[2]) - 1) + '1' + line[position + len(matches[2]):]
                elif cix_count <= self.layers * self.nodes * 2:
                    modified_lines[i] = line[:position] + ' ' * (len(matches[2]) - 1) + str(int(matches[2]) - self.layers + 1) + line[position + len(matches[2]):]
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
                    new_kcix_blocks = str(int(kcix_blocks) - self.layers + 1).ljust(len(kcix_blocks))
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
                line = (line.rstrip() + ' ') * self.layers + '\n'
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
                modified_number = number + add_value if number != 0 else part
                modified_line.append(str(modified_number))
            else:
                raise ValueError("The line contains invalid characters.")
        return ' '.join(modified_line)

    def block_finder(self):
        search_line = False
        block_ind = 0
        block_list = [] 
        for i, line in enumerate(self.lines):
            if 'Sigind ' in line:
                search_line = True
                continue
            if search_line and block_ind < self.nodes:
                j = 0
                for char in line.strip():
                    j += 1
                    if char.isdigit() and char != '0':
                        line_block = line[:j]
                block_list.append(line_block)
                block_ind += 1
        return block_list

    def create_sigind_matrix(self):
        matrix_data = []
        sigind_found = False
        j = 0
        for i, line in enumerate(self.lines):
            if 'Sigind ' in line:
                sigind_found = True
                continue
            if sigind_found and j < self.old_dimensions:
                row = list(map(int, line.split()))
                matrix_data.append(row)
                j += 1
        matrix_array = np.array(matrix_data)
        maximum = np.max(matrix_array)
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
            if 'Sigind ' in line:
                sigind_count += 1
                j = 0
                continue
            if sigind_count == 1 and j < self.old_dimensions:
                if line == '0 ' * (self.old_dimensions - 1) + '0\n':
                    for k in range(self.layers):
                        new_line = line.rstrip() + ' '
                        if k == 0:
                            modified_lines[i] = new_line * self.layers + '\n'
                        else:
                            modified_lines.insert(i + self.old_dimensions, new_line * self.layers + '\n')
                else:
                    for k in range(self.layers):
                        new_line = ''
                        for l in range(self.layers):
                            for m in range(self.nodes):
                                new_line += str(self.create_sigind_matrix()[k][l][j][m]) + ' '
                        if k == 0:
                            modified_lines[i] = new_line + '0 ' * (self.old_dimensions - self.nodes) * self.layers + '\n'
                        else:
                            modified_lines.insert(i + self.old_dimensions + (k - 1) * (j + 1), new_line + '0 ' * (self.old_dimensions - self.nodes) * self.layers + '\n')
                j += 1
            elif sigind_count > 1 and j < self.old_dimensions / self.nodes:
                modified_lines[i + self.old_dimensions * (self.layers - 1)] = self.modify_line(line.rstrip()) + '\n'
                j += 1
        self.modified_content = "".join(modified_lines)
        self.lines = modified_lines

    def generate_diagonal_pattern_matrix(self):
        patterns = {
            2: ([1, 2], 1),
            3: ([1, 4, 3], 2),
            4: ([1, 6, 4, 5], 4),
            5: ([1, 9, 7, 8, 6], 6)
        }
        n = self.layers
        if n not in patterns:
            raise ValueError("Unsupported number of layers.")
        pattern, jump = patterns[n]
        matrix = np.array([[pattern[(i + j) % n] for j in range(n)] for i in range(n)]) - 1
        return matrix, jump

    def process_lines(self, lines):
        value_lines = [line for line in lines if not all(ch == '0' for ch in line.replace(' ', ''))]
        zero_lines = [line for line in lines if all(ch == '0' for ch in line.replace(' ', ''))]
        return value_lines + zero_lines

    def organize_matrix(self):
        line_found = False
        modified_lines = copy.deepcopy(self.lines)
        matrix_lines = []
        j = 0
        for line in self.lines:
            if line_found:
                if j < self.old_dimensions:
                    matrix_lines.append(line)
                    j += 1
                    if j == self.old_dimensions:
                        processed_lines = self.process_lines(matrix_lines)
                        modified_lines[i - self.old_dimensions + 1: i + 1] = processed_lines
                        break
            if ' '.join(line.split()).isdigit():
                line_found = True
                i = self.lines.index(line)
        self.lines = modified_lines
        self.modified_content = "".join(modified_lines)

    def write_modified_file(self):
        if not self.modified_content:
            raise ValueError("No modifications were made to the file.")
        with open(f"new_{self.file_name}", 'w') as file:
            file.write(self.modified_content)

def modify_dmft_file(file_name, layers, nodes):
    modifier = DMFTFileModifier(file_name, layers, nodes)
    modifier.read_file()
    modifier.validate_layers()
    modifier.replace_cix()
    modifier.replace_parameters()
    modifier.independent_components()
    modifier.replace_sigind_matrix()
    modifier.organize_matrix()
    modifier.write_modified_file()

modify_dmft_file("dmft_U12_bz2.indmfl", 2, 4)
