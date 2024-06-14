import re

class DMFTFileModifier:
    def __init__(self, file_name, layers=1, nodes=4):
        self.file_name = file_name
        self.layers = layers
        self.nodes = nodes
        self.lines = []
        self.content = ""
        self.modified_content = ""

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

    def find_and_replace_number(self):
        search_next_line = False
        for line in self.lines:
            if search_next_line:
                match = re.search(r'\d+', line)
                if match:
                    number_blocks = match.group(0)
                    if int(number_blocks) != (self.nodes + 1) * self.layers:
                        raise ValueError('Input of layers and nodes do not match the number of independent kcix blocks.')
                    else:
                        new_number = str(int(number_blocks) - layers + 1)
                        if len(new_number) < len(number_blocks):
                            new_number = new_number.ljust(len(number_blocks))
                        self.modified_content = self.content.replace(number_blocks, new_number, 1)
                else:
                    raise ValueError("No number found on the line following '#'.")
                break
            if line.strip().startswith('#'):
                search_next_line = True

        if not search_next_line:
            raise ValueError("No line starting with '#' found in the file.")

    def write_modified_content(self):
        with open(self.file_name, 'w') as file:
            file.write(self.modified_content)


    def process_file(self):
        self.validate_layers()
        self.read_file()
        self.find_and_replace_number()
        self.write_modified_content()
        print(f"The first number in the file '{self.file_name}' has been replaced.")


if __name__ == "__main__":
    # Define the file name and parameters
    file_name = 'dmft_U12_bz2.indmfl'
    layers = 2
    nodes = 4

    # Create an instance of the DMFTFileModifier class and process the file
    modifier = DMFTFileModifier(file_name, layers, nodes)
    try:
        modifier.process_file()
    except Exception as e:
        print("Error:", e)
