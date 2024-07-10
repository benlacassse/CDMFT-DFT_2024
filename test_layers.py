    # def sigind_matrix(self):
    #     sigind_count = 0
    #     modified_lines = copy.deepcopy(self.lines)  
    #     liste = list(np.hstack((np.arange(self.layers)[::-1], np.arange(1, self.layers))))
    #     for i, line in enumerate(self.lines):
    #         this_line = line.rstrip()
    #         if 'Sigind ' in line:
    #             sigind_count += 1
    #             j=1
    #             continue
    #         if sigind_count == 1 : 
    #             if j <= self.old_dimensions:
    #                 for l in range(self.layers):
    #                     new_line = ''
    #                     for k in range(self.layers):
    #                         new_line = new_line + self.modify_line(this_line, liste[self.layers-1+k-l]) + ' '
                            
    #                     if l == 0:
    #                         modified_lines[i] = new_line + '\n'
    #                     else:
    #                         modified_lines.insert(i + self.old_dimensions+(l-1)*j, new_line + '\n')
    #                 j += 1
    #         elif sigind_count > 1:
    #             if j <= self.old_dimensions/self.nodes:
    #                 this_line = line.rstrip()
    #                 modified_lines[i + self.old_dimensions*(self.layers-1)] = self.modify_line(this_line, self.layers-1) +'\n'
    #                 j += 1


    #     self.modified_content = "".join(modified_lines)
    #     self.lines = modified_lines

    # def generate_diagonal_pattern_matrix(self):
    #     n = self.layers
    #     if n == 2:
    #         # liste = [0, 1]
    #         matrix = np.array([[1, 2],[2, 1]])
    #         jump = 1
    #     elif n == 3:
    #         # liste = [0, 2, 3]
    #         matrix = np.array([[1, 4, 3],[4, 1, 3], [3, 3, 2]])
    #         jump = 2
    #     elif n == 4:
    #         # liste = [0, 2, 4, 5]
    #         matrix = np.array([[1, 6, 4, 5],[6, 1, 5, 4], [4, 5, 2, 3], [5, 4, 3, 2]])
    #         jump = 4
    #     elif n == 5:
    #         # liste = [0, 3, 5, 7, 8]
    #         matrix = np.array([[1, 9, 7, 8, 6],[9, 1, 8, 7, 6], [7, 8, 2, 5, 4], [8, 7, 5, 2, 4], [6, 6, 4, 4, 3]])
    #         jump = 6

    #     for i, row in enumerate(matrix):
    #         for j, coefficient in enumerate(row):
    #             matrix[i][j] -= 1

    #     # # Initialize an n x n matrix with zeros
    #     # matrix = np.zeros((n, n), dtype=int)

    #     # for i in range(n):
    #     #     for j in range(n):
    #     #         diff = abs(i - j)
    #     #         # Determine the value for matrix[i][j]
    #     #         diag_num = min(i, j, n - 1 - i, n - 1 - j) + liste[diff]
    #     #         matrix[i, j] = diag_num
