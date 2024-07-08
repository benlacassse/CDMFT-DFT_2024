def process_line(line):
    # Remove spaces and process the digits
    digits = line.replace(' ', '')
    non_zero_chars = [ch for ch in digits if ch != '0']
    zero_count = digits.count('0')
    
    # Append zeros to the end and re-add spaces between digits
    processed_digits = non_zero_chars + ['0'] * zero_count
    return ' '.join(processed_digits)

def process_lines(lines):
    processed_lines = []
    zero_only_lines = []
    
    for line in lines:
        stripped_line = line.replace(' ', '')
        if all(ch == '0' for ch in stripped_line):
            zero_only_lines.append(line)
        else:
            processed_lines.append(process_line(line))
    
    # Append lines with only zeros at the end
    return processed_lines + zero_only_lines

# Example usage:
lines = [
    "1 0 0 2 0 0 3",
    "0 0 0 0 0 0 0",
    "1 2 0 3 0 5 0",
    "0 0 0 0 0 0 0",
    "4 5 6 7 0 0 0",
    "7 0 0 3 0 0 2"
]

processed_lines = process_lines(lines)

for line in processed_lines:
    print(line)