import re
import os

def custom_sort(file):
    # Example custom comparison function: sort by length of file name, then alphabetically
    matches = re.findall(r'\d+', file)
    number = int(matches[1]) if len(matches) > 1 else float('inf')
    return number

def calculate_average_for_files(directory):
    # List all files in the specified directory
    files = sorted([f for f in os.listdir(directory) if f.startswith('dim')], key=custom_sort)
    
    if not files:
        print("No text files found in the directory.")
        return

    # Iterate through each file
    for file in files:
        file_path = os.path.join(directory, file)
        numbers = []

        # Read the file and extract numbers
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    try:
                        numbers.append(float(line.strip()))
                    except ValueError:
                        print(f"Skipping invalid number in file {file}: {line.strip()}")
        except Exception as e:
            print(f"Error reading file {file}: {e}")
            continue

        # Calculate and display the average
        if numbers:
            average = sum(numbers) / len(numbers)
            #print(f"Average for {file}: {average:.2f}")
            print(file.split('=')[-1], average)
        else:
            print(f"No valid numbers found in {file}.")

if __name__ == "__main__":
    directory = os.getcwd()  # Use the current directory by default
    calculate_average_for_files(directory)
