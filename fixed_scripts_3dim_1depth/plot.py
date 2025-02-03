import re
import os
import sys
import matplotlib.pyplot as plt
import statistics

xs=[]
y_sv=[]
y_zer=[]

fig, axes = plt.subplots(5, 2, figsize=(10, 15))

# Flatten axes array for easier iteration
axes = axes.flatten()
svs=[[],[],[],[],[],[],[],[],[],[]]
zers=[[],[],[],[],[],[],[],[],[],[]]

def custom_sort(file):
    # Example custom comparison function: sort by length of file name, then alphabetically
    matches = re.findall(r'\d+', file)
    number = int(matches[1]) if len(matches) > 1 else float('inf')
    return number

def calculate_average_for_files(directory, which, dim):
    # List all files in the specified directory
    if which != 'zero' and which != 'sv':
        print("Error which value")
        return 0

    files = sorted([f for f in os.listdir(directory) if f.startswith(which+'_dim='+str(dim))], key=custom_sort)
    
    if not files:
        print("No text files found in the directory.")
        return

    # Iterate through each file
    for file in files:
        file_path = os.path.join(directory, file)
        numbers = []

        i=0
        # Read the file and extract numbers
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    try:
                        numbers.append(float(line.strip()))
                        if which == 'sv':
                            svs[i].append(numbers[-1])
                        else:
                            zers[i].append(numbers[-1])
                        i+=1
                    except ValueError:
                        print(f"Skipping invalid number in file {file}: {line.strip()}")
        except Exception as e:
            print(f"Error reading file {file}: {e}")
            continue

        # Calculate and display the average
        if numbers:
            average = sum(numbers) / len(numbers)
            #average=statistics.median(numbers)
            num_instances = len(numbers)
            #print(num_instances)
            #print(f"Average for {file}: {average:.2f}")
            print(file.split('=')[-1], average)
            if which == 'sv':
                xs.append(int(file.split('=')[-1]))
                y_sv.append(average)
            else:
                y_zer.append(average)

        else:
            print(f"No valid numbers found in {file}.")

if __name__ == "__main__":
    directory = os.getcwd()  # Use the current directory by default
    dim = sys.argv[1]
    print('zero')
    calculate_average_for_files(directory, 'zero', dim)
    print('sv')
    calculate_average_for_files(directory, 'sv', dim)
    
    print(xs)
    plt.plot(xs, y_sv, color='red')
    plt.plot(xs, y_zer, color='blue')
    plt.legend(["SV", "ZERO"])

    plt.show()
    kok
    
    """
    print(svs)
    num_instances=sum(1 for s in svs if s)
    print("Num instances:", num_instances)
    for i in range(num_instances):
        #print(xs, svs[i])
        print(xs[:-1])
        axes[i].plot(xs[:len(svs[i])], svs[i], color='red')
        axes[i].plot(xs[:len(svs[i])], zers[i], color='blue')
    plt.show()
    """
