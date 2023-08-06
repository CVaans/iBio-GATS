import sys

from Home import HomeWindow

if __name__ == "__main__":
    filename = str(sys.argv[1])
    print("File name is {0}".format(filename))
    if not filename:
        print("Warning", "Please select a valid .fasta file")

    with open(filename, 'r') as f:
        sequence_target = f.read()
        HomeWindow(sequence_target)
