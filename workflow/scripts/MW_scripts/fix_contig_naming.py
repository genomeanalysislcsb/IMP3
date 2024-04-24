#!/usr/bin/env python3
import sys

# Open the file specified by the first command line argument
with open(sys.argv[1]) as file:
    for line in file:
        if line[0] == ">":
            # Replace "=" with "_" in the line
            modified_line = line.replace("=", "_")
            # Print the modified line without adding a new line character at the end
            sys.stdout.write(modified_line)
        else:
            # Print the line without its trailing newline and without adding a new line at the end
            # .rstrip() removes trailing whitespace including newline characters
            print(line.rstrip())
