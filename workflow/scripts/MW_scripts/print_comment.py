#!/usr/bin/env python3
# This script prints any comment in a structured and pretty way.
import sys

if len(sys.argv) < 3:
    print("Usage: script <comment> <delimiter>")
    sys.exit(1)

comm = sys.argv[1]
delim = sys.argv[2]

print('\n' + delim * 120)

max_len = 90

cut = comm.split(" ")
line = ""
for word in cut:
    if (len(line) + 1 + len(word)) > max_len:
        edge1 = (120 - len(line)) // 2 - 5
        edge2 = 120 - edge1 - len(line) - 10
        print(delim * 5 + " " * edge1 + line + " " * edge2 + delim * 5)
        line = word
    else:
        line = line + " " + word
edge1 = (120 - len(line)) // 2 - 5
edge2 = 120 - edge1 - len(line) - 10
print(delim * 5 + " " * edge1 + line + " " * edge2 + delim * 5)

print(delim * 120 + '\n')
