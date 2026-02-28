#!/usr/bin/env python3
import sys
import math

# format asteptat: "<key> <succ> <N>"
exp_line = sys.argv[1].strip().split()
ekey = int(exp_line[0])
esucc = int(exp_line[1])
N = int(exp_line[2])

# linie de output
out_line = sys.argv[2].strip()

# Exemplu format asteptat de output:
# "Lookup 7: 1 -> 5 -> 10"

if not out_line.startswith(f"Lookup {ekey}:"):
    print(f"FAIL: Cheie gresita in output, se astepta Lookup {ekey}")
    sys.exit(2)

# Extrage ruta
route_str = out_line.split(":")[1].strip()
nodes = [int(x.strip()) for x in route_str.split("->")]

# Verifica succesor final
if nodes[-1] != esucc:
    print(f"FAIL: Succesor gresit. S-a primit {nodes[-1]}, se astepta {esucc}")
    sys.exit(2)

# Prag de complexitate
H_max = math.ceil(math.log2(N)) + 2
hop_count = len(nodes)

if hop_count <= H_max:
    print(f"PASS logN (hop_count={hop_count}, H_max={H_max})")
    sys.exit(0)
else:
    print(f"PASS_WITH_PENALTY O(N) (hop_count={hop_count}, H_max={H_max})")
    sys.exit(3)
