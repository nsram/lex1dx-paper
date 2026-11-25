"""12 values, k = 3: prints the result in the same format as examples/cpp/twelve.cpp and the
Scala examples.Twelve, so the three can be compared line by line (bits are exact)."""
import struct

import lex

x = [5.5, 1.2, 9.8, 2.0, 1.0, 5.0, 10.4, 1.9, 6.1, 1.2, 2.1, 5.5]
r = lex.lex_hybrid(x, 3)
bits = lambda v: struct.unpack("<Q", struct.pack("<d", v))[0]
print("sizes  ", " ".join(str(s) for s in r.sizes))
print("centers", " ".join(f"{c:.6f}" for c in r.centers))
print(f"sse     {r.sse:.6f}  bits {bits(r.sse):016x}")
print("centers bits", " ".join(f"{bits(c):016x}" for c in r.centers))
