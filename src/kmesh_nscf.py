#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys

def main():
    """
    Main function for generating k-points.

    Usage:
        n1 n2 n3 [wan]
        n1  - divisions along 1st recip vector
        n2  - divisions along 2nd recip vector
        n3  - divisions along 3rd recip vector
        wan - omit the k-point weight (optional)
    """
    numargs = len(sys.argv) - 1  # Subtract 1 to exclude the script name
    if numargs < 3 or numargs > 4:
        print("usage: n1 n2 n3 [wan]")
        print("       n1  - divisions along 1st recip vector")
        print("       n2  - divisions along 2nd recip vector")
        print("       n3  - divisions along 3rd recip vector")
        print("       wan - omit the kpoint weight (optional)")
        sys.exit(1)
    n1 = int(sys.argv[1])
    n2 = int(sys.argv[2])
    n3 = int(sys.argv[3])
    if n1 <= 0:
        print("n1 must be > 0")
        sys.exit(1)
    if n2 <= 0:
        print("n2 must be > 0")
        sys.exit(1)
    if n3 <= 0:
        print("n3 must be > 0")
        sys.exit(1)
    totpts = n1 * n2 * n3
    if numargs == 3:
        print("K_POINTS crystal")
        print(totpts)
        for x in range(n1):
            for y in range(n2):
                for z in range(n3):
                    print(f"{x/n1:.8f} {y/n2:.8f} {z/n3:.8f} {1/totpts:.6e}")
    if numargs == 4:
        for x in range(n1):
            for y in range(n2):
                for z in range(n3):
                    print(f"{x/n1:.8f} {y/n2:.8f} {z/n3:.8f}")
    sys.exit(0)
if __name__ == "__main__":
    main()
