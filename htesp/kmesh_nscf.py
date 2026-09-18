#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to generate k-points"""
import sys

def main(argv=None):
    """
    Main function for generating k-points.

    FIX(21): ``main`` is a callable top-level entry point; it accepts an
    explicit argument list and returns an exit code instead of always calling
    ``sys.exit()``, so importing and calling it does not kill the caller.

    Usage:
        n1 n2 n3 [wan]
        n1  - divisions along 1st recip vector
        n2  - divisions along 2nd recip vector
        n3  - divisions along 3rd recip vector
        wan - omit the k-point weight (optional)
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    numargs = len(argv)
    if numargs < 3 or numargs > 4:
        print("usage: n1 n2 n3 [wan]")
        print("       n1  - divisions along 1st recip vector")
        print("       n2  - divisions along 2nd recip vector")
        print("       n3  - divisions along 3rd recip vector")
        print("       wan - omit the kpoint weight (optional)")
        return 1
    n1 = int(argv[0])
    n2 = int(argv[1])
    n3 = int(argv[2])
    if n1 <= 0:
        print("n1 must be > 0")
        return 1
    if n2 <= 0:
        print("n2 must be > 0")
        return 1
    if n3 <= 0:
        print("n3 must be > 0")
        return 1
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
    return 0
if __name__ == "__main__":
    sys.exit(main())
