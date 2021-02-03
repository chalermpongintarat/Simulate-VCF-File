#!/usr/local/bin/python
# Build by Chalermpong Intarat

import os, sys, requests, time, progressbar, Bio

from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein

# Module: Answers for test cases with missing positions
button = ""

def main():
    while True:
        print("")
        print("======================================================================================")
        print(" Answers for test cases with missing positions")
        print("======================================================================================")
        print(" Chalermpong Intarat")
        print(" NBT, NSTDA")
        print(" chalermpong.int@biotec.or.th")
        print(" Apr 13, 2020")
        print("======================================================================================")
        print("")
        print("Main Menu: Cases")
        print("")
        print("\tModule 1. Case: Phased genotypes with missing positions\t\t: Button[1]")
        print("\tModule 2. Case: Unphased genotypes with missing positions\t: Button[2]")
        print("\tModule 3. Exit. \t\t\t\t\t\t: Button[3]")
        print("")
        try:
            button = int(input("\t\tPlease push your Button. \t\t\t\t: "))
        except ValueError:
            print("")
            print("\t\t\t\t\t\t\t\t\t: No Module.")
            print("")
            continue
        else:
            if button == 1:
                print("")
                import phased_genotypes_missing_positions
                phased_genotypes_missing_positions.main_phased_genotypes_missing_positions()
            elif button == 2:
                print("")
                print("\t\t\t\t\t\t\t\t\t: Coming soon!")
                print("")
            elif button == 3:
                print("")
                bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
                for i in range(100):
                    time.sleep(0.1)
                    bar.update(i)
                print(" Exit.")
                print("")
                sys.exit(0)
            else:
                print("")
                print("\t\t\t\t\t\t\t\t\t: No Module.")
                print("")

if __name__ == "__main__":
    main()
