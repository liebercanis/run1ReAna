#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE, call
import pprint


def main(args):
    """ run pulses """
    myEnv = os.environ.copy()
    flist = open("run1GoodRuns.txt", "r")
    files = []
    for x in flist:
        tag = x[0:x.rindex("\n")]
        files.append(tag)



    #print(myEnv)
    
    n = len(files)
    if (n < 1):
        print("\n no files found ")
        return


    print(" all good files",  len(files))
    n = len(files) 
    if (len(sys.argv) > 1):
        n = int(args[0])

    #print(files)
    print(" args ", args, " number of files to run   ", n)

    for i in range(0, n):
        print(" run job ", i, " tag ",files[i])
        process = Popen(['pulses', files[i]], stdout=PIPE,stderr=PIPE, env=myEnv)
        stdout, stderr = process.communicate()
        process.wait()
        print(stdout)
        print(stderr)


if __name__ == '__main__':
    main(sys.argv[1:])
