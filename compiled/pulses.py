#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE, call
import pprint

def is_good(runs,tag):
    for r in runs:
        c = r + '.root'
        #print(" look for ",c," in ", tag)
        if c in tag:
            #print(' tag ', tag, ' for ', r)
            return r 
    return 0

def main(args):
    """ run pulses """
    myEnv = os.environ.copy()
    frunlist = open("run1GoodRuns.txt", "r")
    runs = []
    for x in frunlist:
        r = x[0:x.rindex("\n")]
        runs.append(r)

    files = []
    p = os.listdir('rootData')
    print(" number of files in rootData ",len(p))
    for i in p:
        if os.path.isfile('rootData/'+i):
            #print(" file ", i)
            if( i.endswith("root")  and not i.startswith("ana") ) :
                tag = i[0:i.rindex(".")]
                check = is_good(runs,i)
                if( int(check) > 0):
                    files.append(tag)
                    print("\t  good run ", len(files), " check ", check, " tag " , tag )



    #print(myEnv)
    
    n = len(runs)
    if (n < 1):
        print("\n no runs found ")
        return


    print(" all good runs",  len(runs) , " good files found ", len(files) )
    n = len(files) 
    if (len(sys.argv) > 1):
        n = int(args[0])

    #print(runs)
    print(" args ", args, " number of runs to run   ", n)

    exit();

    for i in range(0, n):
        print(" run job ", i, " tag ",runs[i])
        process = Popen(['pulses', runs[i]], stdout=PIPE,stderr=PIPE, env=myEnv)
        stdout, stderr = process.communicate()
        process.wait()
        print(stdout)
        print(stderr)


if __name__ == '__main__':
    main(sys.argv[1:])
