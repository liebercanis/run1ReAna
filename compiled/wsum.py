#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE, call
import pprint


def main(args):
    """ run convertRaw """
    myEnv = os.environ.copy()
    #print(myEnv)
    files = []
    p = os.listdir('rootData/DS2')
    print(" number of files in rootData ",len(p))
    for i in p:
        if os.path.isfile('rootData/DS2/'+i):
            #print(" file ", i)
            if( i.endswith("root")  and i.startswith("RecirculateRun") ) :
                tag = i[0:i.rindex(".")]
                files.append(tag)
                #print("\n\t data root file ", i," tag ",tag)


    n = len(files)
    if (n < 1):
        print("\n no files found ")
        return


    print(" all files  ", len(p), " will run  ", len(files))
    if (len(sys.argv) > 1):
        n = int(args[0])

    #print(files)
    print(" args ", args, " number of files to run   ", n)

    for i in range(0, n):
        print(" run job ", i, " tag ",files[i])
        process = Popen(['sum', files[i]], stdout=PIPE,stderr=PIPE, env=myEnv)
        stdout, stderr = process.communicate()
        process.wait()
        print(stdout)
        print(stderr)


if __name__ == '__main__':
    main(sys.argv[1:])
