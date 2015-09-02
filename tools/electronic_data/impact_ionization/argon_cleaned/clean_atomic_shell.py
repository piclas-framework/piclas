#---------------------------------------------------------------------------------------
# This script is written by Philip Ortwein Anno Domini 2015-08-19
#
#!/usr/bin/env python
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# import modules
import csv, math, sys
import commands
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
input_filename = 'Ar.csv'
small_filename = 'Database_reduced.csv'
#---------------------------------------------------------------------------------------

def BubbleSort( nsize,gold, told, gnew, tnew):
    for ii in range(0,nsize-1):
        print ii


print '\n Reading from file: ', input_filename, '\n'
input_file   = open(input_filename, 'r')
input_reader = csv.reader(input_file, delimiter=',')

nlines = 0
gOld = list()
Told = list()
for row in input_reader:
  #data.append( [ float(row[0]), float(row[1])] )
  gOld.append( [ float(row[0])] )
  Told.append( [ float(row[1])] )
  #print row[0], gOld[nlines], Told[nlines]
  nlines = nlines + 1
input_file.close()

print ' Found number of entries: ', nlines, '\n'

print ' ... Readin done.'

print ' Sorting entries...'
gNew=list()
TNew=list()
BubbleSort( nlines,gOld,Told,gNew,TNew)



# reduced_file = open(small_filename, 'w')
# for item in data:
#   reduced_file.write(" ".join(map(str, item[0:])) + "\n")
# 
# reduced_file.close()
