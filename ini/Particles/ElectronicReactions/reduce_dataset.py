#---------------------------------------------------------------------------------------
# This small peace of script reduces a large data set into one smaller
# or
# deleting not required values
#
# This script is written by Philip Ortwein Anno Domini 2013-03-11
#
#!/usr/bin/env python
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# list values
value1 = 35
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# import modules
import csv, math, sys
import commands
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
large_filename = 'Database_Ttr_10000.csv'
small_filename = 'Database_Tvib_10000.csv'
#---------------------------------------------------------------------------------------


print '\n Reading large file...\n'
large_file   = open(large_filename, 'r')
large_reader = csv.reader(large_file, delimiter=',')

ii = 0
data = list()
for row in large_reader:
  if ii == 0:
    data.append( [ row[0] , row[value1] ] )
  else:
    data.append( [ float(row[0]), float(row[value1]) ] )
  #print row[0], data[ii]
  ii = ii + 1
large_file.close()

print ' ... done.'


reduced_file = open(small_filename, 'w')
for item in data:
  reduced_file.write(" ".join(map(str, item[0:])) + "\n")

reduced_file.close()
