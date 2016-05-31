#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
from extract_userblock import get_userblock, get_part
import os
import re

parser = argparse.ArgumentParser(description='rebuild code revision from h5-state-file')
parser.add_argument('dir',help='name of empty directory where the rebuild will take place')
parser.add_argument('state',help='h5-state file with userblock')

args = parser.parse_args()

if not os.path.exists(args.dir) :
   os.mkdir(args.dir) 
else :
   if os.listdir(args.dir) != []: # is dir empty?
      print os.listdir(args.dir)
      print "Rebuild directory is not empty => exit!"
      exit(1)


# get svn
userblock = get_userblock(args.state)
git_url = get_part(userblock, "GIT URL")
git_url = git_url.strip()

# checkout code
git_revs = get_part(userblock,"GIT REVISIONS")
print git_revs
gitrevision = git_revs.split("master:")[1].strip().split()[0]
print gitrevision

cmd = "git clone " + git_url + " " + args.dir
subprocess.call(cmd, shell=True)
cmd = "git checkout -b rebuildbranch " + gitrevision
os.chdir(args.dir)
subprocess.call(cmd, shell=True)

# apply gitdiff
git_diff = get_part(userblock, "GIT DIFF")
f = open("diff.patch", 'w')
f.write(git_diff)
f.close()
subprocess.call("patch -p1 < diff.patch", shell=True)



# configure
builddir = "build"
os.mkdir(builddir)
cmake = get_part(userblock, "CMAKE")
f = open(os.path.join(builddir, "config.cmake"), 'w')
f.write(cmake)
f.close()

os.chdir(builddir)
p = subprocess.Popen(["cmake", "-C", "config.cmake" ,"../"])
p.wait()

# make
p = subprocess.Popen(["make"])

# write ini file
os.chdir("..")
os.mkdir("ini")
ini = get_part(userblock, "INIFILE")
f = open(os.path.join("ini", "parameter.ini"), 'w')
f.write(ini)
f.close()

