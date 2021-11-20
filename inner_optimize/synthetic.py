#!/usr/bin/python3

import math
import os
import pickle
import random
import re
import subprocess
import time

# Path of executable
if 'DITWI_EXE' in os.environ:
  DiTwi_bin = os.environ['DITWI_EXE']
else:
  DiTwi_bin = '/home/akva/kode/IFEM/Apps/IFEM-DiTwi/r/bin/DiTwi'

N = 31
MAX = 1e-4

ifem_pipe = 0

def get_lines_until(func, sentinel):
  while True:
    line = func()
    if line == sentinel:
      raise StopIteration
    yield line

def ReadValue():
  output = list(get_lines_until(proc.stdout.readline, b'--- end of iteration ---\n'))
  print(output)
  value = 0.0
  for line in output:
     ma = re.match('.*sol2 = ', str(line))
     if not ma is None:
       arr = re.split(r' +', str(line))
       value = float(arr[7])
       break
  return value


def ObjectiveFunction(x):
  print('Look for s_zz = %f' % (x))
  os.write(ifem_pipe, b'<callbacks><ditwi><new_target>%f</new_target></ditwi></callbacks>' %(x))
  value = ReadValue()
  print('Got %f' % (value))
  return value

# Start ditwi
args = [DiTwi_bin, 'instance.xinp', '-controller']
with subprocess.Popen(args, stdout=subprocess.PIPE) as proc:
    time.sleep(2) # Wait for FIFO to appear
    ifem_pipe = os.open('ifem-control', os.O_NONBLOCK | os.O_WRONLY)
    for i in range(1,N):
        # Draw random number
        current_target = MAX * math.sin(2*math.pi/N*i)
        load = ObjectiveFunction(current_target)
        os.write(ifem_pipe, b'<callbacks><ditwi><save_step/><step_ok/></ditwi></callbacks>')
    os.write(ifem_pipe, b'<callbacks><ditwi><quit/></ditwi></callbacks>')
    proc.wait()
