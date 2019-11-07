#!/usr/bin/python3

# Path of executable
DiTwi_bin = '/home/akva/kode/IFEM/Apps/IFEM-DiTwi/r/bin/DiTwi'
N = 31
MAX = 1e-4

current_target = 3.3402660055362014e-05
ifem_pipe = 0

import math
import os
import pickle
import random
import re
import subprocess
import time

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
  os.write(ifem_pipe, b'<callbacks><ditwi><new_target>%f</new_target></ditwi></callbacks>' %(x))
  value = ReadValue()
  print(value)
  return value

# Start ditwi
args = [DiTwi_bin, 'instance.xinp', '-controller', '-vtf', '1']
with subprocess.Popen(args, stdout=subprocess.PIPE) as proc:
    time.sleep(2) # Wait for FIFO to appear
    ifem_pipe = os.open('ifem-control', os.O_NONBLOCK | os.O_WRONLY)
    for i in range(1,N):
        # Draw random number
        current_target = MAX * math.sin(2*math.pi/N*i)

        print('Look for s_zz = %f' %(current_target))
        load = ObjectiveFunction(current_target)
        os.write(ifem_pipe, b'<callbacks><ditwi><save_step/><step_ok/></ditwi></callbacks>')
    os.write(ifem_pipe, b'<callbacks><ditwi><quit/></ditwi></callbacks>')
    proc.wait()
