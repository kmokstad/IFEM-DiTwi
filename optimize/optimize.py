#!/usr/bin/python3

# Path of executable
DiTwi_bin = '/home/akva/kode/IFEM/Apps/IFEM-DiTwi/r/bin/DiTwi'

current_target = 5e-5
ifem_pipe = 0

import os
import re
import scipy.optimize
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
  value = 0.0
  for line in output:
     ma = re.match('.*sol2 = ', str(line))
     if not ma is None:
       arr = re.split(r' +', str(line))
       value = float(arr[3])
       break
  return value


def ObjectiveFunction(x):
  print(x)
  os.write(ifem_pipe, b'<callbacks><ditwi><new_load>%f</new_load></ditwi></callbacks>' %(x*1e6))
  value = ReadValue()
  print(value)
  return value - current_target

# Start ditwi
args = [DiTwi_bin, 'instance.xinp', '-controller']
with subprocess.Popen(args, stdout=subprocess.PIPE) as proc:
    time.sleep(1) # Wait for FIFO to appear
    ifem_pipe = os.open('ifem-control', os.O_NONBLOCK | os.O_WRONLY)
    res = scipy.optimize.minimize(ObjectiveFunction, 1, tol=1e-10, options={'eps':0.1})
    print(res)
    os.write(ifem_pipe, b'<callbacks><ditwi><step_ok/></ditwi></callbacks>')
    current_target = 1e-4
    res = scipy.optimize.minimize(ObjectiveFunction, res.x, tol=1e-10, options={'eps':0.1})
    print(res)
    os.write(ifem_pipe, b'<callbacks><ditwi><step_ok/></ditwi></callbacks>')
    current_target = 1e-3
    res = scipy.optimize.minimize(ObjectiveFunction, res.x, tol=1e-10, options={'eps':0.1})
    print(res)
    proc.kill()
