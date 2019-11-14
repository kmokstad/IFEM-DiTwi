#!/usr/bin/python3

import os
import pickle
import random
import re
import scipy.optimize
import subprocess
import time

# Path of executable
if 'DITWI_EXE' in os.environ:
  DiTwi_bin = os.environ['DITWI_EXE']
else:
  DiTwi_bin = '/home/akva/kode/IFEM/Apps/IFEM-DiTwi/r/bin/DiTwi'

ifem_pipe = 0

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
       value = float(arr[7])
       break
  return value


def ObjectiveFunction(x):
  os.write(ifem_pipe, b'<callbacks><ditwi><new_target>%f</new_target></ditwi></callbacks>' %(x))
  value = ReadValue()
  print(value)
  return value

# Load dictionary
try:
  p_in = open("cases.pickle", 'rb')
  case_dict = pickle.load(p_in)
except:
  pass
  case_dict = {}

print(case_dict)

# Start ditwi
args = [DiTwi_bin, 'instance.xinp', '-controller']
with subprocess.Popen(args, stdout=subprocess.PIPE) as proc:
    time.sleep(2) # Wait for FIFO to appear
    ifem_pipe = os.open('ifem-control', os.O_NONBLOCK | os.O_WRONLY)
    random.seed()
    start = time.time()
    reqs = 0
    while True:
        if time.time() - start > 1:
            print('Processed %i requests at %f rps' %(reqs, reqs/(time.time()-start)))
            start = time.time()

        dict_changed = True
        # Draw random number
        current_target = random.uniform(0.0, 1e-4) + 1e-5
        reqs = reqs + 1

        # Find close enough key in dictionary?
        for key in case_dict.keys():
            if abs(current_target - key) < 5e-7:
#                 print('Value from cache: %f -> %f' %(key, case_dict[key]))
                dict_changed = False
                break

        if dict_changed:
            print('Look for s_zz = %f' %(current_target))
            load = ObjectiveFunction(current_target)
#             res = scipy.optimize.minimize(ObjectiveFunction, 2.0, tol=3e-5, method='Nelder-Mead')
            #print(res)
#            if res.success:
            case_dict[current_target] = load # res.x[0]
            p_out = open("cases.pickle", "wb")
            pickle.dump(case_dict, p_out)
            p_out.close()
#                 os.write(ifem_pipe, b'<callbacks><ditwi><step_ok/></ditwi></callbacks>')
    proc.kill()
