from scipy import *
from pylab import *
import json
import sys
import re

def center_of_mass(coords,num_parts,num_dim,time_slice):
    cntr = [0]*num_dim
    for dim in range(num_dim):
        for part_indx in range(num_parts):
            print "time_slice = ", time_slice,"part_indx=",  part_indx, "dim = ", dim
            cntr[dim] += coords[time_slice][part_indx][dim]
        cntr[dim] /= num_parts
    return cntr


print "running ", sys.argv[0]
print "number of args = ", len(sys.argv)

filename = "coords_out.json"
num = 0
skip = 1


if False:
   if len(sys.argv) > 1: 
       num = int(sys.argv[1])
       if len(sys.argv) == 3:
	  filename = sys.argv[2]
    
commandline = ' '.join(sys.argv[1:])

matchobj = re.search(r"-file\w*\s(.+?)(\s|$)",commandline)
if (matchobj):
         filename = matchobj.group(1)

print "filename = ", filename

matchobj = re.search(r"-num\w*\s(.+?)(\s|$)",commandline)
if (matchobj):
         num = int(matchobj.group(1))

start =   0
matchobj = re.search(r"-star\w*\s(.+?)(\s|$)",commandline)
if (matchobj):
         start = int(matchobj.group(1))


matchobj = re.search(r"-skip\w*\s(.+?)(\s|$)",commandline)
if (matchobj):
         skip = int(matchobj.group(1))


if re.search(r"-help",commandline):
   print "plot_run -file filename -num run_number -skip 1 (-play | -anim) -auto {limits}  -start time -stop time {slice} -view x|y|z"
   sys.exit()

animate = False
if re.search(r"-anim",commandline) or re.search(r"-play",commandline):
   animate = True

auto = False
if re.search(r"-auto",commandline):
   auto = True

COM = False
if re.search(r"-com",commandline,re.IGNORECASE) or re.search(r"-center",commandline):
   COM = True


whole = False
if re.search(r"-whole",commandline) or re.search(r"-play",commandline):
   whole = True

with open(filename, 'r+') as fi:
     a=json.load(fi)
     print "number of runs = ", len(a) - 1
     p = a[0]
     coords = a[num+1]

print p
num_dim = p["num_dim"]

if num_dim >= 2:
    view = 'y'
    if num_dim == 3:
       matchobj = re.search(r"-vie\w*\s(.+?)(\s|$)",commandline)
       if (matchobj):
          view = matchobj.group(1)

    if view == 'x' or view == 'X':
        dim_x = 1
        dim_y = 2
    elif view == 'y' or view == 'Y':
        dim_x = 0
        dim_y = 1
    elif view == 'z' or view == 'Z':
        dim_x = 0
        dim_y = 2

    num_parts = p["num_parts"]
    x = zeros((num_parts))
    y = zeros((num_parts))

    stop =   p["num_time_slices"]
    num_time_slices =   p["num_time_slices"]
    matchobj = re.search(r"-stop\w*\s(.+?)(\s|$)",commandline)
    if (matchobj):
             stop = int(matchobj.group(1))
    if stop < 0:
       stop =   p["num_time_slices"]


    clist = ['b','g','r','c','m','y','k']
    numc = len(clist)
    if whole:
        mark_spacing = num_time_slices/(numc) + 1
    else:
        mark_spacing = num_time_slices/(2*numc) + 1
    mark_spec = []
    for c in clist:
        mark_spec.append(c+".")

    if COM:
        cntr_start = center_of_mass(coords,num_parts,num_dim,0)
        cntr_end = center_of_mass(coords,num_parts,num_dim,num_time_slices-1)
        for time_slice in range(start,stop,skip):
            for dim in range(num_dim):
                for part_indx in range(num_parts):
                   coords[time_slice][part_indx][dim] -= ((cntr_end[dim]-cntr_start[dim])/(num_time_slices-1))*time_slice



    #for time_slice in range(64,66,skip):
    for time_slice in range(start,stop,skip):
        for part_indx in range(num_parts):
            x[part_indx] = coords[time_slice][part_indx][dim_x]
            y[part_indx] = coords[time_slice][part_indx][dim_y]
        if animate:
           cla()
        xlow = -0.5*p["L"]
        xhigh = 1.5*p["L"]
    #    xlow = 0
    #    xhigh =p["L"]
        if not auto:
            xlim(xlow,xhigh)
            ylim(xlow,xhigh)
        if False:
            plot(x,y,"b^")
        else:
            #plot(x,y,"b.")
            if whole or time_slice < num_time_slices/2:# and time_slice%mark_spacing == 0:
                T = time_slice/mark_spacing
                plot(x,y,mark_spec[T])
            plot(mean(x),mean(y),"y*")
             
        draw()
        if animate:
           pause(0.0005)

elif num_dim == 1:
    if p["num_parts"] != 1:
        print "num_parts must be 1 if num_dim = 1"
        sys.exit()

    n = p["num_time_slices"]
    xarr = zeros(n)

    for i in range(n):
        xarr[i] = coords[i][0][0]

    plot(xarr)

    draw()


show()
