# Script for batch sumbission of Asymptotic Limit Calculator
#!/usr/bin/env python
import sys,os,commands,glob,numpy
from math import log,exp
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-w","--workspace",dest="workspace")
parser.add_option("-m","--mass",dest="mass",default=120.,type='float')
parser.add_option("-n","--npoints",dest="npoints",default=20,type='int')
parser.add_option("-r","--range",dest="range",nargs=2,type='float',default=(0.,2.))
parser.add_option("-l","--log",dest="log",action="store_true",default=False)
parser.add_option("","--runLimit",dest="runLimit",action="store_true",default=False)
parser.add_option("","--nCPU",dest="nCPU",type='int',default=1)
parser.add_option("-o","",dest="outputname",default="")
parser.add_option("-d","--directory",dest="directory",default="",type='str')
parser.add_option("-O","",dest="options",default="")

(options,args)=parser.parse_args()

ws=options.workspace
mass=options.mass
outputname=options.outputname
print options.range
step=(options.range[1]-options.range[0])/options.npoints
points = []
if options.log :
	max = options.range[1]
	min = options.range[0]
    	dx = log(max/min)/(options.npoints-1) if options.npoints > 1 else 0
    	points = [ min * exp(dx*i) for i in range(options.npoints) ]
else: 
	points= numpy.arange(options.range[0],options.range[1]+step,step)if options.npoints>1 else [float(options.range[0])]

workingDir=os.getcwd()
wslocation=workingDir+"/"+ws
outputdir =options.directory

if options.directory: jobname =  "%s/limitgrid_%.1f%s"%(outputdir,mass,outputname)
else :
	jobname =  "limitgrid_%.1f%s"%(mass,outputname)
	outputdir = workingDir
f = open("%s.sh"%(jobname),"w")

# Create job script
f.write("#!/bin/bash\n")
f.write("set -x\n")
f.write("cd %s\n"%workingDir)
f.write("eval `scramv1 runtime -sh`\n")
		
f.write("cd -\n")
f.write("cp -p $CMSSW_BASE/lib/slc5_amd64_gcc472/libHiggsAnalysisCombinedLimit.so . \n" )
f.write("mkdir scratch\n")
f.write("cd scratch\n")
f.write("cp -p %s . \n" % (wslocation) )

# Loop over number of points, and add wait after npoints/nCPU jobs
# if running with nCPU > 1, submit with bsub -n N -R "span[hosts=1]" 
ext = " & " if options.nCPU>1 else " " 
if "/" in ws: wsfile = ws.split("/")[-1]
else: wsfile = ws
for i,po in enumerate(points):
  f.write("combine -M Asymptotic %s -m %.1f --singlePoint %f -n %g %s %s \n" % (wsfile,mass,po,po,options.options,ext) )
  if (i>0 and i%options.nCPU==0): f.write("wait\n")
#if options.nCPU>i: f.write("wait\n")
if options.nCPU>0: f.write("wait\n")
f.write("hadd -f grid_%.1f%s.root higgsCombine* \n"%(mass,outputname))
if options.runLimit:
  f.write("rm higgsCombine* \n")
  f.write("combine %s -M Asymptotic -m %.1f --getLimitFromGrid  grid_%.1f%s.root  \n"%(wsfile,mass,mass,outputname))
  f.write("hadd -f res_%.1f%s.root higgsCombine* \n"%(mass,outputname))
  f.write("cp -f res_%.1f%s.root %s \n"%(mass,outputname,outputdir))

f.write("cp -f grid_%.1f%s.root %s \n"%(mass,outputname,outputdir))
f.write("echo 'DONE' \n")
if options.nCPU>1: print "Running with multiple cores, use bsub -n %d -R 'span[hosts=1]'"%options.nCPU
# make it executable
os.system("chmod 755 %s.sh"%(jobname))

