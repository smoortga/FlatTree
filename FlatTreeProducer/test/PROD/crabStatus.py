import os

indir = os.getcwd()
dirs = [i for i in os.listdir(indir) if "crab_" in i and os.path.isdir(indir+"/"+i)]

for d in dirs:
    os.system("crab status -d %s"%d)