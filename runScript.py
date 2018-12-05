import os
import commands

pwd = commands.getoutput("pwd")
script ="analyzeHexaNoise.py"
path="/home/hep/giulia_plots/noise_91/soldered/"
fl=commands.getoutput("ls "+path+"| grep txt").split('\n')

arg=[os.path.join(path, f) for f in fl]
print arg
for a in arg:
    os.system("python "+script+" "+a)

