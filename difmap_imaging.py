import subprocess
import time
logfile = open("difmap.log", "w")

inputfile = "./data/a17078a.cal.fits"
difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)

difmap.stdin.write(("obs " + inputfile + "\n").encode())
#print(difmap.stdout.readline())
#difmap.stdin.flush()
#time.sleep(1)
#print(difmap.stdout.readline())
difmap.stdin.write("exit\n\n".encode())
#difmap.wait()
