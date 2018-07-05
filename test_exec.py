#!/usr/bin/env python

print("hey it loaded!")
import sys
print(sys.version)
for arg in sys.argv:
    print arg
try:
    import yt
    import trident
    import quasar_scan
    print("import statements worked, including quasar_scan")
except: 
    print("import statements didn't work :(")

