import sys
import subprocess

if len(sys.argv) != 0:
    print("usage: ....")
    sys.exit(1)

filename=sys.argv[1]


#consider adding psf, weight_image, testfile outputs... 
subprocess.run(["sex", filename, 
    "-c", "default.sex", 
    "CATALOG_NAME", cataloguename,
])

