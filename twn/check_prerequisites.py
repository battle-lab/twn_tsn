try:
    import numpy as np
except:
    print('numpy package has to be installed in python.')
    raise Exception('twn:numpynotexist - numpy not installed in python')

try:
    import pandas as pd
except:
    print('pandas package has to be installed in python.')
    raise Exception('twn:pandasnotexist - pandas not installed in python')
try:
    import argparse
except:
    print('argparse package has to be installed in python.')
    raise Exception('twn:argparsenotexist - argparse not installed in python')

print('python installation is OK.')
