import numpy as np 
_value = int.from_bytes(b'\x00\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb\xcc\xdd\xee\xff',byteorder='big',signed=True)
charge = np.empty(39, dtype='uint16')
for i in range(36):
    print(i)