import sys
import numpy as np

num_channel = 4
CIF = 'CIF'
version = 1
precision = 2
# default of to_bytes is 'false' i.e. big endian
unsigned_short_one = 1

# Initialize parameters
if len(sys.argv) != 4:
    print("Usage: python3 cif_and_pos_writer.py intensity_file lane tile")
    sys.exit()
intensity_file = sys.argv[1]
lane = sys.argv[2]
tile = sys.argv[3]
num_cyc = 0
num_clusters = 0
with open(intensity_file, 'r') as f:
    data = f.readlines()
    for line_num, line in enumerate(data):
        if line == '\n':
            num_cyc = line_num
            print('number of cycles:', num_cyc)
            break
    num_clusters = int(len(data) / (num_cyc + 1))
    print('number of clusters', num_clusters)

# the document from Illumina said cif file contains intensity of signed 2 byte
data = np.zeros([num_clusters, num_cyc, num_channel], dtype=np.int16)
# store intensity data and write to position file
with open(intensity_file, 'r') as f, open('s_{}_{}_pos.txt'.format(lane, tile), 'w') as pos:
    read_data = f.readlines()
    cyc = 0
    cluster = 0
    for line in read_data:
        line = line.replace('\n', '')
        intensities = line.split(' ')
        if len(intensities) < 4:
            continue
        if line[0:3] == '1 1':
            cyc = 1
            cluster += 1
            x, y = intensities[2], intensities[3]
            # write to position file
            pos.write('{} {}\n'.format(x, y))
            intensities = intensities[4:]
        for index, i in enumerate(intensities):
            data[cluster - 1][cyc - 1][index] = round(float(i))
        cyc += 1

print('Reading Complete')

unsigned_short_one = num_cyc

# write to cif and position file
with open('s_{}_{}.cif'.format(lane, tile), 'wb') as cif:
    # write the starting 13 bytes information
    cif.write(bytes([ord(i) for i in CIF]))
    cif.write(version.to_bytes(1, 'little'))
    cif.write(precision.to_bytes(1, 'little'))
    cif.write(cyc.to_bytes((2), 'little'))
    cif.write(unsigned_short_one.to_bytes(2, 'little'))
    cif.write(num_clusters.to_bytes((4), 'little'))
    for cyc in range(1, num_cyc + 1):
        # write intensities
        for i in range(num_channel):
            intensities = data[:, cyc - 1, i]
            for x in intensities:
                cif.write(int(x).to_bytes(2, 'little', signed=True))
