import sys

skip_list = ['NC_001960.1', 'KT223110.1_SDY']
for line in sys.stdin:

    if line.startswith('##contig'):
        if 'NW_' in line or 'KT223110.1_SDY' in line or 'NC_001960.1' in line:
            continue

    if line.startswith('NW_') or line.startswith('KT223110.1_SDY') or line.startswith('NC_001960.1'):
        continue

    print(line.rstrip())

