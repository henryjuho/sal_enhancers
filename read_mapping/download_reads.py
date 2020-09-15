from qsub import q_sub


def main():
    
    counter = 0
    cmds = ['cd /scratch/project_2002047/barson_mapping_v2/reads']

    for line in open('PRJEB10744.txt'):
        
        counter += 1

        if line.startswith('study'):
            continue
        
        reads = line.split('\t')[9].split(';')
        cmds += ['wget -c ftp://' + reads[0], 'wget -c ftp://' + reads[1]]
        
        if counter % 10 == 0:
            q_sub(cmds, out='/scratch/project_2002047/barson_mapping_v2/reads/read_download' + str(counter), t=48, scheduler='SLURM')
            cmds = ['cd /scratch/project_2002047/barson_mapping_v2/reads']
            

    q_sub(cmds, out='/scratch/project_2002047/barson_mapping_v2/reads/read_download' + str(counter), t=48, scheduler='SLURM')

if __name__ == '__main__':
    main()
