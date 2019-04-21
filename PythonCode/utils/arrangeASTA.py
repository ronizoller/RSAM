path  = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

COGS_names = ['COG2856']
for name in COGS_names:
    FASTA_a = open(path+name+'/FASTA_alignedroot.txt', 'r').read()
    res = ''
    with open(path+name+'/FASTA_a.txt','r') as fp:
        for line in fp:
            if line[0] == '>':
                i = 0
                while line[i] != ' ':
                    temp = line[i]
                    temp = temp.replace('Z', '')
                    temp = temp.replace('z', '')
                    temp = temp.replace('X', '')
                    temp = temp.replace('x', '')
                    res += temp
                    i += 1
            else:
                res += line
    file = open(path + name+ "/fixed_FASTA_a.fa", 'w')
    file.write(str(res))
    file.close()
