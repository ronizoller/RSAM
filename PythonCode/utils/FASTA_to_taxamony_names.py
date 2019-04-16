cogs_names = ['COG4679','COG3549','COG4737','COG2856','COG3550']
path = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'
for name in cogs_names:
    with open('/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'+name+'/FASTA.txt','r') as fp:
        res = ''
        for line in fp:
            if line[0] == '>':
                start = line.find('[')+1
                end = line.find(']')
                while start < end:
                    if line[start] == ' ':
                        res += ' '
                    else:
                        res += line[start]
                    start += 1
                res += '\n'
    file = open(path + name + "/taxa_names.txt", 'w')
    file.write(str(res))
    file.close()
