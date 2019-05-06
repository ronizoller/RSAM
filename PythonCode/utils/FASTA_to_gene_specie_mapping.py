
path  = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

COGS_names = ['COG1396']
for name in COGS_names:
    res = '{ '
    S_text = open(path+name+'/phyliptree(binary,deltaepsilon).phy', 'r').read()
    input = open(path+name+"/old_new_names.txt", 'r')
    old_new_names = []
    for line in input:
        old_new_names.append(eval(line))
    old_new_names = old_new_names[0]

    with open(path+name+'/FASTA_bacteria.txt','r') as fp:
        for line in fp:
            if line[0] == '>':
                start = 1
                end = line.find(' ')
                gene = ''
                while start < end :
                    gene = gene + line[start]
                    start += 1
                gene = gene.replace("_"," ")
                start  = line.find('[')+1
                end  = line.find(']')
                specie = ''
                while start < end :
                    specie = specie + line[start]
                    start += 1
                if (S_text.find(old_new_names[specie].replace(' ','')) != -1):
                    res += "'"+gene+"' : '"+old_new_names[specie].replace(' ','')+"',"

    res += ' }'
    file = open(path + name+ "/0/sigma0.0.txt", 'w')
    file.write(str(res))
    file.close()
