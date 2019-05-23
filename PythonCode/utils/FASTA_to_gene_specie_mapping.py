
path  = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

COGS_names = ['COG2602']
def main (path,COGS_names,exte):
    for name in COGS_names:
        res = '{ '
        S_text = open(path+name+'/S_'+exte+'.txt', 'r').read()
        input = open(path+name+"/old_new_names.txt", 'r')
        old_new_names = []
        for line in input:
            old_new_names.append(eval(line))
        old_new_names = old_new_names[0]

        with open(path+name+'/FASTA_'+exte+'.txt','r') as fp:
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
        file = open(path + name+ "/sigma.txt", 'w')
        file.write(str(res))
        file.close()
if __name__ == "__main__":
    main(path,COGS_names,'')
