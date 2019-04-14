with open('/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COG2856(the_awsome)/FASTA.txt','r') as fp:
    for line in fp:
        if line[0] == '>':
            start = line.find('[')+1
            end  = line.find(']')
            while start < end :
                if line[start] == ' ':
                    print(' ',end='')
                else:
                    print(line[start],end='')
                start += 1
            print('\n')
