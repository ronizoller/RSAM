path  = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

name = 'COG1396'
number_of_domains = 2

if name == 'COG2602':
    list_of_interesting_substrings = ['MecR1','cl28898','BlaR1','Peptidase_M56']
elif name == 'COG3093':
    list_of_interesting_substrings = ['Peptidase_M78','cl01076']
elif name == 'COG1396':
    pattern = "(['HT'], None, False)_Single-Mode"
    list_of_interesting_substrings = ['cupin','Cupin_2']
    to_check = 'hitdata_all'
elif name == 'COG3550':
    pattern = "(['D'], None, False)_Single-Mode"
    list_of_interesting_substrings = []
    to_check = 'hitdata_1th_solution'

domain_presense = {}

with open(path + name + '/saved_data/results/'+pattern+'/'+to_check+'.txt', 'r') as fp:
    for line in fp:
        start = line.find('[') + 1
        end = line.find(']')
        gene_start = line.find('>') + 1
        specie = ''
        index = gene_start
        gene = ''
        while index < start - 2:
            gene = gene + line[index]
            index += 1
        while start < end:
            specie = specie + line[start]
            start += 1
        if specie != '':
            if list_of_interesting_substrings != []:
                if specie not in domain_presense:
                    domain_presense.update({specie:False})
                for inter in list_of_interesting_substrings:
                    if line.find(inter) != -1:
                        domain_presense.update({specie:True})
            else:
                if gene not in domain_presense:
                    domain_presense.update({gene:1})
                else:
                    temp = domain_presense[gene]
                    domain_presense.update({gene: temp + 1})
all = 0
all_Trues = 0
for specie,presence in domain_presense.items():
    all += 1
    if list_of_interesting_substrings != []:
        if presence:
            all_Trues +=1
        else:
            print('%s' % str(specie))
    else:
        if presence < number_of_domains:
            all_Trues +=1
            print('%s has < %s domains' % (str(specie),str(number_of_domains)))


print('\n\nall: %s\nwith Trues:%s' % (str(all),str(all_Trues)))