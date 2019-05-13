path  = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

name = 'COG3093'
number_of_domains = 2

if name == 'COG2602':
    list_of_interesting_substrings = ['MecR1','cl28898','BlaR1','Peptidase_M56','BlaR']
    pattern = "((['HT'], 'red', True)_(['S', 'D', 'HT'], 'black', False))_Double-Mode"
    to_check = 'hitdata_all'
elif name == 'COG3093':
    pattern = "((['D'], None, False)_(['HT'], None, False))_Double-Mode"
    list_of_interesting_substrings = ['COG3093','HTH_XRE','antidote_HigA']
    to_check = 'hitdata_all'
elif name == 'COG1396':
    pattern = "(['HT'], None, True)_Single-Mode"
    #list_of_interesting_substrings = ['cupin','Cupin_2']
    list_of_interesting_substrings = ['Peptidase_M78','cl01076','ImmA']
    to_check = 'hitdata_all'
elif name == 'COG3550':
    pattern = "(['D'], None, False)_Single-Mode"
    list_of_interesting_substrings = []
    to_check = 'hitdata_all'

domain_presense = {}

def check_counting (counter):
    all = 0
    all_Trues = 0
    for gene,count in counter.items():
        all += 1
        if count < 2:
            print('%s has < %s domains' % (str(gene), str(number_of_domains)))
            all_Trues += 1
    return all,all_Trues

def extract_from_to (line,end):
    while not line[end].isdigit():
        end += 1
    while not line[end] == '	':
        end += 1
    fro = ''
    end += 1
    while line[end].isdigit():
        fro += line[end]
        end += 1
    end += 1
    to = ''
    while line[end].isdigit():
        to += line[end]
        end += 1
    return fro,to

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
        if gene != '':
            if list_of_interesting_substrings != []:
                found = False
                for inter in list_of_interesting_substrings:
                    if line.find(inter) != -1:
                        found = True
                if not found:
                        domain_presense.update({gene:line})
            else:
                fro,to = extract_from_to(line,end)
                if gene not in domain_presense:
                    temp = []
                    to_append = [int(fro), int(to)]
                    domain_presense.update({gene: [(to_append)]})
                else:
                    temp = list(domain_presense[gene])
                    to_add = temp.append([int(fro), int(to)])
                    domain_presense.update({gene:temp})

for gene,line in domain_presense.items():
    print(line)
