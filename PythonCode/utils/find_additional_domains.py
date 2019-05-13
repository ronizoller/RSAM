path  = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

name = 'COG3093'
number_of_domains = 2

if name == 'COG3093':
    pattern = "((['D'], None, False)_(['HT'], None, False))_Double-Mode"
    list_of_interesting_substrings = ['COG3093','HTH_XRE','antidote_HigA']
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
            found = False
            for inter in list_of_interesting_substrings:
                if line.find(inter) != -1:
                    found = True
            if not found:
                    domain_presense.update({gene:line})

res = {}
for gene,line in domain_presense.items():
    line = line[::-1]
    index = 0
    counter = 0
    while counter != 2:
        index += 1
        if line[index] == '\t':
            counter += 1
    index += 1
    additional_domain = ''
    while line[index] != '\t':
        additional_domain += line[index]
        index += 1
    additional_domain = additional_domain[::-1]
    if additional_domain in res:
        temp = res[additional_domain]
        res.update({additional_domain:temp+1})
    else:
        res.update({additional_domain:1})

for additional_domain,count in res.items():
    print('%s:  %s' % (str(additional_domain),str(count)))
