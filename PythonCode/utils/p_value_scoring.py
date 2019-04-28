path  = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

COGS_names = ['COG3550']

domain_presense = {}

with open(path + name + '/FASTA.txt', 'r') as fp:
    for line in fp:
        start = line.find('[') + 1
        end = line.find(']')
        specie = ''
        while start < end:
            specie = specie + line[start]
            start += 1
        if specie not in domain_presense:
            domain_presense.update({specie:False})
        if line.find('MecR1') != -1 or line.find('cl21491') != -1 or line.find('BlaR1') != -1:
            domain_presense.update({specie:True})

all_specie = 0
all_Trues = 0
for specie,presence in domain_presense.items():
    all_specie += 1
    if presence:
        all_Trues +=1
    else:
        print('%s has no BlaR1' % str(specie))