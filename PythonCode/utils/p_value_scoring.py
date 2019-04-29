path  = '/Users/ronizoller/Google Drive (ronizo@post.bgu.ac.il)/COGS/'

name = 'COG2602'
pattern = "((['HT'], 'red', True)_(['S', 'D', 'HT'], 'black', None))_Double-Mode"
to_check = 'hitdata_all'

domain_presense = {}

with open(path + name + '/saved_data/results/'+pattern+'/'+to_check+'.txt', 'r') as fp:
    for line in fp:
        start = line.find('[') + 1
        end = line.find(']')
        specie = ''
        while start < end:
            specie = specie + line[start]
            start += 1
        if specie != '':
            if specie not in domain_presense:
                domain_presense.update({specie:False})
            if line.find('MecR1') != -1 or line.find('cl28898') != -1 or line.find('BlaR1') != -1 or line.find('Peptidase_M56') != -1 or line.find('cl29776') != -1:
                domain_presense.update({specie:True})
all_specie = 0
all_Trues = 0
for specie,presence in domain_presense.items():
    all_specie += 1
    if presence:
        all_Trues +=1
    else:
        print('%s' % str(specie))

print('\n\nall: %s\nwith BlaR1:%s' % (str(all_specie),str(all_Trues)))