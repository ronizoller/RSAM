with open('/Users/ronizoller/Dropbox/ComparativeOperonAnalysis/Proposal/Test/BetaLactamaseCOGs/COG4219(315 genes in 173 species)/tax_report.txt','r') as fp:
    for line in fp:
        line = line.replace('	 ', "")
        line = ''.join(line.split())

        index = 0
        i = 0
        while i < 3 and index < len(line):
            if (line[index] == '|'):
                i = i+1
            index = index + 1
        print (line[index:])
