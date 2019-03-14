with open('/Users/ronizoller/PycharmProjects/TreeReconciliation/trees/COG2602/tax_report.txt','r') as fp:
    for line in fp:
        line = line.replace('	 ', "")
        line = ''.join(line.split())

        index = 0
        i = 0
        while i < 1 and index < len(line):
            if (line[index] == '|'):
                i = i+1
            index = index + 1
        start = index
        while i < 2 and index < len(line):
            if (line[index] == '|'):
                i = i+1
            index = index + 1
        mid = index
        while i < 3 and index < len(line):
            if (line[index] == '|'):
                i = i+1
            index = index + 1
        end = index-1
        if (end != mid) and (line[start:mid-1] != line[mid:end]):
            print (line[start:end])
