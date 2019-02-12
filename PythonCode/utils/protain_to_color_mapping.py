print('{',end=' ')
with open('protain-annotaions.txt','r') as fp:
    for line in fp:
        line = line.replace('\t','&')
        index = 0
        i = 0
        while i < 2 and index < len(line):
            if (line[index] == '&'):
                i = i+1
            index = index + 1
        print("'",end='')
        while index < len(line)-1:
            if (line[index] != ' '):
                print(line[index],end='')
            index = index + 1
        print("' :", end=' ')
        if (not line.find("Beta") == -1) or (not line.find("beta") == -1):
            print("'red'",end='')
        else:
            print("'black'",end='')
        print(",\n")
print('}',end=' ')

