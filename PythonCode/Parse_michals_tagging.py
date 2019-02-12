#red : Human/animal
#black :

## what is AH?

print('{',end=' ')
with open('habitat.txt','r') as fp:
    for line in fp:
        if line.find('(') != -1:
            index = 0
            print("'",end='')
            while index < len(line) and line[index] != '(' :
                if (line[index] != ' '):
                    print(line[index],end='')
                index = index + 1
            print("' :", end=' ')
            index = line.find("Habitat:")
            if index == -1:
                print('** No color **')
            else:
                index = index + 8
                if (line[index:].find("HP;") != -1) or (line[index:].find("AH;") != -1) or (line[index:].find("AH ") != -1) or (line[index:].find("HNP;") != -1) or (line[index:].find("HO;") != -1) or (line[index:].find("C;") != -1) or (line[index:].find("AP;") != -1) or (line[index:].find("ANP;") != -1) or (line[index:].find("HP ") != -1) or (line[index:].find("HNP ") != -1) or (line[index:].find("HO ") != -1) or (line[index:].find("C ") != -1) or (line[index:].find("AP ") != -1) or (line[index:].find("ANP ") != -1):
                    print("'red'",end='')
                else:
                    print("'black'",end='')
                print(",\n")
print('}',end=' ')

