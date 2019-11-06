#red : Human/animal
#black :

## what is AH?

with open('habitat.txt','r') as fp:
    for line in fp:
        if line.find('(') != -1:
            index = 0
            print(";",end='')
            while index < len(line) and line[index+1] != '(' :
                if (line[index] != ' '):
                    print(line[index],end='')
                else:
                    print('_',end='')
                index += 1

            print(";", end='')
            index = line.find("Habitat:")
            if index == -1:
                print('')
            else:
                index = index + 8
                if (line[index:].find("HP;") != -1) or (line[index:].find("AH;") != -1) or (line[index:].find("AH ") != -1) or (line[index:].find("HNP;") != -1) or (line[index:].find("HO;") != -1) or (line[index:].find("C;") != -1) or (line[index:].find("AP;") != -1) or (line[index:].find("ANP;") != -1) or (line[index:].find("HP ") != -1) or (line[index:].find("HNP ") != -1) or (line[index:].find("HO ") != -1) or (line[index:].find("C ") != -1) or (line[index:].find("AP ") != -1) or (line[index:].find("ANP ") != -1):
                    print("human",end='')
                else:
                    print("soil",end='')
                print("\n")

