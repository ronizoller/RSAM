def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start


def main(list_of_ids_species,COG_to_find):
    list_of_ids = list(list_of_ids_species.keys())
    with open('species.mappings.v11.0.txt', 'r') as fp:
        specie_cogs = {}
        for line in fp:
            id = line[0:find_nth(line,'\t',0)]
            COG = line[find_nth(line,'\t',0)+1:find_nth(line,'\t',2)]
            if id in specie_cogs and id in list_of_ids:
                specie_cogs[id].append(COG)
            elif id in list_of_ids:
                specie_cogs.update({id:[COG]})
    colors = {}
    for id,COGS in specie_cogs.items():
        if COG_to_find in COGS:
            colors.update({list_of_ids_species[id]:'red'})
        else:
            colors.update({list_of_ids_species[id]:'black'})
    return colors

if __name__ == '__main__':
    main()