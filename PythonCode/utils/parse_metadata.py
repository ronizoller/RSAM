class genome:
    def __init__(self,vars):
        self.genome_id = vars[0]
        self.genome_name = vars[1].replace('"""','')
        self.organism_name = vars[2]
        self.taxon_id = vars[3]
        self.genome_status = vars[4]
        self.strain = vars[5]
        self.serovar = vars[6]
        self.biovar = vars[7]
        self.pathovar = vars[8]
        self.mlst = vars[9]
        self.other_typing = vars[10]
        self.culture_collection = vars[11]
        self.type_strain = vars[12]
        self.completion_date = vars[13]
        self.publication = vars[14]
        self.bioproject_accession = vars[15]
        self.biosample_accession = vars[16]
        self.assembly_accession = vars[17]
        self.genbank_accessions = vars[18]
        self.refseq_accessions = vars[19]
        self.sequencing_centers = vars[20]
        self.sequencing_status = vars[21]
        self.sequencing_platform = vars[22]
        self.sequencing_depth = vars[23]
        self.assembly_method = vars[24]
        self.chromosomes = vars[25]
        self.plasmids = vars[26]
        self.contigs = vars[27]
        self.sequences = vars[28]
        self.genome_length = vars[29]
        self.gc_content = vars[30]
        self.patric_cds = vars[31]
        self.brc1_cds = vars[32]
        self.refseq_cds = vars[33]
        self.isolation_site = vars[34]
        self.isolation_source = vars[35]
        self.isolation_comments = vars[36]
        self.collection_date = vars[37]
        self.isolation_country = vars[38]
        self.geographic_location = vars[39]
        self.latitude = vars[40]
        self.longitude = vars[41]
        self.altitude = vars[42]
        self.depth = vars[43]
        self.other_environmental = vars[44]
        self.host_name = vars[45]
        self.host_gender = vars[46]
        self.host_age = vars[47]
        self.host_health = vars[48]
        self.body_sample_site = vars[49]
        self.body_sample_subsite = vars[50]
        self.other_clinical = vars[51]
        self.antimicrobial_resistance = vars[52]
        self.antimicrobial_resistance_evidence = vars[53]
        self.gram_stain = vars[54]
        self.cell_shape = vars[55]
        self.motility = vars[56]
        self.sporulation = vars[57]
        self.temperature_range = vars[58]
        self.optimal_temperature = vars[59]
        self.salinity = vars[60]
        self.oxygen_requirement = vars[61]
        self.habitat = vars[62]
        self.disease = vars[63]
        self.comments = vars[64]
        self.additional_metadata = vars[65]

    def find_name(self,name):
        return self.genome_name.replace(' ','').replace('_','').find(name.replace(' ','').replace('_','')) != -1 or\
               name.replace(' ','').replace('_','').find(self.genome_name.replace(' ','').replace('_','')) != -1

path = '/Users/ronizoller/PycharmProjects/RSAM/PythonCode/data/COG2602(class_D,single,all)/'

input1 = open(path + '/taxa_names_array.txt', 'r')
to_tag = []
for line in input1:
    to_tag.append(eval(line))
to_tag = to_tag[0]

genomes = []
counter = 0
with open('genome_metadata', 'r') as fp:
    res = ''
    for line in fp:
        vars = line.split('\t')
        if len(vars) == 66:
            genomes.append(genome(line.split('\t')))

names = [genome.genome_name for genome in genomes]

filtered_genomes = []
for genome in genomes:
    flag = False
    index = 0
    while not flag and index < len(to_tag):
        name = to_tag[index]
        if genome.find_name(name):
            filtered_genomes.append(genome)
            flag = True
        index += 1

tags_from_patric = ''
for genome in filtered_genomes:
    tags_from_patric += ';'+genome.genome_name.replace(' ','_')+';'+genome.gram_stain+' '+genome.isolation_site+' '+genome.host_name+' '+genome.habitat+' '+genome.other_environmental+'\n'

print('was found in patric: '+str(len(filtered_genomes))+'/'+str(len(to_tag)))
file = open(path + "/tags_from_patric.txt", 'w')
file.write(str(tags_from_patric))
file.close()