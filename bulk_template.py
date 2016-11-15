'''
Created on 3 Aug 2016

@author: ajones7
'''
from datetime import datetime
from dateutil.relativedelta import relativedelta
from dateutil.parser import parse


class decipher():

    def __init__(self):
        self.inputfile = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160715 Decipher Bulk Upload\\161115_db_out.txt"
        #self.inputfile = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160715 Decipher Bulk Upload\\prenatal_db_out.txt"

        self.outputfile = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160715 Decipher Bulk Upload\\script_output\\161115_bulk_upload.txt"
        #self.outputfile = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160715 Decipher Bulk Upload\\script_output\\prenatal_cleaned.txt"

        self.mean_ratio_dict = {"x0": "-2", "x0~1": "-2", "x1": "-1",
                                "x1~2": "-1", "x2": "0.58", "x2~3": "0.58", "x3": "0.58", "x4": "1", "x2~4": "1"}
        self.pathogenicity = {"Pathogenic / abnormal result (class 5)":"Definitely pathogenic","Abnormal result (class 5)": "Definitely pathogenic", "Abnormal result (retrospectively assigned)": "Definitely pathogenic", "Not in use: Pathogenic (retrospectively assigned)": "Definitely pathogenic", "Pathogenic, likely (retrospectively assigned)": "Probably pathogenic", "Likely to be pathogenic (class 4)": "Probably pathogenic",
                              "Unknown (retrospectively assigned)": "Uncertain", "Uncertain clinical significance (class 3)": "Uncertain", "Unlikely to be pathogenic (class 2)": "Likely benign", "Benign, likely (retrospectively assigned)": "Likely benign"}

        self.hpo_translation_file = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160715 Decipher Bulk Upload\\hpo translation.txt"
        self.hpo_dict = {}
        self.all_results = {}

    def translate(self):
        # populate hpo_dict
        hpo_translation_file = open(self.hpo_translation_file, 'r')
        for i, line in enumerate(hpo_translation_file):
            if i != 0:
                line = line.rstrip()
                splitline = line.split("\t")
                phenoformID = splitline[0]
                HPO_ID = splitline[1]
                self.hpo_dict[phenoformID] = HPO_ID

        # open input file
        read_input = open(self.inputfile, 'r')
        # loop through line at a time
        for number, line in enumerate(read_input):
            if number != 0:
                #split line on tab
                splitline = line.split("\t")
                PRU = splitline[0]
                chr = splitline[1]  # changed to text in function below
                start = splitline[2]
                end = splitline[3]
                genome_assembley = splitline[4]
                copy_number = splitline[5]  # changed to ratio below
                # not used - calculated using is_male/is_female
                Chromosomal_sex = splitline[6]
                Openaccess_consent = splitline[7]
                DOB = splitline[8]
                Note = splitline[9]
                Mother_is = splitline[10]
                Father_is = splitline[11]
                Inheritance = splitline[12]
                Pathogenicity = splitline[13]  # changed using dict below
                Phenotypes = splitline[14]  # not used - phenoformID used below
                Responsible_contact = splitline[15]
                #print splitline[16]
                Requested_date = parse(
                    splitline[16])
                is_male = splitline[17]
                is_female = splitline[18]
                is_unknown = splitline[19]
                # used to create hpo numbers
                phenoformID = splitline[20].rstrip()
                # print splitline[21], len(splitline[21])
                gestation = None
                if len(splitline[21]) > 3:
                    gestation = splitline[21].rstrip()
                # print gestation
                
                #exclude aberation if has these pathogenicities
                if Pathogenicity in ("Below array resolution", "<1Mb (targeted array,not pathogenic,not reported)"):
                    pass
                else:
                    #ignore if the gender is unknown. Prenatals without a gender are excluded as part of the query
                    if is_female == "0" and is_male == "0" and is_unknown == "0" and chr != "24":
                        pass
                    else:
                        # catch prenatal females 
                        # NB for prenatals is_male == is_female == is_unknown
                        if is_male in ("F F", "F", "f"):
                            sex="46XX"
                        # catch males
                        elif is_male in ("M") or is_male == "-1":
                            sex = "46XY"
                            #catch female post natals
                        elif is_female == "-1" or is_unknown == "-1":
                            sex = "46XX"
                        # if unknown but on Y can say they are male
                        elif chr == "24":
                            sex = "46XY"
                        # else raise error
                        else: 
                            raise ValueError("Unknown gender: is female: "+is_female+" is_male: "+is_male+" is_unknown: "+is_unknown+"\n"+line)
                        #print sex
                        
                        # convert DOB into date
                        if DOB != '':
                            DOB = datetime.strptime(splitline[8], '%d/%m/%Y')
    
                        # get chromosome_text
                        chromosome = decipher().get_chromosome(chr)
    
                        # set inheritance for mosaics
                        if "~" in copy_number:
                            Inheritance = "De novo mosaic"
                        else:
                            pass
    
                        # create a name for imbalance
                        chr_start_end = chromosome + ":" + start + "-" + end
    
                        # pass copy number to get_ratio
                        cleaned_ratio = decipher().get_ratios(
                            copy_number, is_male, is_female, chromosome)
    
                        # calculate the age in years subtracting DOB from date of
                        # test
                        if DOB is not None:
                            rdelta = relativedelta(Requested_date, DOB)
                            Age = rdelta.years
                        else:
                            Age = ""
    
                        # set prenatal age to empty string
                        Prenatal_age = ""
                        # set predefined gestations based on sample type
                        if gestation == "Chorion":
                            Prenatal_age = "12"
                            #overwrite age
                            Age = "Prenatal"
                            # set all other sample types to 20 weeks
                        elif gestation in ("Amnio", "DNA for Cyto", "Fetal Blood", "No Specimen", "Blood", "DNA"):
                            Prenatal_age = ""
                            Age = "Prenatal"
                        # else is postnatal so has no gestation
                        elif gestation is None:
                            Prenatal_age = ""
                        else:
                            raise ValueError("unknown sample type:\t" + gestation)
    
                        # hardcode the responsible contact
                        Responsible_contact = "Viapath Genetics Laboratories"
    
                        # get HPO term from dict
                        if phenoformID in self.hpo_dict:
                            hpo_term = self.hpo_dict[phenoformID]
                        else:
                            hpo_term = ''
    
                        # decipher pathogenicity categories
                        if Pathogenicity in self.pathogenicity:
                            decipher_pathogenicity = self.pathogenicity[
                                Pathogenicity]
                        else:
                            raise ValueError(
                                "unknown pathogenicity:\t" + Pathogenicity)
    
                        to_add_to_dict = [PRU, chromosome, start, end, genome_assembley, copy_number, cleaned_ratio, sex, Openaccess_consent,
                                          Age, Prenatal_age, Note, Mother_is, Father_is, Inheritance, decipher_pathogenicity, [hpo_term], Responsible_contact, Requested_date]
    
                        # check if patient is already in dict
                        if PRU in self.all_results:
                            # create a list to hold all imbalances
                            list_of_imbalances = []
                            # for each imbalance for that patient add id
                            # (chr_start_end)to list
                            for imbalance in self.all_results[PRU]:
                                list_of_imbalances.append(imbalance)
                            for imbalance in list_of_imbalances:
                                # check if this imbalance is the same imbalance as one
                                # already in the dict(does it overlap? and same copy
                                # number)
                                if chromosome == self.all_results[PRU][imbalance][1] and end >= self.all_results[PRU][imbalance][2] and start <= self.all_results[PRU][imbalance][3] and copy_number == self.all_results[PRU][imbalance][5]:
                                    # now we know it's the same imbalance:
    
                                    # 1. the most likely cause is there are multiple
                                    # HPO terms so we can append the HPO term to the
                                    # list
                                    if hpo_term in self.all_results[PRU][imbalance][16]:
                                        pass
                                    else:
                                        self.all_results[PRU][imbalance][
                                            16].append(hpo_term)
    
                                    # 2. we also need to check which breakpoints to use we will use the most recent coordinates
                                    # if requested date is older or same as the
                                    # existing imbalance pass
                                    if Requested_date <= self.all_results[PRU][imbalance][18]:
                                        pass
                                    else:
                                        # if newer use the newer coords
                                        self.all_results[PRU][imbalance][2] = start
                                        self.all_results[PRU][imbalance][3] = end
    
                                # it's a different imbalance in same patient
                                else:
                                    self.all_results[PRU][
                                        chr_start_end] = to_add_to_dict
    
                        else:
                            # add to dictionary a record for this PRU
                            # the value is a dictionary coordinates as a key and the
                            # value a list
                            self.all_results[PRU] = {chr_start_end: to_add_to_dict}

    def get_ratios(self, copy_number, is_male, is_female, chr):
        ''' using dictionary turn the copy number into a predefined log ratio'''
        # if male and on Sex chrom (is_unknown only has 46XX so no need to use
        # this module )
        if chr in ("Y", "X") and is_male == -1:
            # set amended log ratio
            if copy_number in self.mean_ratio_dict:
                ratio_cleaned = self.mean_ratio_dict[copy_number]
                return ratio_cleaned
            elif copy_number == "other":
                pass
            else:
                raise ValueError("unknown copy number:\t" + copy_number)
        else:
            # continue as normal
            if copy_number in self.mean_ratio_dict:
                ratio_cleaned = self.mean_ratio_dict[copy_number]
                return ratio_cleaned
            elif copy_number == "other":
                pass
            else:
                raise ValueError("unknown copy number:\t" + copy_number)

    def get_chromosome(self, chr):
        '''returns the text chromosome name eg X and Y not 23 and 24'''
        if chr == '23':
            return "X"
        elif chr == '24':
            return"Y"
        else:
            return chr

    def write_output(self):
        outputfile = open(self.outputfile, 'w')
        # for each patient
        for i in self.all_results:
            # check there is at least 1 aberation for that patient
            if len(self.all_results[i]) > 0:
                # loop through aberations for that patient
                for j in self.all_results[i]:
                    # capture values
                    PRU = self.all_results[i][j][0]
                    chromosome = self.all_results[i][j][1]
                    start = str(self.all_results[i][j][2])
                    end = str(self.all_results[i][j][3])
                    genome_assembley = self.all_results[i][j][4]
                    copy_number = self.all_results[i][j][5]
                    cleaned_ratio = str(self.all_results[i][j][6])
                    sex = self.all_results[i][j][7]
                    Openaccess_consent = self.all_results[i][j][8]
                    Age = str(self.all_results[i][j][9])
                    Prenatal_age = self.all_results[i][j][10]
                    Note = self.all_results[i][j][11]
                    Mother_is = self.all_results[i][j][12]
                    Father_is = self.all_results[i][j][13]
                    Inheritance = self.all_results[i][j][14]
                    decipher_pathogenicity = self.all_results[i][j][15]
                    hpo_term = self.all_results[i][j][16]
                    Responsible_contact = self.all_results[i][j][17]
                    # tab
                    tab = "\t"
                    # start a string to concatenate HPO terms
                    hpo_string = ''

                    # loop through the list of hpo terms
                    for k in range(len(hpo_term)):
                        if len(hpo_term[k]) > 5:
                            hpo_string = hpo_string + hpo_term[k] + ", "

                    hpo_string = hpo_string.rstrip(", ")

                    outputfile.writelines(PRU + tab + chromosome + tab + start + tab + end + tab + genome_assembley + tab + cleaned_ratio + tab + sex + tab + Openaccess_consent + tab + Age +
                                          tab + Prenatal_age + tab + Note + tab + Mother_is + tab + Father_is + tab + Inheritance + tab + decipher_pathogenicity + tab + hpo_string + tab + Responsible_contact + "\n")

        outputfile.close()

if __name__ == "__main__":
    a = decipher()
    a.translate()
    a.write_output()
    print "done"
