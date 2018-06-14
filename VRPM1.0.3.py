#!/usr/bin/env python3.6
from collections import Counter
import os

# Program Parsing Through Ensemble VEP Results (txt output), performs some statistical
# analysis and merges key results back with the original VCF file submitted

# Was written in various stages and much of the code is brute force coding style
# Program became a kind of frankstein program that probably can be simplified, 
# but it works and over run time is very fast. 

# Parts of Program
# 1. Analyzes the number of genes that are found in a VEP file by gene symbol,
#       ensemble gene ID and protein ID
# 2. Analyzes the type of consequences found in a VEP file (quantifies counts) and
#       determines how far downstream variants are from a gene (verify extended 3'UTR)
# 4. Merge VEP results with Original VCF file analyzed

# Note: Program is kind of a Frankenstein/Brute force program that where different parts were
# developed at different times and then all functionalities were merged into one program to simplify
# a graduate students life. Program runs pretty quickly, so further cleaning up and development
# was postponed to another day.

def parsing_input_parameter_file():
    """ 
    Parses through the paramter file to return a dictionary
    of the input parameters
    :param none: automatically opens file
    :return dictionary: dictionary of file names for analysis
    """

    input_file = open('VRPM_Parameter_File.txt', 'r')

    print ("Parsed Lines from User")
    
    for line in input_file:
        if line.startswith("--"):
            line=line.rstrip('\n')
            # just printing user inpt for validation
            print (line)
            parsed_parameters=line.split("--")

            for x in range(1, len(parsed_parameters)):
                inputs = parsed_parameters[x].split(" ")
                
                if inputs[0] == "VCF_Input_File":
                    vcf_file_name = inputs[1] 

                elif inputs[0] == "VEP_Input_File":
                    vep_input_file = inputs[1]

    # Skip anything else
    else:
        pass

    input_file.close

    #Printing 
    print ("Name of VCF input File is: ", vcf_file_name)
    print ("Name of VEP Input File is: ", vep_input_file)
    
    return{'vcf_file_name':vcf_file_name, 'vep_input_file':vep_input_file}


def average(list_values):
    """
    Finds the average of a list (no tricks)
    :param list_values: list of values to average
    :return average_numb: average number
    """
    list_sum = 0
    for value in list_values:
        list_sum = list_sum + int(value)
    average_numb = list_sum/len(list_values)
    return average_numb


def max_value(list_values):
    """
    Finds the max value from a list of values
    :param list_values: list of values
    :return max_numb: returns the max number of a list
    """
    # need to initialize starting value
    max_numb = int(list_values[0])

    for value in list_values:

        if max_numb > int(value):
            pass
        else:
            max_numb = int(value)

    return max_numb


def min_value(list_values):
    """
    Finds the min value from a list of values
    :param list_values: list of values
    :return min_numb: returns the min number of a list
    """
    min_numb = int(list_values[0])

    for value in list_values:

        if min_numb < int(value):
            pass
        else:
            min_numb = int(value)

    return min_numb


def parsing_vep_consequences(filename):
    """
    Parsing VEP consequence output and prints a summary Report
    of the consquences. Similiar to overview diagrams of consquences
    found on their website
    :param filename: name of file being analyzed
    :return dictionary: returns a dictionary of upstream and downstream variant results

    """

    # opening the file
    vep_file = open(filename, 'r')

    # creating new output files of upstream/downstream results
    upstream_file = open('upstream_variants_LOG.txt', 'w')
    downstream_file = open('downstream_variants_LOG.txt', 'w')

    # Creating Summary Report of Data
    summary_report = open('summary_report.txt', 'w')

    # Counters for Types of Consequences (all)
    all_variants_counter = 0
    synonymous_variant_counter = 0
    downstream_gene_variant_counter = 0
    three_prime_utr_variant_counter = 0
    intron_variant_counter = 0
    missense_variant_counter = 0
    upstream_gene_variant_counter = 0
    intergenic_variant_counter = 0
    five_prime_utr_variant_counter = 0
    splice_region_variant_counter = 0
    others_counter = 0
    non_coding_transcript_exon_variant = 0


    # Creating lists to store values of distances
    upstream_bin = []
    downstream_bin = []

    for line in vep_file:
        if line.startswith('#'):
            # getting the header of the file
            header_vep = line

            # print header to new files
            upstream_file.write(header_vep)
            downstream_file.write(header_vep)

        else:
            # Counting number of variants in file
            all_variants_counter += 1

            # Parsing actual data line
            parsed_line = line.rstrip().split('\t')

            # Setting up various filters/counters of the data 
            # getting variant count
            if parsed_line[2] == 'synonymous_variant':
                synonymous_variant_counter += 1
                continue

            # Getting variant counts and writing to file
            elif parsed_line[2] == 'downstream_gene_variant':
                downstream_file.write(line)

                # getting downstream distance values
                downstream_bin.append(parsed_line[19])

                # adding to counter
                downstream_gene_variant_counter += 1

            elif parsed_line[2] == '3_prime_UTR_variant':
                three_prime_utr_variant_counter += 1
                continue

            elif parsed_line[2] == 'intron_variant':
                intron_variant_counter += 1
                continue

            elif parsed_line[2] == 'missense_variant':
                missense_variant_counter += 1
                continue

            elif parsed_line[2] == 'upstream_gene_variant':
                upstream_file.write(line)

                # getting upstream distance values
                upstream_bin.append(parsed_line[19])

                # adding to counter
                upstream_gene_variant_counter += 1

            elif parsed_line[2] == 'intergenic_variant':
                intergenic_variant_counter += 1
                continue

            elif parsed_line[2] == '5_prime_UTR_variant':
                five_prime_utr_variant_counter += 1
                continue

            elif parsed_line[2].startswith('splice_region_variant'):
                splice_region_variant_counter += 1
                continue

            elif parsed_line[2].startswith('non_coding_transcript_exon_variant'):
                non_coding_transcript_exon_variant += 1
                continue

            else:
                others_counter += 1
                continue

    # Summary Report section of code
    summary_report.write("Summary Report of VEP Data\n")
    summary_report.write("\n")
    summary_report.write("Consequences\tTotal\tPercentage\n")
    summary_report.write("Synonymous variants\t" + str(synonymous_variant_counter)
                         + "\t" + str(round(synonymous_variant_counter / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("Downstream gene variants\t" + str(downstream_gene_variant_counter)
                         + "\t" + str(round(downstream_gene_variant_counter / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("3' prime UTR variants\t" + str(three_prime_utr_variant_counter)
                         + "\t" + str(round(three_prime_utr_variant_counter / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("Intron variants\t" + str(intron_variant_counter)
                         + "\t" + str(round(intron_variant_counter / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("Missense variants\t" + str(missense_variant_counter)
                         + "\t" + str(round(missense_variant_counter / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("Upstream gene variants\t" + str(upstream_gene_variant_counter)
                         + "\t" + str(round(upstream_gene_variant_counter / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("Intergenic variants\t" + str(intergenic_variant_counter)
                         + "\t" + str(round(intergenic_variant_counter / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("5' prime UTR variants\t" + str(five_prime_utr_variant_counter)
                         + "\t" + str(round(five_prime_utr_variant_counter / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("Splice region variants\t" + str(splice_region_variant_counter)
                         + "\t" + str(round(splice_region_variant_counter / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("Non-coding transcript exon variants\t" + str(non_coding_transcript_exon_variant)
                         + "\t" + str(round(non_coding_transcript_exon_variant / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("'Other' variants\t" + str(others_counter)
                         + "\t" + str(round(others_counter / all_variants_counter * 100, 2))
                         + "%" + '\n')
    summary_report.write("Total\t" + str(all_variants_counter) + '\n')

    summary_report.write("\n")
    summary_report.write("\n")

    vep_file.close()
    upstream_file.close()
    downstream_file.close()
    summary_report.close()

    return{'upstream_bin': upstream_bin, 'downstream_bin': downstream_bin}


def binning_data(name_of_binning, data_values):
    """
    Function bins the upstream and downstream variant results
    Bins are based on upto 5,000 bps becauses VEPs limit for considering
    a neighboring gene is this far
    :param name_of_binning: upstream or downstream variants
    :param data_values: values being analyzed
    :return NONE 
    """
    bin_one = 0
    bin_two = 0
    bin_three = 0
    bin_four = 0
    bin_five = 0
    for value in data_values:
        if 0 <= int(value) < 1000:
            bin_one += 1
        elif 1000 <= int(value) < 2000:
            bin_two += 1
        elif 2000 <= int(value) < 3000:
            bin_three += 1
        elif 3000 <= int(value) < 4000:
            bin_four += 1
        elif 4000 <= int(value) <= 5000:
            bin_five += 1
        else:
            print("Input error, double check list")

    # printing results
    summary_report = open('summary_report.txt', 'a')
    summary_report.write("Binning "+name_of_binning+" Variants\n")
    summary_report.write("Bins\tCounts\n")
    summary_report.write("(0 - 1000 nt)\t" + str(bin_one) + "\n")
    summary_report.write("(1000 - 2000 nt)\t" + str(bin_two) + "\n")
    summary_report.write("(2000 - 3000 nt)\t" + str(bin_three) + "\n")
    summary_report.write("(3000 - 4000 nt)\t" + str(bin_four) + "\n")
    summary_report.write("(4000 - 5000 nt)\t" + str(bin_five) + "\n")

    average_distance = round(average(data_values), 2)
    minimum_distance = min_value(data_values)
    maximum_distance = max_value(data_values)

    summary_report.write("\n")
    summary_report.write("Average distance from genes (nt): "+str(average_distance)+"\n")
    summary_report.write("Distance Range from genes (nt): "+str(minimum_distance)+"-"+str(maximum_distance))

    summary_report.write("\n")
    summary_report.write("\n")

    summary_report.close()


def print_data_to_file(splitting_data, file_name):
    """
    Function analyzes the counter results and prints the results
    to the file of interest
    :param splitting_data: counter data will all junk removed
    :param file_name: name of file being printed to:
    :return NONE
    """
    for x in range(len(splitting_data)):
        all_data_values = splitting_data[x]
        split_data = all_data_values.split(":")
        gene_name = split_data[0]
        gene_count = split_data[1]
        returned_data = gene_name + "\t" + gene_count + "\n"
        file_name.write(returned_data)


# Counter module adds in tons of junk
def replacing_counter_junk(data):
    """
    Count module adds in tons of junk into 
    the analysis, function removes ALL it (tons of steps)
    :param data: data to get rid of junk
    :return replacement_eight: cleaned up data
    """

    replacement_one = data.replace('{', '')
    replacement_two = replacement_one.replace('}', '')
    replacement_three = replacement_two.replace('{', '')
    replacement_four = replacement_three.replace("'", '')
    replacement_five = replacement_four.replace("Counter", '')
    replacement_six = replacement_five.replace("(", '')
    replacement_seven = replacement_six.replace(")", '')
    replacement_eight = replacement_seven.replace(" ", '')
    return (replacement_eight)


def parsing_vep_results(file_name, y):
    """
    Function parses through a specific column of VEP Results
    to create a list of results which is then tallied with the 
    Counter module (note: tons of junk gets added to tally 
        with this module)
    :param file_name: name of VEP file being analyzed
    :param y: specific column of VEP txt file being parsed_line
    :return counts: counts from that column of interest  

    """

    with open(file_name) as data_file:
        # creating an empty list for storing data
        values = []
        # loop through the file
        for line in data_file:
            # code removes the header files from vcf file
            if line.startswith('#'):
                continue
            else:
                parse_line = line.split("\t")
                data = (parse_line[y])
                if data == '-':
                    values.append("no_ID")
                else:
                    values.append(data)
    counts = str(Counter(values))
    return counts


def merge_final_results(vcf_input_file, vep_input_file):
    """
    Functions opens the vep file (again...), retrieves important variant effect information
    from the file, this information then gets merged back with the original vcf file
    :param vcf_input_file: original vcf file submitted to VEP
    :param vep_input_file: vep output file
    return NONE
    """

    #Creating an empty dictionary to store information
    vep_results_dict ={}

    #Use brute force and get results
    # opening the file
    vep_file = open(vep_input_file, 'r')

    # Parse through VEP results and get all the important data
    for line in vep_file:
        if line.startswith('#'):
            continue
        else:
            # Parsing actual data line
            parsed_line = line.rstrip().split('\t')
            variant_name = parsed_line[1]
            consequence = parsed_line[2]
            impact = parsed_line[3]
            symbol = parsed_line[4]
            gene = parsed_line[5]

            vep_results_dict.update({variant_name: {'variant_name': variant_name,
                                                    'consequence': consequence, 
                                                    'impact': impact,
                                                    'symbol': symbol,
                                                    'gene': gene}})
    #Close the file
    vep_file.close()


    #Now print results to a new file
    annota_vcf_file_name = open("annotated_"+vcf_input_file, 'w')
    
    vcf_file = open (vcf_input_file, 'r')

    # Parse through VEP results and get all the important data
    for line in vcf_file:
        if line.startswith('#'):
            parsed_line = line.rstrip().split('\t')
            header_gene_info = parsed_line[:5]
            header_rest = parsed_line[5:]

            #Convert to header to a string (tab seperated)
            header_gene_info = ('\t'.join(map(str,header_gene_info)))
            
            header_rest = ('\t'.join(map(str,header_rest)))

            annota_vcf_file_name.write(header_gene_info
                                       + "\tConsequence\tImpact\tSymbol\tGene\t"
                                       + header_rest + "\n")
            
        else:
            #Now matching up the variant from the vcf file and merging annotation documentation
            
            parsed_line = line.rstrip().split('\t')

            vcf_variant_name = parsed_line[0] + ":" + parsed_line[1] + "-" + parsed_line[1]
            
            variant_gene_info = parsed_line[:5]
            variant_rest = parsed_line[5:]

            variant_gene_info = ('\t'.join(map(str,variant_gene_info)))
            variant_rest = ('\t'.join(map(str,variant_rest)))

            try:
                variant_consequence = str(vep_results_dict[vcf_variant_name]['consequence'])
                variant_impact = str(vep_results_dict[vcf_variant_name]['impact'])
                variant_symbol = str(vep_results_dict[vcf_variant_name]['symbol'])
                variant_gene = str(vep_results_dict[vcf_variant_name]['gene'])

                annota_vcf_file_name.write(variant_gene_info
                                        + "\t" + variant_consequence + "\t"
                                        + variant_impact + "\t"
                                        + variant_symbol + "\t"
                                        + variant_gene + "\t"
                                        + variant_rest + "\n")
            except KeyError:
                variant_consequence = 'No Entry'
                variant_impact = 'No Entry'
                variant_symbol = 'No Entry'
                variant_gene = 'No Entry'

                annota_vcf_file_name.write(variant_gene_info
                                        + "\t" + variant_consequence + "\t"
                                        + variant_impact + "\t"
                                        + variant_symbol + "\t"
                                        + variant_gene + "\t"
                                        + variant_rest + "\n")


    # Closing the files
    annota_vcf_file_name.close()
    vcf_file.close()


def main():

    #########################################################
    # Opens the parameter file to get file name
    print("Opening the parameter file")
    parameter_stuff = parsing_input_parameter_file()

    # Retrieving the parameters from the parameter file parsing output
    vcf_input_file = parameter_stuff['vcf_file_name']
    vep_input_file = parameter_stuff['vep_input_file']

    
    # Creating the output files for data analysis (three different txt files)
    ensemble_file = open('ensemble_gene_IDs.txt', "w")
    gene_file = open('genes.txt', "w")
    protein_file = open('protein_IDs.txt', "w")

    #########################################################
    # Ensemble_Gene_IDS
    # Creating Header for the output file
    file_header = "Ensemble Gene ID\tCounts\n"
    ensemble_file.write(file_header)

    # Parse the vcf files and counts the  Ensemble_Gene_IDs
    counted_data = parsing_vep_results(vep_input_file, 5)

    # Calls the function to convert counter module output
    # into a more appropriate parsing form (get rid of dictionary junk
    data_no_junk = replacing_counter_junk(counted_data)

    # Splitting the final data for printing results
    splitting_data = data_no_junk.split(",")

    # Call the function to parse through the results and print in a file
    print_data_to_file(splitting_data, ensemble_file)

    # Closing the gene file
    ensemble_file.close()

    #########################################################
    # Gene Symbol
    # Creating Header for the output file
    file_header = "Gene ID\tCounts\n"
    gene_file.write(file_header)

    # Parse the vcf files and counts the Gene Symbols
    counted_data = parsing_vep_results(vep_input_file, 4)

    # Calls the function to convert counter module output
    # into a more appropriate parsing form (get rid of dictionary junk
    data_no_junk = replacing_counter_junk(counted_data)

    # Splitting the final data for printing results
    splitting_data = data_no_junk.split(",")

    # Call the function to parse through the results and print in a file
    print_data_to_file(splitting_data, gene_file)

    # Closing the gene file
    gene_file.close()

    #########################################################
    # Protein_Ensemble_IDs
    # Creating Header for the output file
    file_header = "Protein Ensemble ID\tCounts\n"
    protein_file.write(file_header)

    # Parse the vcf files and counts the protein ensemble IDs
    counted_data = parsing_vep_results(vep_input_file, 7)

    # Calls the function to convert counter module output
    # into a more appropriate parsing form (get rid of dictionary junk
    data_no_junk = replacing_counter_junk(counted_data)

    # Splitting the final data for printing results
    splitting_data = data_no_junk.split(",")

    # Call the function to parse through the results and print in a file
    print_data_to_file(splitting_data, protein_file)

    # Closing the gene file
    protein_file.close()

    #########################################################
    # Analyzing Distance of Variants form Genes
    print("Running Analysis of Variant Distances from Genes")
    data_results = parsing_vep_consequences(vep_input_file)

    upstream_bin = data_results['upstream_bin']
    downstream_bin = data_results['downstream_bin']

    # Dealing with upstream and downstream data
    binning_data('Upstream_Analysis', upstream_bin)
    binning_data('Downstream_Analysis', downstream_bin)
    print("Finished Variant Analysis of Distances")

    #########################################################
    # Merging VEP Results with VCF File
    print ("Merging VEP Results with VCF Files")
    merge_final_results(vcf_input_file, vep_input_file)

    print ("Program Done Running")
    

main()
