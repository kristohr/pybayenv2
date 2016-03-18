#!/local/bin/python

import sys


################## Class Locus ######################
class Locus:
    
    def __init__(self, name, pops):
        
        self.name = name.strip()
        self.pop = [[0 for x in xrange(2)] for x in xrange(pops)]
        self.freqs = [[0 for x in xrange(1)] for x in xrange(pops)]
        self.is_monomorphic = False #The locus is monomorphic for one allele.
        self.is_consensus = False #True if freq differences are pops < cut off.
        self.freq_diff = 0 #Difference in allele freq among pop.

    #Updates the allele frequencies for each population.
    def update_freqs(self, al_type, n):
        if (al_type == HOMOZYGOTE_1):
            self.pop[n][0] += 2
        elif (al_type == HOMOZYGOTE_2):
            self.pop[n][1] += 2
        elif (al_type == MISSING_DATA):
            pass
        else:
            self.pop[n][0] += 1
            self.pop[n][1] += 1
        #print self.pop

    def get_name(self):
        return self.name

    #Prints the name and allele frequency for each locus
    def to_string(self):
        if (self.is_monomorphic_locus()):
            print "Removing monomorphic locus " + self.name + " with the following allele counts:"
            for i in range(0, len(self.pop)):
                curr_pop = self.pop[i]
                print curr_pop

    #Adding the allele freqs to bayenv format (null model)   
    def freqs_to_lines(self):
        line_1 = ""
        line_2 = ""
        if (not self.is_monomorphic_locus()):
            for i in range(0, len(self.pop)):
                line_1 += str(self.pop[i][0]) + "\t"
                line_2 += str(self.pop[i][1]) + "\t"
            line_1 += "\n" + line_2 + "\n"
        return line_1

        #Adding the allele freqs to bayenv format for test purposes
    def freqs_to_lines2(self):
        line_1 = ""
        line_2 = ""
        if (not self.is_monomorphic_locus()):
            line_1 += self.name.replace(":","-") + "\t"
            for i in range(0, len(self.pop)):
                line_1 += str(self.pop[i][0]) + "\t"
                line_2 += str(self.pop[i][1]) + "\t"
            line_1 += "\n" + line_2 + "\n"
        return line_1
    
##      def to_string(self):
##         if (not self.equal_freqs_among_pops()):
##             print self.name
##             for i in range(0, len(self.pop)):
##                 curr_pop = self.pop[i]
##                 print curr_pop


    #Returns true if all pops are zero in one of the alleles.
    def is_monomorphic_locus(self):
        prev_freq = None
        if (self.pop[0][0] == 0): #First allele is zero
            for i in range(0, len(self.pop)):
                curr_freq = self.pop[i][0]
                if (prev_freq is not None and curr_freq != 0):
                    return self.is_monomorphic
                prev_freq = curr_freq
            self.is_monomorphic = True
        elif (self.pop[0][1] == 0): #Second allele is zero
            for i in range(0, len(self.pop)):
                curr_freq = self.pop[i][1]
                if (prev_freq is not None and curr_freq != 0):
                    return self.is_monomorphic
                prev_freq = curr_freq
            self.is_monomorphic = True
        else:
            pass
        return self.is_monomorphic

    #Returns true if the allele frequecy is less than x.
    def equal_freqs_among_pops(self):
        min_diff = 0.01
        prev_freq = None
        is_equal_freq = True
        for i in range(0, len(self.pop)):
            curr_pop = self.pop[i]
            curr_freq = get_allele_freq(curr_pop[0], curr_pop[1])
            if (prev_freq != None):
                diff = abs(curr_freq-prev_freq)
            else:
                diff = 0
            if (prev_freq is not None and diff > min_diff):
                is_equal_freq = False
                break
            prev_freq = curr_freq
        return is_equal_freq
    
    #Calculates the allele freqs for the locus
    def set_freqs(self):
        for i in range(0, len(self.pop)):
            curr_pop = self.pop[i]
            #print self.name
            ####
            self.freqs[i] = get_allele_freq(curr_pop[0], curr_pop[1])
            if self.freqs[i] == 99: #Marking the SNP for deletion
                self.is_consensus = True
                self.is_monomorphic = True
                
            
    #Returns the frequency array for the population
    def get_freqs(self):
        return self.freqs


    #Returns true if the allele frequecy is less than x.
    def is_consensus_locus(self, min_diff):
        is_consensus = True
        for i in range(1, len(self.pop)):
            #print self.freqs
            diff = abs(self.freqs[0]-self.freqs[i])
            if (diff > min_diff):
                is_consensus = False
                break
        return is_consensus

    #Sets the MAFD
    def set_max_freq_diff(self,):
        max_diff = 0
        k = 0
        while (k < len(self.pop)-1):
            for i in range(k+1, len(self.pop)):
                diff = abs(self.freqs[k]-self.freqs[i])
                if (diff > max_diff):
                    max_diff = diff
                    self.freq_diff = diff
            k += 1
        #print "Maximum difference is: " + str(self.freq_diff)

    def get_max_freq_diff(self):
        return self.freq_diff



def convert_fileformat(in_file, testdata):
    #global NULL_FILE
    #global ITERATIONS
    #global NUM_POP
    #global NULL_SIZE
    global TEST_SIZE
    global CUT_OFF
    #global NUM_ENV
    #global NUM_TESTS

    out_file = in_file.split(".")[0] + ".bnv"
    num_pops = str(NUM_POP)
    #Generate loci objects
    locus_list = gen_loci(in_file, num_pops)

            
    freqs_list = []

    TEST_SIZE = len(locus_list)
   
    if (DEBUG):
        debug(locus_list)       
    else:
        for i in range(0, len(locus_list)):
            locus_list[i].set_freqs()
            locus_list[i].set_max_freq_diff()
            #locus_list[i].to_string()
            freqs_list.append(locus_list[i].get_max_freq_diff())
    
        freqs_data = np.array(freqs_list)
            
        count = 0
        lines1 = ""
        for i in range(0, TEST_SIZE):
            lines1 += locus_list[i].freqs_to_lines()
            count += 1
            
        print "\nTotal loci in null model = "+  str(count)

        count = 0
        lines2 = ""
        if (CUT_OFF < EPSILON):   
            for i in range(0, TEST_SIZE):
                lines2 += locus_list[i].freqs_to_lines2()
                count += 1
            print "\nTotal SNPs to test = "+  str(count)
        elif (testdata):  
            print "Computing cutoff data.."
            cut_off = stats.scoreatpercentile(freqs_data, CUT_OFF)
            print "Maximum allele frequency cutoff: " + str(cut_off)
            for i in range(0, TEST_SIZE):
                if (locus_list[i].get_max_freq_diff() > cut_off):
                    lines2 += locus_list[i].freqs_to_lines2()
                    count += 1
            print "Total test snps = "+  str(count)
            print "Total excluded snps = "+  str(TEST_SIZE - count)
        else:  
            cut_off = stats.scoreatpercentile(freqs_data, CUT_OFF)
            for i in range(0, TEST_SIZE):
                if (not locus_list[i].is_consensus_locus(cut_off)):
                    lines2 += locus_list[i].freqs_to_lines2()
                    count += 1
            print "Total test snps = "+  str(count)

        FILE = open(out_file, 'w')
        FILE.write(lines1)
        FILE.close()

        new_out_file = "2"+out_file
        FILE = open(new_out_file, 'w')
        FILE.write(lines2)
        FILE.close()
        TEST_SIZE = count
        return out_file, new_out_file








def main(in_file):


    #convert_fileformat(in_file, True)
    print in_file

    





if __name__ == '__main__':
    
    # Terminate if too few arguments
    if len(sys.argv) < 2:
        print 'usage: %s <ped-file> ' % sys.argv[0]
        sys.exit(-1)
    main(sys.argv[1])
