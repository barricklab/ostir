#!/usr/bin/python

import os.path
import os, subprocess, time
import tempfile

current_dir = os.path.dirname(os.path.abspath(__file__)) + "/tmp"
if not os.path.exists(current_dir): os.mkdir(current_dir)

debug=0

#Class that encapsulates all of the functions from NuPACK 2.0
class ViennaRNA(dict):

    debug_mode = 0
    RT = 0.61597 #gas constant times 310 Kelvin (in units of kcal/mol)
    param_file = "-P rna_turner1999.par "


    def __init__(self,Sequence_List,material = "rna37"):

        self.ran = 0

        import re
        import random
        import string

        exp = re.compile('[ATGCU]',re.IGNORECASE)

        for seq in Sequence_List:
            if exp.match(seq) == None:
                error_string = "Invalid letters found in inputted sequences. Only ATGCU allowed. \n Sequence is \"" + str(seq) + "\"."
                raise ValueError(error_string)

        self["sequences"] = Sequence_List
        self["material"] = material

        random.seed(time.time())
        with tempfile.NamedTemporaryFile(delete=False) as temp_file_maker:
            self.prefix = temp_file_maker.name

    def mfe(self, strands, constraints, Temp , dangles, outputPS = False):

        self["mfe_composition"] = strands
        
        Temp = float(Temp)
        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        seq_string = "&".join(self["sequences"])

        if constraints is None:
            input_string = seq_string + "\n" 
        else: 
            input_string = seq_string + "\n" + constraints + "\n"

        #write sequence file 
        handle = open(self.prefix,"w")
        handle.write(input_string)
        handle.close()

        #Set arguments
        material = self["material"]

        param_file = "-P rna_turner1999.par "

        if dangles is "none":
            dangles = " -d0 "
        elif dangles is "some":
            dangles = " -d1 "
        elif dangles is "all":
            dangles = " -d2 "

        if outputPS:
            outputPS_str = " "
        else:
            outputPS_str = " --noPS "
            
        if constraints is None:
            args = outputPS_str + dangles + param_file + self.prefix

        else:
            args = outputPS_str + dangles + "-C " + param_file + self.prefix


        #Call ViennaRNA C programs
        cmd = "RNAfold"
        #print(cmd + args)

        output = subprocess.Popen(cmd + args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines = True) #newlines argument was added because subprocess.popen was appending b's to the list output (byte stuff? I dont totally understand)
        std_out = output.communicate()[0]
        #output.tochild.write(input_string)

        while output.poll() is None:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        if debug == 1: print(output.stdout.read())

        #Skip the unnecessary output lines
        #line = output.stdout.read()
        #print(std_out)

        line = std_out
        words = line.split()
        #print("These are the mfe value for " + str(words))
        bracket_string = words[1]
        #print("This is the bracket string " + str(bracket_string))
        (strands,bp_x, bp_y) = self.convert_bracket_to_numbered_pairs(bracket_string)
        #print(strands, bp_x, bp_y)

        energy = float(words[len(words)-1].replace(")","").replace("(","").replace("\n",""))

        self["program"] = "mfe"
        self["mfe_basepairing_x"] = [bp_x]
        self["mfe_basepairing_y"] = [bp_y]
        self["mfe_energy"] = [energy]
        self["totalnt"]=strands
        try:
            self._cleanup()
        except FileNotFoundError:
            print("file already removed")

        #print "Minimum free energy secondary structure has been calculated."

    def subopt(self, strands, constraints, energy_gap, Temp = 37.0, dangles = "all", outputPS = False):

        self["subopt_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")
        
        param_file = "-P rna_turner1999.par "

        seq_string = "&".join(self["sequences"])
        
        if constraints is None:
            input_string = seq_string + "\n" 
        else: 
            input_string = seq_string + "\n" + constraints + "&........." "\n"

        handle = open(self.prefix,"w")
        handle.write(input_string)
        handle.close()

        #Set arguments
        material = self["material"]

        if dangles is "none":
            dangles = " -d0 "
        elif dangles is "some":
            dangles = " -d1 "
        elif dangles is "all":
            dangles = " -d2 "

        if outputPS:
            outputPS_str = " "
        else:
            outputPS_str = " -noPS "

        #Call ViennaRNA C programs
        cmd = "RNAsubopt"
        if constraints is None:
            args = " -e " + str(energy_gap) + dangles + param_file + " < " + self.prefix
        else:
            args = " -e " + str(energy_gap) + dangles + "-C " + param_file + " < " + self.prefix

        #print(cmd + args)
        output = subprocess.Popen(cmd + args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines = True) #newlines argument was added because subprocess.popen was appending b's to the list output (byte stuff? I dont totally understand)
        std_out = output.communicate()[0]
        #output.tochild.write(input_string)

        while output.poll() is None:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        if debug == 1: print(output.stdout.read())

        #Skip unnecessary line
        #line = output.std_out.readline()
        line = std_out
        #print(std_out)

        self["subopt_basepairing_x"] = []
        self["subopt_basepairing_y"] = []
        self["subopt_energy"] = []
        self["totalnt"]=[]
        counter=0

        if len(line) > 0:
            words = line.split()
            for i in range(3,len(words)):
                counter+=1
                if i % 2 == 0:
                    energy = float(words[i].replace("\n",""))
                    #print(energy)
                    self["subopt_energy"].append(energy)

                else:
                    bracket_string = words[i]
                    #print(bracket_string)
                    (strands,bp_x, bp_y) = self.convert_bracket_to_numbered_pairs(bracket_string)
                    self["subopt_basepairing_x"].append(bp_x)
                    self["subopt_basepairing_y"].append(bp_y)

        self["subopt_NumStructs"] = counter

        self._cleanup()
        self["program"] = "subopt"

        #print "Minimum free energy and suboptimal secondary structures have been calculated."

    def energy(self, strands, base_pairing_x, base_pairing_y, Temp = 37.0, dangles = "all"):

        self["energy_composition"] = strands

        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")

        seq_string = "&".join(self["sequences"])
        strands = [len(seq) for seq in self["sequences"]]
        bracket_string = self.convert_numbered_pairs_to_bracket(strands,base_pairing_x,base_pairing_y)
        input_string = seq_string + "\n" + bracket_string + "\n"

        #print(seq_string,strands, bracket_string, input_string)

        handle = open(self.prefix,"w")
        handle.write(input_string)
        handle.close()

        #Set arguments
        material = self["material"]
        if dangles is "none":
            dangles = " -d0 "
        elif dangles is "some":
            dangles = " -d1 "
        elif dangles is "all":
            dangles = " -d2 "

        #Call ViennaRNA C programs
        cmd = "RNAeval"
        args = dangles + self.prefix
        #print(cmd + args)
        output = subprocess.Popen(cmd + args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines = True) #newlines argument was added because subprocess.popen was appending b's to the list output (byte stuff? I dont totally understand)

        while output.poll() is None:
            try:
                output.wait()
                time.sleep(0.001)
            except:
                break

        #if debug == 1: print output.fromchild.read()
        std_out = output.communicate()[0]
        #print(std_out)

        self["energy_energy"] = []        

        line = std_out
        words = line.split()
        #print(line, words)
        if words:
            energy = float(words[-1].replace("(","").replace(")","").replace("\n",""))
        else:
            energy = 0  # TODO: Handle this better

        self["program"] = "energy"
        self["energy_energy"].append(energy)
        self["energy_basepairing_x"] = [base_pairing_x]
        self["energy_basepairing_y"] = [base_pairing_y]
        self._cleanup()

        return energy

    def convert_bracket_to_numbered_pairs(self,bracket_string):

        all_seq_len = len(bracket_string)
        #print(all_seq_len)
        bp_x = []
        bp_y = []
        strands = []

        for y in range(bracket_string.count(")")):
            bp_y.append([])
            #print(bp_y)

        last_nt_x_list = []
        counter=0
        num_strands=0
        #print(bracket_string)
        for (pos,letter) in enumerate(bracket_string[:]):
            if letter is ".":
                counter += 1

            elif letter is "(":
                bp_x.append(pos-num_strands)
                last_nt_x_list.append(pos-num_strands)
                counter += 1

            elif letter is ")":
                nt_x = last_nt_x_list.pop() #nt_x is list of "(" except last entry
                #print('this is the last_nt_x_list ' + str(last_nt_x_list.pop()))
                nt_x_pos = bp_x.index(nt_x) 
                bp_y[nt_x_pos] = pos-num_strands
                counter += 1

            elif letter is "&":
                strands.append(counter)
                counter=0
                num_strands+=1

            else:
                print("Error! Invalid character in bracket notation.")
    

        if len(last_nt_x_list) > 0:
            print("Error! Leftover unpaired nucleotides when converting from bracket notation to numbered base pairs.")

        strands.append(counter)
        if len(bp_y) > 1:
            bp_x = [pos+1 for pos in bp_x[:]] #Shift so that 1st position is 1
            bp_y = [pos+1 for pos in bp_y[:]] #Shift so that 1st position is 1
            #print("subopt bp_x " + str(bp_x))
            #print("subopt bp_y " + str(bp_y))
        
        return (strands,bp_x, bp_y)

    def convert_numbered_pairs_to_bracket(self,strands,bp_x,bp_y):

        bp_x = [pos-1 for pos in bp_x[:]] #Shift so that 1st position is 0
        bp_y = [pos-1 for pos in bp_y[:]] #Shift so that 1st position is 0

        bracket_notation = []
        counter=0
        for (strand_number,seq_len) in enumerate(strands):
            if strand_number > 0: bracket_notation.append("&")
            for pos in range(counter,seq_len+counter):
                if pos in bp_x:
                    bracket_notation.append("(")
                elif pos in bp_y:
                    bracket_notation.append(")")
                else:
                    bracket_notation.append(".")
            counter+=seq_len

        return "".join(bracket_notation)

    def _cleanup(self):

        if os.path.exists(self.prefix): 
            os.remove(self.prefix) 
        return

if __name__ == "__main__":

    sequences = ["TCTGGCAGGGACCTGCACACGGATTGTGTGTGTTCCAGAGATGATAAAAAAGGAGTTAGTCTTGGTATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTG"] #,"acctcctta"]
    len(sequences) #,"acctcctta"]
    const_str = None

    test = ViennaRNA(sequences, material = "rna37")
    print(test)
    #dangles = "none"
    #constraints = none

    test.mfe([1], constraints = const_str , dangles = "none", Temp = 37.0)

    bp_x = test["mfe_basepairing_x"][0]
    bp_y = test["mfe_basepairing_y"][0]
    strands = test["totalnt"]
    bracket_string = test.convert_numbered_pairs_to_bracket(strands,bp_x,bp_y)
    print(bracket_string)

    (strands,bp_x, bp_y) = test.convert_bracket_to_numbered_pairs(bracket_string)

    print("Strands = ", strands)
    print("bp_x = ", bp_x)
    print("bp_y = ", bp_y)

    print(test.energy(strands, bp_x, bp_y, dangles = "all"))
    test.subopt(strands,const_str,0.5,dangles = "all")
    print(test)

#    print bracket_string
#    print test.convert_numbered_pairs_to_bracket(strands,bp_x,bp_y)
