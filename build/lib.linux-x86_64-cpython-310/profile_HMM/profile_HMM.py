import pickle as pkl
from collections import defaultdict
#custom imports
from HMMSTR_utils.HMMSTR_utils import rev_comp
from importlib_resources import files

class ProfileHMM:
    def __init__(self, name,alphabet,prefix,repeat,suffix,background=None,out=".",transitions=None, mismatch_probs=None, repeat_probs =None,flanking_size=30):
        self.out = out
        self.name = name
        self.alphabet = alphabet
        self.flanking_size = flanking_size
        #initialize forward target characteristics
        self.prefix = prefix.upper()[-flanking_size:]
        self.repeat = repeat.upper()
        self.suffix = suffix.upper()[:flanking_size]
        #initialize reverse target characteristics
        self.rev_prefix = rev_comp(prefix.upper())[:flanking_size] #opposite subsetting for revcomp
        self.rev_repeat = rev_comp(repeat.upper())
        self.rev_suffix = rev_comp(suffix.upper())[-flanking_size:]
        self.get_states()
        
        #If no background inputted, default to uniform distribution
        est_072022_transitions_path = files('profile_HMM.data').joinpath('est_072022_transitions.pkl')
        est_072022_emissions_path = files('profile_HMM.data').joinpath('est_072022_emissions.pkl')
        if background is None:
            self.background = {'A':0.25,'T':0.25,'C':0.25,'G':0.25,'':0}
        else:
            self.background = background
        if transitions is None:
            #read default pkl file
            with open(est_072022_transitions_path,"rb") as file: #assumes defaults are in same directory as scripts, may need to change this
                self.transitions = pkl.load(file)
        else:
            self.transitions = transitions
        if mismatch_probs is None:
            with open(est_072022_emissions_path, "rb") as file:
                self.emissions = pkl.load(file)
        else:
            self.emissions = mismatch_probs
        if repeat_probs is None:
            self.repeat_probs = None
        else:
            self.repeat_probs = repeat_probs

    def get_states(self):
        '''
        Function to calculate the number of match states needed
        Returns:
            states (int): number of match states

        '''
        states = 0
        states += len(self.prefix) + len(self.repeat) + len(self.suffix)
        self.states = states

    def build_profile(self, strand):
        '''
        Function to initialize a Profile HMM structure

        Args:
            strand: which strand to build model with reference to
        Returns:
            None, writes HMM file in format needed by umdhmm-v1.02

        '''
        #set prefix, suffix, repeat according to strand
        if strand == "forward":
            prefix = self.prefix
            repeat = self.repeat
            suffix = self.suffix
            if self.repeat_probs is not None:
                repeat_probs = self.repeat_probs
            else:
                repeat_probs = None
        else:
            prefix = self.rev_suffix
            repeat = self.rev_repeat
            suffix = self.rev_prefix
            if self.repeat_probs is not None:
                    conversion = {'A':'T','G':'C','T':'A','C':'G','':''}
                    conversion_states = dict(zip(list(self.repeat_probs.keys()),list(self.repeat_probs.keys())[::-1]))
                    repeat_probs = {}
                    for key in self.repeat_probs.keys():
                        repeat_probs[key] = {}
                        for bp in self.repeat_probs[key]:
                            repeat_probs[key][bp] = self.repeat_probs[conversion_states[key]][conversion[bp]]
            else:
                repeat_probs = None
        #Initialize empty matrices A, E, and I
        A = defaultdict(dict)
        E = defaultdict(dict)
        I = {}
        hidden_states = []

        #initialize 'G' state at start to absorb bases before prefix starts
        hidden_states.append('G1')

        #initialize prefix states
        for i in range(1,len(prefix)+1):
            hidden_states.append('PD'+str(i))
            hidden_states.append('PI'+str(i))
            hidden_states.append('PM'+str(i))

        #initialize repeat states
        for i in range(1,len(repeat)+1):
            hidden_states.append('RD'+str(i))
            hidden_states.append('RI'+str(i))
            hidden_states.append('RM'+str(i))

        #initialize suffix states
        for i in range(1,len(suffix)+1):
            #we don't want an insert in the last position
            hidden_states.append('SD'+str(i))
            if i < len(suffix):
                hidden_states.append('SI'+str(i))
            hidden_states.append('SM'+str(i))
        #initialize 2nd 'G' state
        hidden_states.append('G2')

        hidden_states.append('E') #end node

        # Initialize all transitions to 0
        for item in hidden_states:
            for next_item in hidden_states:
                A[item][next_item] = 0.0

        #initialize start state
        for item in hidden_states:
            I[item] = 0.0

        #need to test probabilities here
        I['G1'] = 1.0

        #initialize emission probabilities
        E = self.get_emission(E, prefix, repeat, suffix,repeat_probs)
        A = self.get_transitions(A, prefix, repeat, suffix)
        if strand == "reverse":
            out = self.out + "_" + self.name + "_revcomp"
            self.rev_hidden_states = hidden_states
            self.rev_A = A
            self.rev_E = E
            self.rev_I = I
        else:
            out = self.out + "_" + self.name
            self.hidden_states = hidden_states
            self.A = A
            self.E = E
            self.I = I
        #call with current A,E,I since then we don't have to change the call based on strand but we have them saved correctly
        self.write_HMM(hidden_states, A, E, I, out)

        return   
    def get_emission(self,E, prefix, repeat, suffix,repeat_probs):
        '''
        Function to initialize emission probability matrix

        Args:
            E (defaultdict): empty default dictionary to fill
            prefix (str): prefix from reference genome to match to
            repeat (str): repeat motif to match
            suffix (str): suffix from reference to match to
            repeat_probs (dictionary): custom repeat probability matrix from user input (if applicable)
        Returns:
            E (defaultdict): filled emissions matrix

        '''
        #initialize 'G' states to have background nucleotide frequencies
        E['G1'] = self.background
        E['G2'] = self.background

        #i don't know if I need this but I'm filling this out as silent just in case
        for bp in self.alphabet:
            E['E'][bp] = 0.0
        #initialize prefix emissions based on prefix sequence
        for i in range(1, len(prefix)+1):
            for bp in self.alphabet:
                E['PD' + str(i)][bp] = 0.0 #deletion states don't emit a base
                E['PI' + str(i)][bp] = self.background[bp]
                if prefix[i-1] == bp:
                    E['PM' + str(i)][bp] = self.emissions[bp][bp]#probability of match
                elif bp == "":
                    E['PM' + str(i)][bp] = 0.0
                else:
                    E['PM' + str(i)][bp] = self.emissions[prefix[i-1]][bp] #probability of mismatch error
            E['PD' + str(i)][""] = 0.25 #try with background freqs

        #initialize repeat emissions based on repeat sequence
        #if custom probabilities given, add dictionary to full emissions dictionary, else continue as usual with emission probabilities
        if repeat_probs is None:
            for i in range(1,len(repeat)+1):
                for j,bp in enumerate(self.alphabet):
                    E['RD' + str(i)][bp] = 0.0 #deletion states don't emit a base
                    E['RI' + str(i)][bp] = self.background[bp] #might need to update this when I update background frequencies

                    if repeat[i-1] == 'N':
                        if bp in ['A','C','T','G']:
                            E['RM' + str(i)][bp] = 0.25 #set ptobaility to 0.25 for all bases since it can be any of them
                        else:
                            E['RM' + str(i)][bp] = 0.0
                    elif repeat[i-1] == bp:
                        E['RM' + str(i)][bp] = self.emissions[bp][bp] #probability of match
                    elif bp == "":
                        E['RM' + str(i)][bp] = 0.0
                    else:
                        E['RM' + str(i)][bp] = self.emissions[repeat[i-1]][bp] #probability of mismatch error
                E['RD' + str(i)][""] = 0.25 #background freqs
        else: #load dictionary and proceed with all other state types as usual
            E.update(repeat_probs) #add ALL match state probs for the repeat region
            for i in range(1,len(repeat)+1):
                for j,bp in enumerate(self.alphabet):
                    E['RD' + str(i)][bp] = 0.0 #deletion states don't emit a base
                    E['RI' + str(i)][bp] = self.background[bp] #might need to update this when I update background frequencies
                E['RD' + str(i)][""] = 0.25 #try with background freqs
        #initialize the suffix emissions based on the suffix sequence
        for i in range(1,len(suffix)+1):
            for j,bp in enumerate(self.alphabet):
                E['SD' + str(i)][bp] = 0.0 #deletion states don't emit a base
                if i < len(suffix): #skip last insertion states since its replaced by 'G2'
                    E['SI' + str(i)][bp] = self.background[bp] #might need to update this when I update background frequencies
                if suffix[i-1] == bp:
                    E['SM' + str(i)][bp] = self.emissions[bp][bp] #probability of match
                elif bp == "":
                    E['SM' + str(i)][bp] = 0.0
                else:
                    E['SM' + str(i)][bp] = self.emissions[suffix[i-1]][bp] #probability of mismatch error
            E['SD' + str(i)][""] = 0.25 #try with background freqs
        
        return E
    def get_transitions(self,A, prefix, repeat, suffix):
        '''
        Function to fill out all non-zero transitions
        Args:
            A (defaultdict): initialized transition matrix
            prefix (str): prefix from reference genome to match to
            repeat (str): repeat motif to match
            suffix (str): suffix from reference to match to


        Returns:
            A (defaultdict): filled out transition matrix

        '''
        #initialize 'G' transitions ***THIS IS A GUESS, MAY NEED TO STILL DECREASE
        A['G1']['G1'] = 0.50 
        A['G1']['PD1'] = 0.15
        A['G1']['PM1'] = 0.35

        #really unsure about these probabilities ***MAY NEED TO STILL CHANGE
        A['G2']['G2'] = 0.50
        A['G2']['E'] = 0.50

        #prefix transitions
        for i in range(1, len(prefix)+1): #all but the last Transitions
            #transitions into insertion states
            A['PI' + str(i)]['PI' + str(i)] = self.transitions["P_ii"]
            A['PD' + str(i)]['PI' + str(i)] = self.transitions["P_di"]
            A['PM' + str(i)]['PI' + str(i)] = self.transitions["P_mi"]
            #transition into deletion states
            if i > 1:
                A['PD'+str(i-1)]['PD'+str(i)] = self.transitions["P_dd"]  # from previous deletion
                A['PM'+str(i-1)]['PD'+str(i)] = self.transitions["P_md"]
                A['PI'+str(i-1)]['PD'+str(i)] = self.transitions["P_id"]
            #transition into match states
            if i > 1:
                A['PD'+str(i-1)]['PM'+str(i)] = self.transitions["P_dm"] # from previous deletion
                A['PM'+str(i-1)]['PM'+str(i)] = self.transitions["P_mm"]
                A['PI'+str(i-1)]['PM'+str(i)] = self.transitions["P_im"]
        #prefix border transitions
        A['PD' + str(len(prefix))]['RM1'] = self.transitions["P_dm"]
        A['PD' + str(len(prefix))]['RD1'] = self.transitions["P_dd"]
        A['PD' + str(len(prefix))]['PI' + str(len(prefix))] = self.transitions["P_di"]

        A['PI' + str(len(prefix))]['PI' + str(len(prefix))] = self.transitions["P_ii"]
        A['PI' + str(len(prefix))]['RM1'] = self.transitions["P_im"]
        A['PI' + str(len(prefix))]['RD1'] = self.transitions["P_id"]

        A['PM' + str(len(prefix))]['RM1'] = self.transitions["P_mm"]
        A['PM' + str(len(prefix))]['RD1'] = self.transitions["P_md"]
        A['PM' + str(len(prefix))]['PI' + str(len(prefix))] = self.transitions["P_mi"]

        #repeat transitions
        for i in range(1, len(repeat)+1): #all but the last Transitions
            #transitions into insertion states
            A['RI' + str(i)]['RI' + str(i)] = self.transitions["P_ii"]
            A['RD' + str(i)]['RI' + str(i)] = self.transitions["P_di"]
            A['RM' + str(i)]['RI' + str(i)] = self.transitions["P_mi"]
            #transition into deletion states
            if i > 1:
                A['RD'+str(i-1)]['RD'+str(i)] = self.transitions["P_dd"]  # from previous deletion
                A['RM'+str(i-1)]['RD'+str(i)] = self.transitions["P_md"]
                A['RI'+str(i-1)]['RD'+str(i)] = self.transitions["P_id"]
            #transition into match states
            if i > 1:
                A['RD'+str(i-1)]['RM'+str(i)] = self.transitions["P_dm"] # from previous deletion
                A['RM'+str(i-1)]['RM'+str(i)] = self.transitions["P_mm"]
                A['RI'+str(i-1)]['RM'+str(i)] = self.transitions["P_im"]
        #border and looping transitions (I haven't looked at this yet and it probably follows some distribution but thats unclear rn)
        A['RM' + str(len(repeat))]['SM1'] = 0.42
        A['RM' + str(len(repeat))]['SD1'] = 0.01
        A['RM' + str(len(repeat))]['RI' + str(len(repeat))] = 0.015
        #loop back from end of repeat
        A['RM' + str(len(repeat))]['RM1'] = 0.54
        A['RM' + str(len(repeat))]['RD1'] = 0.015
        #from deletion state
        A['RD' + str(len(repeat))]['SD1'] = 0.002
        A['RD' + str(len(repeat))]['SM1'] = 0.43
        A['RD' + str(len(repeat))]['RI' + str(len(repeat))] = 0.004
        #loop back from deletion
        A['RD' + str(len(repeat))]['RD1'] = 0.004
        A['RD' + str(len(repeat))]['RM1'] = 0.56
        #from insertion
        A['RI' + str(len(repeat))]['RI' + str(len(repeat))] = 0.004
        A['RI' + str(len(repeat))]['SM1'] = 0.43
        A['RI' + str(len(repeat))]['SD1'] = 0.002
        #loop back from insertion
        A['RI' + str(len(repeat))]['RM1'] = 0.56
        A['RI' + str(len(repeat))]['RD1'] = 0.004

        #suffix transitions
        for i in range(1, len(suffix)+1): #all but the last Transitions
            #transitions into insertion states
            if i != len(suffix): # there is no last insertion state
                A['SI' + str(i)]['SI' + str(i)] = self.transitions["P_ii"]
                A['SD' + str(i)]['SI' + str(i)] = self.transitions["P_di"]
                A['SM' + str(i)]['SI' + str(i)] = self.transitions["P_mi"]
            #transition into deletion states
            if i > 1:
                A['SD'+str(i-1)]['SD'+str(i)] = self.transitions["P_dd"]# from previous deletion
                A['SM'+str(i-1)]['SD'+str(i)] = self.transitions["P_md"]
                A['SI'+str(i-1)]['SD'+str(i)] = self.transitions["P_id"]
            #transition into match states
            if i > 1:
                A['SD'+str(i-1)]['SM'+str(i)] = self.transitions["P_dm"] # from previous deletion
                A['SM'+str(i-1)]['SM'+str(i)] = self.transitions["P_mm"]
                A['SI'+str(i-1)]['SM'+str(i)] = self.transitions["P_im"]
        #border suffix transitions
        A['SM' + str(len(suffix))]['G2'] = 1.0
        A['SD' + str(len(suffix))]['G2'] = 1.0
        return A   
    def write_HMM(self, hidden_states, A, B, pi, out):
        '''
        Function to write out Profile HMM in format needed by umdhmm
        Args:
            alphabet (list): alphabet characters for the model
            A (defaultdict(defaultdict)): transition matrix of HMM
            B (defaultdict(defaultdict)): emission probability matrix of HMM
            pi (defaultdict): initial probability matrix of HMM
            out (str): output file prefix to save to

        Returns:
            None, saves to out.hmm and out.hidden_states.txt
        '''
        M = len(self.alphabet)#length of alphabet
        N = self.states*3 + 2 #number of states including 'G', 'E', insertions and deletions
        with open(out + ".hidden_states.txt",'w') as outfile:
            outfile.write(".".join(hidden_states))
        with open(out + ".hmm", 'w') as outfile:
            outfile.write("M= " + str(M) + "\n")
            outfile.write("N= " + str(N) + "\n")
            outfile.write("A:\n")
            #loop across A to get transition matrix
            #use hidden states order since states were added in order of graph
            for curr_key in hidden_states:
                for curr_item in hidden_states: #loop across all keys in the inner dictionary
                    outfile.write(str(A[curr_key][curr_item])+ " ")
                outfile.write("\n")
            #write out the emission probability matrix
            outfile.write("B:\n")
            for curr_key in hidden_states:
                for bp in self.alphabet:
                    outfile.write(str(B[curr_key][bp])+ " ")
                outfile.write("\n")
            #write out pi
            outfile.write("pi:\n")
            for curr_key in hidden_states:
                outfile.write(str(pi[curr_key]) + " ")
            outfile.write("\n")
        return
    
        