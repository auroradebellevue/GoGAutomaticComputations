#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:31:45 2020

@author: aurora
"""
import subprocess
import Multiplier as mult
import FSA_manipulations as FSA


class graph_of_groups(object):
    """
    This is a graph of groups.
    """

    def __init__(self, v, e, a):
        self.vertices = v  # v should be a list of vertex objects
        self.edges = e  # e should be a list of edge objects
        self.adjacency = a  # directed adjacency matrix of underlying graph
        self.spanning_tree = []
        self.inflated_wa = ""
        self.inflated_wa_file = "InflatedHigginsNF.wa"
        self.wa = ""
        self.wa_file = "Pi1GraphofGroupsWordAcceptor.wa"
        self.gm = ""
        self.gm_file = "somenamethatmakessense.gm" #TODO
        self.mult_list=[]

    def cascade(self, w, file_dir, kbmag_ftn_dir):
        """
        Given a word w, perform the cascade algorithm as outlined in
        many sources :)
        This algorithm must be used in combination with the free reduction 
        function for gog.wa's finite state automaton.

        Parameters
        ----------
        w : string

        Returns
        -------
        u : a string which represents a group element in the fundamental group
        of the graph of groups.

        """
        ####################################################################
        # WARNING!!!!!!!!!!!!!!!!!!
        # THIS CODE ONLY WORKS FOR HNN EXTENSIONS
        ####################################################################
        
        #check for "IdWord"
        if w =="IdWord" or w == str():
            return str()
        
        # determine the number of stable letters
        n = 0  # n is number of stable letters; initialize to 0
        stable_indices = []  # indices of the letter t or T
        for i in range(len(w)):
            if w[i] == "t" or w[i] == "T":
                n += 1
                stable_indices.append(i)
            else:
                n += 0
        #print("input word", w)
        #print("number of stable letters", n)
        #print("list of stable indices", stable_indices)

        # split the word w into a list of words
        u = []
        nums = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
        # if there are stable letters in w
        if len(stable_indices) > 0:
            u.append(w[:stable_indices[0]])
        # else there are no stable letters in w
        else:
            u.append(w)
        for i in range(len(stable_indices)):
            j = stable_indices[i]

            if i == len(stable_indices)-1:
                temp = w[j:]
                stable_subword = str()
                vx_subword = str()
                for i in range(len(temp)):
                    if temp[i] == "t" or temp[i] == "T":
                        stable_subword = stable_subword + temp[i]
                    elif temp[i] in nums:
                        stable_subword = stable_subword + temp[i]
                    else:
                        vx_subword = vx_subword + temp[i]

                u.append(stable_subword)
                u.append(vx_subword)
            else:
                k = stable_indices[i+1]
                temp = w[j:k]
                stable_subword = str()
                vx_subword = str()
                for i in range(len(temp)):
                    if temp[i] == "t" or temp[i] == "T":
                        stable_subword = stable_subword + temp[i]
                    elif temp[i] in nums:
                        stable_subword = stable_subword + temp[i]
                    else:
                        vx_subword = vx_subword + temp[i]

                u.append(stable_subword)
                u.append(vx_subword)
            #print(j, " j ", k, " k ", u)
        #print("list of subwords of w, split at stable letters", u)

        # working right to left, simplify w by using u
        original_word = u.copy()
        word_diff = str()
        rev_index = list(range(len(u)))
        rev_index.reverse()
        for i in rev_index:
            #print("\n", "Index of cascade work", str(i), "\n")
            #print("Current cascaded word", u)
            if "t" not in u[i] and "T" not in u[i]:
                # determine algorithm mechanics based on the edge and stable letter info
                word_diff = u[i]
                #print("Current word in cascade step", word_diff)
                if i != 0:
                    stab = u[i-1]
                    if "t" in u[i-1]:
                        e_num = int(stab.replace("t", ""))
                        # set vx as the terminal vertex for the edge
                        vx = self.edges[e_num].tv
                        stab = "t"
                        gp_file = self.edges[e_num].fg.supergp_file
                    else:
                        e_num = int(stab.replace("T", ""))
                        # set vx as the initial vertex for the edge
                        vx = self.edges[e_num].iv
                        stab = "T"
                        gp_file = self.edges[e_num].rg.supergp_file
                else:
                    vx = 0  # the starting vertex is always 0
                    stab = "IdWord"
                    e_num = None
                    gp_file = self.vertices[vx].group.file['vx']

                # use the assigned gp_file to reduce the word
                new_diff = word_reduce(
                    word_diff, gp_file, file_dir, kbmag_ftn_dir)
                #print("new word difference", new_diff)
        #######################################################################
        # End of the cascading function
        #######################################################################
                # We've completed every cascading reduction if i=0
                if i == 0:

                    # Update u[0] and check for "IdWord" in new_diff
                    if new_diff != "IdWord":
                        u[0] = new_diff
                    # new_diff == "IdWord" and len(u) == 0: #new_diff = IdWord and u has other entries
                    else:
                        u[0] = "IdWord"
                    # convert u (a list of strings) to a single string
                    string_word = list_to_string(u)
                    
                    return string_word

                # Else: we need to find the prefix which is the subgroup word

                # determine which subgroup to be contained in
                subgp_word = str()
                supgp_word = str()
                if stab == "t":
                    sbgp_alpha = self.edges[e_num].fg.gens + \
                        self.edges[e_num].fg.invgens

                else:  # stab == "T"
                    sbgp_alpha = self.edges[e_num].rg.gens + \
                        self.edges[e_num].rg.invgens

                # find prefix of of new word diff that is in the subgroup
                if new_diff != "IdWord":
                    found_end = False
                    index = int(0)

                    while found_end != True:
                        # print("Error line","new_diff", new_diff, "new_diff[index]",
                        #new_diff[index],"index", index,"type of index",
                        #type(index),"subgroup alphabet",
                        # sbgp_alpha,"FoundEnd", found_end)
                        if new_diff[index] in sbgp_alpha:
                            if index != len(new_diff)-1:
                                # split new diff into subgp_word * supgp_word by adding
                                # a letter each loop to subgp and deleting one from supgp
                                subgp_word = subgp_word + new_diff[index]
                                supgp_word = new_diff[index+1:]
                               #print("subgroup word ", subgp_word,"supergroup word ", supgp_word)
                                index += 1
                            else:
                                # all of new_diff is a subgp_word
                                subgp_word = subgp_word + new_diff[index]
                                supgp_word = str()
                                found_end = True
                        else:
                            found_end = True
                            # whole word is in supergroup
                            supgp_word = new_diff[index:]

                else:
                    subgp_word = str()
                #print("Subgroup word", subgp_word, "Supergroup word", supgp_word)

                # push subgroup word into the next piece of u for the next
                # cascade step and leave the remainder in u[i]
                u[i-1] = u[i-1]+subgp_word
                u[i] = supgp_word
            # Case where t or T is in u[i]
            else:
                #print("use the stable function to flip subgroup words")
                # print(u[i])
                # print(self.edges[0].isom)

                temp = u[i]
                # delete the number - it's always 0 for now
                # TODO: HNN only implementation
                for inum in range(len(nums)):
                    temp = temp.replace(nums[inum], str())
                if "t" in u[i]:
                    # use the forward isomorphism to flip the suffix of
                    # u[i] across the stable letter and into u[i-1]
                    temp = temp.replace("t", str())
                    #print("Subgroup word part 2 \nSubword without stable letter", temp)
                    new_temp = temp.translate(self.edges[0].isom)
                    #print("Subword after isomorphism", new_temp)
                else:  # T in u[i]
                    # use the backwards isomorphism to flip the suffix of
                    # u[i] across the stable letter and into u[i-1]
                    temp = temp.replace("T", str())
                    #print("Subgroup word part 2 \nSubword without stable letter", temp)
                    new_temp = temp.translate(self.edges[0].isom)
                    #print("Subword after isomorphism", new_temp)
                # update u[i] and u[i-1]
                u[i-1] = u[i-1]+new_temp
                u[i] = original_word[i]
            # end of else
        # end of for loop - return statement is at line 135


class vertex(object):
    """
    This is a vertex of a graph of groups.
    """

    def __init__(self, label, group):
        self.label = label
        self.group = group
        self.subgp = "e0gp-ind" #This subgroup is only defined for the edge group going into v0


class edge(object):
    """
    This is a edge and the reverse edge of a graph of groups.
    """

    def __init__(self, label, iv, tv, fG, rG):
        self.label = label
        self.iv = iv  # initial vertex
        self.tv = tv  # terminal vertex
        self.fg = fG  # forward oriented subgroup
        self.rg = rG  # reverse oriented subgroup
        self.stab_ltr = ["t" + label, "T" + label]
        self.coset_block_files = []
        self.lL_subword_file = []
        self.isom = {}  # isomorphism
        # self.inv_isom=self.inv_iso() #inverse isomorphism

    """
    The inverse isomorphism can be stored in the same translation dictionary
    in Python. 
    
    def inv_iso(self):
        print("todo inverse isoms")
        #TODO figure this out on Linux
        """


class group(object):
    """
    This is a group.
    """

    def __init__(self, gens, invgens, eqns, file):
        self.gens = gens
        self.invgens = invgens
        self.eqns = eqns
        self.ordering = "shortlex"
        self.file = {"vx": file}
        self.file_lines = {"vx": str()}
        self.wa = ""
        self.wa_file = file+".wa"
        self.gm = ""
        self.gm_file = file+".gm"
        self.mult_list = []

    def get_from_file(self, file, file_dir):
        f = open(file_dir+file, "r")
        f_lines = f.readlines()
        self.file_lines = f_lines
        # find important lines in group kbmag file
        gen_line = [i for i in f_lines if 'generatorOrder' in i]
        inv_line = [i for i in f_lines if 'inverses' in i]
        eqn_lines = [i for i in f_lines if '[' in i and ']' in i and '*' in i]

        self.gens = process_line(gen_line, "letters")
        self.invgens = process_line(inv_line, "letters")
        for line in eqn_lines:
            self.eqns.append(process_line(line, "equations"))

        f.close()

    def get_label(self):
        file_name = self.file["vx"]
        index = file_name.find("gp")
        return file_name[1:index]  # typical file name: v0gp


def word_reduce(w, file, file_dir, kbmag_ftn_dir):
    """
    Calls kbmag's word reduce function on the word w using the 
    RWS in file+file_dir

    Parameters
    ----------
    w : string

    Returns
    -------
    u : string

    """
    fmtd_word = str()
    # if l!="IdWord":
    #    new_word=new_word+l.swapcase()
    if w != "IdWord":
        for i in range(len(w)):
            if fmtd_word != str():
                fmtd_word = fmtd_word+"*"+w[i]
            else:
                fmtd_word = w[i]
    else:
        fmtd_word = "IdWord"
    # if r!="IdWord":
    #    if new_word!=str():
    #        new_word=new_word+"*"+r
    #    else:
    #        new_word=r
    # if new_word==str():
    #    new_word="IdWord"

    fmtd_word = fmtd_word + ';'
    #print("formatted word", fmtd_word, '\n')
    #print(file_dir+self.file["vx"], '\n')
    # print(kbmag_ftn_dir+"wordreduce")
    output = subprocess.run([kbmag_ftn_dir+"wordreduce",
                             file_dir + file],
                            input=fmtd_word.encode(),
                            capture_output=True)

    txt_out = output.stdout.decode("utf-8")
    # print(txt_out)
    if "reduces to:" in txt_out:
        index1 = txt_out.index("reduces to:")
        suffix = txt_out[index1:]
        suffix = suffix.replace("\n", "")
        suffix = suffix.replace(" ", "")
        suffix = suffix.replace("reducesto:", "")
        # print("Reduced word:", suffix, ". Variable Type", type(suffix),
        #      ". Is a space in the saved reduced word?", bool(" " in suffix))
        new_word = str()
        if suffix == "IdWord":
            new_word = "IdWord"
            suffix = str()
        while suffix != str():
            if "*" in suffix:
                index3 = suffix.index("*")
                temp = suffix[:index3]
                #print(suffix,"   ", temp)
                if "^" in temp:
                    index4 = temp.index("^")
                    letter = temp[:index4]
                    # print(letter)
                    number = int(temp[index4+1:])
                    for _ in range(number):
                        new_word = new_word+letter
                else:
                    new_word = new_word+temp
                suffix = suffix[index3+1:]
                #print("new suffix", suffix)
            else:
                temp = suffix
                if "^" in suffix:
                    index4 = temp.index("^")
                    letter = temp[:index4]
                    # print(letter)
                    number = int(temp[index4+1:])
                    for _ in range(number):
                        new_word = new_word+letter
                    suffix = str()
                else:
                    new_word = new_word+suffix
                    suffix = str()
    else:
        #print("kbmag output:", txt_out, '\n')
        new_word = 'ERROR: word reduce function failed'
    #print("new word", new_word)
    return(new_word)


def list_to_string(list_wd):
    str_wd = str()
    for i in range(len(list_wd)):
        str_wd = str_wd+list_wd[i]
    str_wd = str_wd.replace("IdWord", str())
    if len(str_wd) == 0:
        str_wd == str()
    return str_wd


class subgroup(object):
    """
    This is a subgroup of a group.
    """

    def __init__(self, gens, file, super_file):
        self.gens = gens
        self.invgens = []
        self.file = file
        self.supergp_file = super_file
        self.coswa = ""
        self.coswa_file = file.replace(".sub", "")+".cos.wa"
        self.cosgm = ""
        self.cosgm_file = file.replace(".sub", "")+".cos.gm"
        self.mult_list = []

    def get_from_file(self, file, file_dir):
        f = open(file_dir+file, "r")
        file_lines = f.readlines()
        # find important lines in group kbmag file
        gen_line = [i for i in file_lines if 'subGenerators' in i]
        self.gens = process_line(gen_line, "letters")
        self.invgens = self.get_inv()
        f.close()

    def get_inv(self):
        inv_list = []
        for i in self.gens:
            inv_list.append(i.swapcase())
        return inv_list

    def get_label(self):
        file_name = self.file
        index1 = file_name.find("gp")
        return file_name[1:index1]  # typical file name: e0gp-f.sub

    # Broken as written
    # def get_isom_sgp(self, isom):
    #    old_gens=self.gens
    #    new_gens=old_gens
    #    for i in range(len(old_gens)):
    #        new_gens[i]=old_gens[i].translate(isom)
    #    return subgroup(new_gens)


def make_gog(vfiles, efiles, init_term_vs, isoms, adj, file_dir):
    """
    Full GOG ready
    Returns
    -------
    Graph of Groups object

    """
    vertices = []
    edges = []
    # populate vertex objects

    for vfile in vfiles:
        gv = group([], [], [], vfile)
        gv.get_from_file(vfile, file_dir)
        label = gv.get_label()
        v = vertex(label, gv)
        vertices.append(v)

    # populate edge objects
    for efile in efiles:
        # efile is a list
        # [forward edge super group file, forward edge file,
        # reverse edge super group file, reverse edge file]
        ge = subgroup([], efile[1], efile[0])
        ge_rev = subgroup([], efile[3], efile[2])
        ge.get_from_file(efile[1], file_dir)
        ge_rev.get_from_file(efile[3], file_dir)
        label = ge.get_label()
        e = edge(label, init_term_vs[efiles.index(efile)][0],
                 init_term_vs[efiles.index(efile)][1], ge, ge_rev)
        # get isomorphisms
        index0 = efiles.index(efile)
        e.isom = isoms[index0]
        edges.append(e)

    # create Graph of Groups Object
    gog = graph_of_groups(vertices, edges, adj)
    return gog


def other_orderings(gog):
    # TODO: this function is not needed anymore since the user should create the
    # additional files needed to describe the different orderings
    for e in gog.edges:
        # find terminal and initial vertices
        for v in gog.vertices:
            if int(v.label) == e.tv:
                term_vx = v
            if int(v.label) == e.iv:
                init_vx = v
        # create new file for terminal vertex with an ordering that has the
        # e letters before the other letters
        f_sub_all = e.fg.gens+e.fg.invgens
        gp_all = term_vx.group.gens.copy()
        sub_gens_list = []
        index_list = []
        for i in gp_all:
            if i in f_sub_all:
                sub_gens_list.append(i)
                index_list.append(gp_all.index(i))
        print("forward subgroup gens:", f_sub_all, "group gens:", gp_all)
        index_list.reverse()
        # pop in the opposite order
        for i in index_list:
            gp_all.pop(i)
        fsub_gen_first = f_sub_all+gp_all
        print(fsub_gen_first)
        # create new file for intial vertex with an ordering that has the
        # e reverse letters before the other letters
        r_sub_all = e.rg.gens+e.rg.invgens
        gp_all = init_vx.group.gens.copy()
        sub_gens_list = []
        index_list = []
        for i in gp_all:
            if i in r_sub_all:
                sub_gens_list.append(i)
                index_list.append(gp_all.index(i))
        print("reverse subgroup gens:", r_sub_all, "group gens:", gp_all)
        index_list.reverse()
        # pop in the opposite order
        for i in index_list:
            gp_all.pop(i)
        rsub_gen_first = r_sub_all+gp_all
        print(rsub_gen_first)
        # copy the list into a string and make the string for the inverse generators
        fsub_gen_str = str()
        fsub_invgen_str = str()
        for i in range(len(fsub_gen_first)):
            if i != len(fsub_gen_first)-1:
                fsub_gen_str = fsub_gen_str+fsub_gen_first[i]+","
                fsub_gen_str = fsub_gen_str.swapcase()+fsub_gen_first[i]+","
            else:
                fsub_gen_str = fsub_gen_str+fsub_gen_first[i]
                fsub_gen_str = fsub_gen_str+fsub_gen_first[i].swapcase()
        print(fsub_gen_str)

        # change the appropriate lines in filelines
        fsub_gen_filelines = term_vx.group.file_lines.copy()
        print(fsub_gen_filelines)
        gen_index = str()
        for line in fsub_gen_filelines:
            if "generatorOrder" in line:
                gen_index = fsub_gen_filelines.index(line)

    return "hi"


def process_line(line, flag):
    if flag == "letters":
        line = line[0].split(":=")
        line.pop(0)
        line[0] = line[0].replace(" ", "")
        line[0] = line[0].replace("\n", "")
        line[0] = line[0].replace("[", "")
        line[0] = line[0].replace("]", "")
        line = line[0].split(",")
        line.pop(-1)
    if flag == "equations":
        line = line.replace(" ", "")
        line = line.replace("\n", "")
        line = line.replace("[", "")
        line = line.replace("]", "")
        line = line.split(",")
    return line

#############################################################################
# Testing area
#Gv=group(['a', 'A', 'b', 'B'], ['A', 'a', 'B', 'b'], ['abAB'], "file")
#v=vertex("0", Gv)
#Gfe=group(['a', 'A'], ['A', 'a'], [], "file")
#Gre=group(['b', 'B'], ['B', 'b'], [], "file")
#e=edge("0", v, v, Gfe, Gre)
# i=[]
#gog=graph_of_groups(v, e, i)


# w="bt0abT0bat0"
# gog.cascade(w)
