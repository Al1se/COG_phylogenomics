"""
Module for manipulation with the PDB files
------- Version: 1.6
        1.0  * Basic implementation of <Atom>, <Ligand> and <PDB_file> classes
        1.5  * <PDB_file> can now be obtained from the CIF-formatted file
        1.6  * <auth_asym_id> is now correctly used as chain_id, not <label_asym_id>        
"""
import sys, os, re

class Atom:
    """
    Class describes a single atom in a PDB or CIF format
    """    
    def __init__(self, atom_string, column_names = None):
        if column_names == None:
            self.init_from_PDB_format(atom_string)
        else:
            try:
                self.init_from_CIF_format(atom_string, column_names)
            except:
                print("Failed to read the following string:\n'%s'" % atom_string)
                raise
        self.pdb_string = atom_string

    def init_from_CIF_format(self, atom_string, column_names):
        """
        Method requires an <atom_string> in CIF format to define an <Atom>. Example:
        ATOM   1      O  "O5'" . C   A  1  1    ? 183.60600 272.57800 306.73800 1.000 55.01000 ? 1    C   L5 "O5'" 1 
        The following attributes will be created from it:
        * record_type: 'HETATM'       _atom_site.group_PDB
        * name       : 'O'            _atom_site.label_atom_id
        * res_name   : 'HOH'          _atom_site.label_comp_id        
        * chain_id   : 'L5'           _atom_site.auth_asym_id (FIXED in v. 1.6: WRONG: * chain_id : 'A'_atom_site.label_asym_id) 
        * res_number : 171 (integer)  _atom_site.label_seq_id
        * x          : 38.133 (float) _atom_site.Cartn_x
        * y          : 12.471 (float) _atom_site.Cartn_y
        * z          : 18.173 (float) _atom_site.Cartn_z
        * element    : 'O'            _atom_site.type_symbol
        * pdb_string : <atom_string>  
        """        
        fields = re.split("\s+", atom_string)
        for i in range(len(column_names)):
            if column_names[i] == "group_PDB":
                self.record_type = fields[i]
            if column_names[i] == "label_atom_id":
                self.name = fields[i]
            if column_names[i] == "label_comp_id":
                self.res_name = fields[i]
            #if column_names[i] == "label_asym_id": #FIX v. 1.6: wrong column was taken
            #    self.chain_id = fields[i]
            if column_names[i] == "auth_asym_id":
                self.chain_id = fields[i]
            if column_names[i] == "label_seq_id":
                try:
                    self.res_number = int(fields[i])
                except ValueError:
                    self.res_number = 0
            if column_names[i] == "Cartn_x":
                self.x = float(fields[i])
            if column_names[i] == "Cartn_y":
                self.y = float(fields[i])
            if column_names[i] == "Cartn_z":
                self.z = float(fields[i])
            if column_names[i] == "type_symbol":
                self.element = fields[i]

    def init_from_PDB_format(self, atom_string):
        """
        Method requires an <atom_string> in PDB format to define an <Atom>. Example:
        HETATM 1305  O   HOH A 171      38.133  12.471  18.173  1.00 11.12           O
        The following attributes will be created from it:
        * record_type: 'HETATM'
        * name       : 'O'
        * res_name   : 'HOH'
        * chain_id   : 'A'
        * res_number : 171 (integer)
        * x          : 38.133 (float)
        * y          : 12.471 (float)
        * z          : 18.173 (float)
        * element    : 'O'
        * pdb_string : <atom_string>
        """        
        self.record_type = atom_string[0:6].strip()
        #serial = int(atom_string[6:11].strip())
        self.name = atom_string[12:16].strip()
        #alt_loc = atom_string[16].strip()
        self.res_name = atom_string[17:20].strip()
        self.chain_id = atom_string[21].strip()
        #insert_cod = atom_string[26].strip()
        self.res_number = int(atom_string[22:26].strip())
        self.x = float(atom_string[30:38].strip())
        self.y = float(atom_string[38:46].strip())
        self.z = float(atom_string[46:54].strip())
        #self.occupancy = float(atom_string[54:60].strip())
        #self.temp_factor = float(atom_string[60:66].strip())
        self.element = atom_string[76:78].strip()
        #charge = atom_string[78:79].strip()        

    def print_data(self):
        print ("Data for the atom:")
        print ("record_type: '%s'" % self.record_type)
        print ("name       : '%s'" % self.name)
        print ("res_name   : '%s'" % self.res_name)
        print ("chain_id   : '%s'" % self.chain_id)
        print ("res_number : '%s'" % self.res_number)
        print ("x          : '%s'" % self.x)
        print ("y          : '%s'" % self.x)
        print ("z          : '%s'" % self.x)
        print ("element    : '%s'" % self.element)

class Ligand:
    """
    Class defines a single ligand instance (non-protein, non-RNA, non-DNA and non-water
    set of atoms)
    """
    def __init__(self, ligand_id, name):
        """
        Method creates the following attributes of the <Ligand>:
        * ligand_id : <ligand_id>, taken from HETNAM record of a PDB file
        * name      : <name>, taken from HETNAM record of a PDB file
        * formula   : None (should be a chemical formula from FORMUL record)
        * syn       : empty list, will be filled with synonymous names (HETSYN record)
        * instances : empty dict of lists of <Atom> instances. Key is formed like: 2CGP_CMP_621 
        """
        self.ligand_id = ligand_id
        self.name = name
        self.formula = None
        self.syn = list()       
        self.instances = dict()
        
def work_out_fields(fields, column_names, log):
    if len(fields) != len(column_names):
        log.write("[WARNING] Length of space-separated fields in a string #%s and in column list is different!\n" % n)
        log.write(column_names)
        log.write(fields)
        return None
    ligand = None
    ligand_id = None
    ligand_type = None
    ligand_name = None
    ligand_synonyms = None
    ligand_formula = None
    for i in range(len(column_names)):
        if column_names[i] == "id":
            ligand_id = fields[i]
        if column_names[i] == "type":
            ligand_type = fields[i]
        if column_names[i] == "name":
            ligand_name = fields[i]
        if column_names[i] == "pdbx_synonyms":
            ligand_synonyms = fields[i]
        if column_names[i] == "formula":
            ligand_formula = fields[i]
    if ligand_type != "L-peptide linking":
        ligand = Ligand(ligand_id, ligand_name)
        ligand.formula = ligand_formula
        ligand.syn.append(ligand_synonyms)
    return ligand

def work_out_cif_loop(ligand_loop_data, column_names, log):
    ligands = list()
    fields = list()
    multiline_started = False
    for string in ligand_loop_data:
        if (string[0] == ";"): # Multiline property mark
            if multiline_started == False:
                multiline_started = True
                fields.append("")
            string = string.strip(";")
            fields[-1] += string
            if len(string) == 0: # This is the end of a multiline property
                multiline_started = False
        else:
            if len(fields) == len(column_names):
                #print("Working out fields:")
                #print(fields)
                ligand = work_out_fields(fields, column_names, log)
                if ligand != None:
                    ligands.append(ligand)
                fields = list()
            if len(fields) > len(column_names):
                print("WARNING: more fields than column_names!")
                print(fields)
                print(column_names)
                return list()

            reading_long_field = False
            long_field_symbol = None
            prev_symbol = None
            for s in string:
                if (s == "'") and (reading_long_field == False): # Starting of the space-containing field
                    reading_long_field = True
                    long_field_symbol = "'"
                    fields.append("")
                    continue
                if (s == '"') and (reading_long_field == False): # Starting of the space-containing field
                    reading_long_field = True
                    long_field_symbol = '"'
                    fields.append("")
                    continue
                if reading_long_field == True:                   
                    if (s == long_field_symbol): # Ending of the space-containing field
                        reading_long_field = False
                    else:
                        fields[-1] += s
                else:
                    if (prev_symbol == " ") or (prev_symbol == None):
                        if (s != " "):
                            fields.append("")
                            fields[-1] += s
                    else:
                        if (s != " "):
                            fields[-1] += s
                prev_symbol = s

    if len(fields) == len(column_names):
        #print("Working out fields:")
        #print(fields)
        ligand = work_out_fields(fields, column_names, log)
        if ligand != None:
            ligands.append(ligand)

    return ligands

class PDB_file:
    """
    Class describes a single PDB file
    """
    def __init__(self, filename, log = None, cif_format = False):
        """
        Defines templates for the following attributes of a class:
        * header     : '' (short description of the PDB-file from the HEADER record)
        * title      : '' (extended description of the PDB-file from the TITLE record)
        * pdbid      : '' (PDB ID from the HEADER record)
        * org        : None (name of the organism)
        * org_taxid  : None (taxonomy ID of the organism)
        * ligands    : empty dict of <Ligand> objects (key is the <Ligand>.ligand_id)
        * atoms      : empty list of all protein atoms (<Atom> objects)
        * misc_atoms : empty list of all non-protein and non-ligand atoms (<Atom> objects)
        Runs <self.obtain_data> method to fill these attributes
        #FIX: v. 1.5 - if <cif_format> is True, runs another method: <self.obtain_data_from_cif> 
        If <log_file> is not None, all warnings are appended to it, otherwise - to STDOUT
        """
        self.header = ""
        self.title = ""
        self.pdbid = ""
        self.org = None
        self.org_taxid = None
        self.ligands = dict()
        self.atoms = list()
        self.misc_atoms = list()
       
        if log == None:
            log = sys.stdout
        if cif_format == False:
            self.obtain_data(filename, log)
        else:
            self.obtain_data_from_cif(filename, log)

    def print_ligand_data(self):
        for ligand_id in self.ligands.keys():
            ligand = self.ligands[ligand_id]
            print ("%s\t'%s'" % (ligand_id, ligand.name))
            print ("Synonyms:")
            for s in ligand.syn:
                print s
            print ("------------------------------------")             

    def obtain_data_from_cif(self, filename, log):
        """
        Actually defines all attributes of an object with the data provided in a
        CIF-file <filename>. Writes warnings to the <log>: could be STDOUT or a file
        """
        pdb_file = open(filename, "r")
        curr_table_columns = list()
        ligand_loop_data = list()
        column_names = list()
        table_columns_started = False
        atom_table_started = False
        ligand_table_started = False
        n = 0
        for string in pdb_file:
            string = string.strip()
            n += 1
            if len(string) == 0:
                continue
            global_string_parts = re.split("\s+", string)
            if (string[0] == "_") and (len(global_string_parts) == 2):
                if global_string_parts[0] == "_entry.id":
                    self.pdbid = global_string_parts[1]
            if string == "loop_": # Starting of a table
                curr_table_columns = list()
                ligand_loop_data = list() # all strings reffering to this loop
                table_columns_started = True
                continue
            if table_columns_started == True:
                if string[0] == "_":
                    parts = string.split(".")
                    curr_table_columns.append((parts[0], parts[1]))
                else:
                    table_columns_started = False
                    column_names = list()
                    table_type = None
                    for column_values in curr_table_columns:
                        column_names.append(column_values[1])
                        if table_type == None:
                            table_type = column_values[0]
                        else:
                            if column_values[0] != table_type:
                                print ("FATAL ERROR: Table contains different values: '%s' and '%s' (string #%s)" % (table_type, column_values[0], n))
                                sys.exit()
                    if table_type == "_atom_site": # This is a table descriding atom positions
                        atom_table_started = True
                    if table_type == "_chem_comp": # This is a table with ligands
                        ligand_table_started = True
            if (atom_table_started == True) and (not string[0] in ["H", "A"]): # ATOM and HETATM parts are finished
                atom_table_started = False
            if (ligand_table_started == True) and (string[0] == "#"): # ligand data finished
                curr_ligands = work_out_cif_loop(ligand_loop_data, column_names, log)
                for curr_ligand in curr_ligands:
                    if not curr_ligand.ligand_id in self.ligands:
                        self.ligands[curr_ligand.ligand_id] = curr_ligand
                    else:
                        print("WARNING: ligand '%s' occured more than once in file '%s'!" % (ligand_id, self.pdbid))                
                ligand_table_started = False          
  
            if atom_table_started == True:
                fields = re.split("\s+", string)
                if len(fields) != len(column_names):
                    log.write("Length of space-separated fields in a string #%s and in column list is different!" % n)
                    log.write(column_names)
                    log.write(fields)
                    sys.exit()
                atom_type = None
                PDB_type = None
                name = None
                for i in range(len(column_names)):                   
                    if column_names[i] == "label_comp_id":
                        name = fields[i]
                    if column_names[i] == "group_PDB":
                        PDB_type = fields[i]
                if PDB_type == "ATOM":
                    atom_type = "PROTEIN"    
                    if (name == "DA") or (name == "DT") or (name == "DC") or (name == "DG"):
                        atom_type = "DNA"
                    if (name == "A") or (name == "T") or (name == "C") or (name == "G"):
                        atom_type = "RNA"  
                if PDB_type == "HETATM":
                    if name == "MSE": # Selenocystein
                        atom_type = "PROTEIN"
                    elif (name == "HOH") or (name == "DOD"):
                        atom_type = "WATER"
                    else:
                        atom_type = "LIGAND"
                if atom_type != None:
                    curr_atom = Atom(string, column_names)
                    if atom_type == "PROTEIN":
                        self.atoms.append(curr_atom)
                    if (atom_type == "DNA") or (atom_type == "RNA") or (atom_type == "WATER"):
                        self.misc_atoms.append(curr_atom)
                    if atom_type == "LIGAND":
                        curr_key = "%s_%s_%s" % (self.pdbid, curr_atom.res_name, curr_atom.res_number)
                        if not curr_atom.res_name in self.ligands: # Such ligand type exists
                            log.write("[WARNING] in pdb %s: ligand '%s' found in HETATOM, but not before!\n" % (self.pdbid, curr_atom.res_name))
                            self.ligands[curr_atom.res_name] = Ligand(curr_atom.res_name, "Unknown")
                        if not curr_key in self.ligands[curr_atom.res_name].instances: # This particular ligand exists
                            self.ligands[curr_atom.res_name].instances[curr_key] = list()      
                        self.ligands[curr_atom.res_name].instances[curr_key].append(curr_atom)
                else:
                    log.write("[WARNING] in pdb %s: This atom was not considered as atom of supported types!\n")
                    log.write("%s\n" % string)
            if ligand_table_started == True:
                ligand_loop_data.append(string)                
        pdb_file.close()

    def obtain_data(self, filename, log):
        """
        Actually defines all attributes of an object with the data provided in a
        PDB-file <filename>. Writes warnings to the <log>: could be STDOUT or a file
        """
        pdb_file = open(filename, "r")
        for string in pdb_file:
            string = string.strip()
            if len(string) == 0:
                continue

            if string[0:6] == "HEADER":
                #HEADER    HYDROLASE                               20-SEP-14   4RDY
                fields = re.split(" +", string[10:].strip())
                self.pdbid = fields.pop(-1)
                fields.pop(-1)
                self.header = " ".join(fields)

            if string[0:5] == "TITLE":
                #TITLE     CRYSTAL STRUCTURE OF THE DELTA-PYRROLINE-5-CARBOXYLATE DEHYDROGENASE  
                #TITLE    2 FROM MYCOBACTERIUM TUBERCULOSIS                                      
                self.title += " " + string[10:].strip()

            if string[0:6] == "HETNAM":
                #HETNAM     1ET 6,6'-{[5-(AMINOMETHYL)BENZENE-1,3-DIYL]                          
                #HETNAM   2 1ET  DIETHANE-2,1-DIYL}BIS(4-METHYLPYRIDIN-2-AMINE)                  
                #HETNAM   3 1ET                                                                  

                #HETNAM     KCX LYSINE NZ-CARBOXYLIC ACID
                fields = re.split(" ", string[10:].strip(), 1)
                if not fields[0] in self.ligands:
                    self.ligands[fields[0]] = Ligand(fields[0], fields[1])
                else:
                    if len(fields) == 2:
                        self.ligands[fields[0]].name += fields[1]

            if string[0:6] == "HETSYN":
                #HETSYN     3M5 N-(3-OXODECANOYL)-L-HOMOSERINE LACTONE
                fields = re.split(" +", string[10:].strip(), 1)
                #print fields
                ligand_id = fields[0]
                if ligand_id in self.ligands:
                    parts = fields[1].split(" ;")
                    #print parts
                    if string[9] != " ": # Additional string                        
                        self.ligands[ligand_id].syn[-1] += parts[0]
                        if len(parts) > 1: 
                            self.ligands[ligand_id].syn.extend(parts[1:])
                    else:                       
                        self.ligands[ligand_id].syn.extend(parts)
                else:
                    log.write("[WARNING] in pdb %s: ligand '%s' (%s) has synonym name, but no normal name!\n" % (
                              self.pdbid, ligand_id, fields[1]))
                          
            if string[0:6] == "FORMUL":
                #FORMUL   1  KCX    2(C7 H14 N2 O4)
                fields = re.split(" +", string[12:].strip(), 1)
                ligand_id = fields[0]
                if len(fields) == 1:
                    log.write("[WARNING] in pdb %s: formula is missing! String: '%s'\n" % (self.pdbid, string))
                else:
                    formula_search = re.search("^\**\d*\(*([^\)]+)\)*$", fields[1])
                    if formula_search == None:
                        log.write("[WARNING] in pdb %s: formula was not detected! Formula field: '%s', string: '%s'\n" % (self.pdbid, fields[1], string))
                    else:
                        formula = formula_search.group(1)
                        if ligand_id in self.ligands:
                            if self.ligands[ligand_id].formula == None: 
                                self.ligands[ligand_id].formula = formula
                            else:
                                log.write("[WARNING] in pdb %s: ligand '%s' has several formula: '%s' and '%s'\n" % (                      
                                          self.pdbid, ligand_id, self.ligands[ligand_id].formula, formula))                          
                        else:
                            if (ligand_id != "HOH") and (ligand_id != "DOD"): # Water never has name
                                log.write("[WARNING] in pdb %s: ligand '%s' (%s) has formula, but no normal name!\n" % (
                                          self.pdbid, ligand_id, formula))     

            atom_type = None # By default, this is not a coordinate string
            if (string[0:4] == "ATOM"):
                #ATOM     35  P    DG L   2      35.975 -35.601  -5.885  1.00 30.00           P
                atom_type = "PROTEIN"
                name = string[17:20].strip()
                if (name == "DA") or (name == "DT") or (name == "DC") or (name == "DG"):
                    atom_type = "DNA"
                if (name == "A") or (name == "T") or (name == "C") or (name == "G"):
                    atom_type = "RNA"  
            if string[0:6] == "HETATM":
                #HETATM 1305  O   HOH A 171      38.133  12.471  18.173  1.00 11.12           O
                if string[17:20] == "MSE": # Selenocystein
                    atom_type = "PROTEIN"
                elif (string[17:20] == "HOH") or (string[17:20] == "DOD"):
                    atom_type = "WATER"
                else:
                    atom_type = "LIGAND"

            if ((string[0:6] == "HETATM") or (string[0:4] == "ATOM")) and (atom_type == None):
                log.write("[WARNING] in pdb %s: This atom was not considered as atom of supported types!\n")
                log.write("%s\n" % string) 

            if atom_type != None:
                curr_atom = Atom(string)
                if atom_type == "PROTEIN":
                    self.atoms.append(curr_atom)
                if (atom_type == "DNA") or (atom_type == "RNA") or (atom_type == "WATER"):
                    self.misc_atoms.append(curr_atom)
                if atom_type == "LIGAND":
                    curr_key = "%s_%s_%s" % (self.pdbid, curr_atom.res_name, curr_atom.res_number)
                    if not curr_atom.res_name in self.ligands: # Such ligand type exists
                        log.write("[WARNING] in pdb %s: ligand '%s' found in HETATOM, but not before!\n" % (self.pdbid, curr_atom.res_name))
                        self.ligands[curr_atom.res_name] = Ligand(curr_atom.res_name, "Unknown")
                    if not curr_key in self.ligands[curr_atom.res_name].instances: # This particular ligand exists
                        self.ligands[curr_atom.res_name].instances[curr_key] = list()      
                    self.ligands[curr_atom.res_name].instances[curr_key].append(curr_atom)
                                
        pdb_file.close()
        self.title = self.title.strip()
        #self.print_ligand_data()