#!/usr/bin/env python
import sys, os, argparse, time
import udav_pdb
sys.path.append("D:\\UdavBackup\\_Complete_genomes\\_scripts")
import udav_soft

#========================================================================================
curr_version = 1.0
parser = argparse.ArgumentParser(description = 
"This script will read all '.ent'/'.pdb' or '.cif' files in a given directory and will \
document all required ligands in required chains. \
Current version is %s" % curr_version 
)
parser.add_argument("-i", help = "Name of the directory with PDB files", required = True, dest = "input_dir")
parser.add_argument("-o", help = "Prefix for the output files", required = True, dest = "output")
parser.add_argument("-l", help = "List of required ligand IDs", required = True, dest = "ligands")
parser.add_argument("-c", help = "List of required chains (containing e.g. '2CGP_A' etc)", required = True, dest = "chains")
myargs = parser.parse_args()

req_ligands = udav_soft.read_ordered_list(myargs.ligands)
(COG_to_chains, chains) = udav_soft.read_two_column_assignment(myargs.chains)
chain_to_COG = dict()
COG_count = dict()
for COG in COG_to_chains.keys():
    for chain in COG_to_chains[COG]:
        if not COG in COG_count:
            COG_count[COG] = 0
        COG_count[COG] += 1
        if not chain in chain_to_COG:
            chain_to_COG[chain] = COG
        else:
            print("WARNING: Chain '%s' is assigned to different COGs, all except for the first assignment are ignored!" % (chain))

pdb_log = open("%s.log" % myargs.output, "w")
pdb_log.write("#This log file was created by <find_ligands_in_PDB.py> script at: %s\n" % time.ctime())
pdb_log.write("#Input directory                       : %s\n" % myargs.input_dir)
pdb_log.write("#File with the list of required ligands: %s\n" % myargs.ligands)
pdb_log.write("#Total ligands in the list             : %i\n" % len(req_ligands))
if len(req_ligands) < 20:
    pdb_log.write("#List of required ligand IDs           : %s\n" % req_ligands)
pdb_log.write("#File with the list of required chains : %s\n" % myargs.chains)
pdb_log.write("#Total chain groups found              : %i\n" % len(COG_to_chains))
pdb_log.write("#-------------------------------------------------------------------------------\n")

n = 0
t_start = time.time()
COG_to_ligands = dict() # Dictionary of dictionaries with all data
file_list = os.listdir(myargs.input_dir)
for curr_file in file_list:
    n += 1
    print ("Working with the file '%s' (%i of %i)" % (curr_file, n, len(file_list)))
    pdb_log.write("Working with the file '%s' (%i of %i)\n" % (curr_file, n, len(file_list)))
    curr_pdb = None
    extension = curr_file.split(".")[-1]
    file_path = os.path.join(myargs.input_dir, curr_file)
    try:
        if extension in ["ent", "pdb"]:
            curr_pdb = udav_pdb.PDB_file(file_path, pdb_log)
        if extension == "cif":
            curr_pdb = udav_pdb.PDB_file(file_path, pdb_log, cif_format = True)
        if curr_pdb != None:
            for ligand_id in curr_pdb.ligands.keys():
                ligand = curr_pdb.ligands[ligand_id]
                for l_instance in ligand.instances.keys():                 
                    pdb_with_chain = "%s_%s" % (curr_pdb.pdbid.lower(), ligand.instances[l_instance][0].chain_id)
                    if pdb_with_chain in chain_to_COG:
                        curr_COG = chain_to_COG[pdb_with_chain]
                        if ligand_id in req_ligands:                    
                            if not curr_COG in COG_to_ligands:
                                COG_to_ligands[curr_COG] = dict()
                            if not ligand_id in COG_to_ligands[curr_COG]:
                                COG_to_ligands[curr_COG][ligand_id] = dict()
                            if not pdb_with_chain in COG_to_ligands[curr_COG][ligand_id]:
                                COG_to_ligands[curr_COG][ligand_id][pdb_with_chain] = 0
                            COG_to_ligands[curr_COG][ligand_id][pdb_with_chain] += 1
        else:
            pdb_log.write("[ERROR] The PDB file '%s' red without exceptions, but returned 'None'" % curr_file)
    except:
        pdb_log.write("[ERROR] Failed to read PDB file '%s'" % curr_file)
t_end = time.time()
time_worked = t_end - t_start
pdb_log.write("#-------------------------------------------------------------------------------\n")
pdb_log.write("DONE in %.1f seconds (or %.1f minutes)!" % (time_worked, time_worked/60))
pdb_log.close()

output = open("%s.ligands.txt" % myargs.output, "w")
output_num = open("%s.ligands.num.txt" % myargs.output, "w")
output.write("#COG\tChain Num.")
output_num.write("#COG\tChain Num.")
for ligand in req_ligands:
    output.write("\t%s" % ligand)
    output_num.write("\t%s" % ligand)
output.write("\n")
output_num.write("\n")

for COG in COG_count.keys():
    output.write("%s\t%i" % (COG, COG_count[COG]))
    output_num.write("%s\t%i" % (COG, COG_count[COG]))
    for ligand in req_ligands:
        chain_data = None
        n_lig = 0
        if COG in COG_to_ligands:
            if ligand in COG_to_ligands[COG]:
                chain_data = ", ".join(COG_to_ligands[COG][ligand].keys())
                chain_data = chain_data.strip()
                n_lig += len(list(COG_to_ligands[COG][ligand].keys()))
        if chain_data == None:
            chain_data = ""
        output.write("\t%s" % chain_data)
        output_num.write("\t%i" % n_lig)
    output.write("\n")
    output_num.write("\n")
output.close()
output_num.close()