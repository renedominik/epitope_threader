# /////////////////////////////////////////////////////////////////////////
# //
# //  This file is part of the EpitopeThreader (ET) program
# //
# //  Copyright (C) 2013 by Rene Staritzbichler
# //  renedominik@yahoo.de
# //
# //  Dez 16, 2013
# /////////////////////////////////////////////////////////////////////////

#!/usr/bin/python

import sys, os


#####  ADJUST #####

exec_path = ""
#exec_path = "PATH/epitope_threader/"

##### ADJUST END #####



if len( sys.argv) < 6:
    print "\n\nUSAGE: \n\n", sys.argv[0], "  IN_MODELS  IN_EPITOPES  EPITOPE_LENGTH  OUT_PATH  OUT_SCAN_SCORES\n\n" #  (optional:PREFIX) \n\n"
    print "brief description of arguments:\n\n"
    print "IN_MODELS:  \nfile containing list of models without .pdb ending, one per line, (PATH only if necessary) e.g.:\n"
    print " PATH/1afo\n PATH/1zif\n\n"
    print "IN_EPITOPES:  \nfile containing list of epitopes to be scanned, one epitope per line, file needs to start with number of epitopes contained, e.g.:\n"
    print " 3\n ABCDEFGH\n IJKLMNOP\n QRSTUVXY\n\n"  # CHANGE !!!!
    print "EPITOPE_LENGTH:  \nnumber of amino acids of the peptide in models, has to match length of epitope sequences in list (8 in previous example)\n\n"
    print "OUT_PATH:  \npath for output files\n\n"
    print "OUT_SCAN_SCORES:  \nfile for output of scores, the actual result of algorithm\n\n"
#    print "PREFIX: \noptional flag for distinguishing different runs\n\n"
    exit(1)

model_file = sys.argv[1]
epitope_file = sys.argv[2]
epitope_length = sys.argv[3]
out_path = sys.argv[4]
out_scores = sys.argv[5]

if out_path[-1] != "/":
    out_path += "/"


et = exec_path + "exe/epitope_threader.exe  "
top = exec_path + "toppar/top_all27_prot_lipid.rtf "
par = exec_path + "toppar/par_all27_prot_lipid.prm "
default_aa = out_path + "default_amino_acids.et "



prefix = ""

if len( sys.argv) > 6:
    prefix = sys.argv[6]
    print "prefix set to:", prefix

if not os.path.exists( out_path):
    print "create output directory:", out_path
    os.makedirs( out_path)

print "\nextract default amino acids from CHARMM topology and parameter files:"
cmd = et + " -default_amino_acids " + par + top + default_aa
print cmd
sys.stdout.flush()
os.system( cmd)

matrix_file = out_path + out_scores + "_scoring_matrices.txt"

w_scores = open( matrix_file, 'w')


print "\nmain part 1) build score matrices:"
with open( model_file) as f:
    for m in f:
        m = m.strip()
        if len(m) == 0:
            continue
        mol = out_path + prefix
        ind = m.find('/')
        if ind > 0:
            mol += m[ind+1:]
        else:
            mol += m

        # file names
        inpdb = m + ".pdb"
        cleanpdb = mol + ".cleaned.pdb"
        out_et = mol + ".all.et"
        out_ca = mol + "_ca.et"
        out_sorted = mol + "_sorted.et"
        psf = mol + ".psf"
        out_score = mol + ".scm"

        print "\nresort structure file such that epitope is chain A, others B, C, ..."
        print "remark: the length of the epitope needs to be passed (adjust in top section of script)"
        cmd = et + " -sort_chains " + inpdb + " " + epitope_length + " " + out_path   # write chains 
        print cmd
        sys.stdout.flush()
        os.system( cmd )                                                                                                                             
        print "\nwrite command file for calling psfgen tool from NAMD to clean PDB file:"
        # clean pdb with psfgen tool from NAMD
        w = open( psf, 'w') 
        w.write( "/Users/rene/NAMD_2.9_MacOSX-x86_64-multicore/psfgen << ENDMOL \n")
        w.write( "topology " + top + "\n")
        w.write( "pdbalias residue HIS HSE \n")
        w.write( "pdbalias atom ILE CD1 CD \n")
        # find chain pdbs and assign segments
        
        list = ("A", "B", "C", "D", "E", "F", "G", "H")
        check = True
        count = 0
        while( check):
            chain = list[count]
            file = mol + "_" + chain + ".pdb"

            if os.path.exists( file):
                w.write( "segment " + chain + " {pdb " + file + "}\n") # adjust
                w.write( "coordpdb " + file + " " + chain + "\n")
                count += 1
            else:
                check = False

#            w.write( "segment B {pdb " + mol + "_b.pdb} \n")
#            w.write( "coordpdb " + mol + "_b.pdb B \n")
  
        w.write( "guesscoord \n")
        w.write( "writepdb " + cleanpdb + "\n")
        w.write( "ENDMOL \n")
        w.close()
        print "found", (count-1), "chains"

        print "\nchange permissions of command file to executable and call it:"
        cmd = "chmod u+x " + psf + "; ./" + psf
        print cmd
        sys.stdout.flush()
        os.system( cmd )

        print "\nconvert PDB into internal format containing force field information:"
        cmd = et + " -pdb2mol " + par + top + " " + cleanpdb + " " + out_et
        print cmd
        sys.stdout.flush()
        os.system( cmd )

        print "\nmove all atoms of amino acid to its C alpha position (centroid)" 
        cmd = et + " -ca_centroid " + out_et + " " + out_ca
        print cmd
        sys.stdout.flush()
        os.system( cmd )

        print "\nbuild scoring matrix:"
        cmd = et + " -score " + out_ca + " A " + default_aa + " " + out_score 
        print cmd
        sys.stdout.flush()
        os.system( cmd )

        w_scores.write( out_score + "\n")

w_scores.close()



print "\n\nmain part 2) scan epitopes:"
cmd = et + " -scan " + matrix_file + " " + epitope_file + " " + out_scores
print cmd
sys.stdout.flush()
os.system( cmd )
