"""The program is used for noncompliance with requirements in pdb format file a and store the data in b.
    
    synx a: input file
    synx b: output dealed file
"""
def clean_file(a,b):
    """The program is used for noncompliance with requirements in pdb format file a,
    and store the data in b.
    
    synx a: input file
    synx b: output dealed file
    """
    with open(a) as f_in: 
        with open(b, "w") as f_out: 
            for line in f_in: 
                if line[:4] == "ATOM": 
                    f_out.write(line)
 
#print celan_file.__doc__
#print 'nihao'

