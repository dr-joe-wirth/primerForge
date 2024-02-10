""" provides the nearest neighbor thermodynamics data
"""
from Bio.SeqUtils.MeltingTemp import DNA_NN4
import copy

# start with the original table
NN_TABLE = copy.deepcopy(DNA_NN4)

# add data from https://doi.org/10.1021/bi9825091
NN_TABLE.update({'AA/TA': (1.2,1.7),
                 'CA/GA': (-0.9,-4.2),
                 'GA/CA': (-2.9,-9.8),
                 'TA/AA': (4.7,12.9),
                 'AC/TC': (0.0,-4.4),
                 'CC/GC': (-1.5,-7.2),
                 'GC/CC': (3.6,8.9),
                 'TC/AC': (6.1,16.4),
                 'AG/TG': (-3.1,-9.5),
                 'CG/GG': (-4.9,-15.3),
                 'GG/CG': (-6.0,-15.8),
                 'TG/AG': (1.6,3.6),
                 'AT/TT': (-2.7,-10.8),
                 'CT/GT': (-5.0,-15.8),
                 'GT/CT': (-2.2,-8.4),
                 'TT/AT': (0.2,-1.5)})

# add data from https://doi.org/10.1093/nar/26.11.2694
NN_TABLE.update({'AC/TT': (0.7,0.2),
                 'AT/TC': (-1.2,-6.2),
                 'CC/GT': (-0.8,-4.5),
                 'CT/GC': (-1.5,-6.1),
                 'GC/CT': (2.3,5.4),
                 'GT/CC': (5.2,13.5),
                 'TC/AT': (1.2,0.7),
                 'TT/AC': (1.0,0.7)})

# add data from https://doi.org/10.1021/bi962590c
NN_TABLE.update({'AG/TT': (1.0,0.9),
                 'AT/TG': (-2.5,-8.3),
                 'CG/GT': (-4.1,-11.7),
                 'CT/GG': (-2.8,-8.0),
                 'GG/CT': (3.3, 10.4),
                 'GG/TT': (5.8,16.3),
                 'GT/CG': (-4.4,-12.3),
                 'GT/TG': (4.1,9.5),
                 'TG/AT': (-0.1,-1.7),
                 'TG/GT': (-1.4,-6.2),
                 'TT/AG': (-1.3,-5.3)})

# add data from https://doi.org/10.1021/bi9803729
NN_TABLE.update({'AA/TC': (2.3,4.6),
                 'AC/TA': (5.3,14.6),
                 'CA/GC': (1.9,3.7),
                 'CC/GA': (0.6,-0.6),
                 'GA/CC': (5.2,14.2),
                 'GC/CA': (-0.7,-3.8),
                 'TA/AC': (3.4,8.0),
                 'TC/AA': (7.6,20.2)})

# add data from https://doi.org/10.1021/bi9724873
NN_TABLE.update({'AA/TG': (-0.6,-2.3),
                 'AG/TA': (-0.7,-2.3),
                 'CA/GG': (-0.7,-2.3),
                 'CG/GA': (-4.0,-13.2),
                 'GA/CG': (-0.6,-1.0),
                 'GG/CA': (0.5,3.2),
                 'TA/AG': (0.7,0.7),
                 'TG/AA': (3.0,7.4)})
