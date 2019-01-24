# This program can convert an outdat file of MiStImm into three separate csv files. 
# run example:
# $ python Outdat2CSV.py outdat-1518016285-1956-respERS

import sys
IN_name=sys.argv[1]

OUT_no_elements=open(IN_name+'-no_elements.csv', 'w')
print(IN_name+'-no_elements.csv')
OUT_no_elements.write('t,self_cells,foreign_antigens,B_cells,antibodies,Thelpers,danger_signals_and_ILs,bone_marrow_cells\n')


OUT_Thelpers=open(IN_name+'-Thelpers.csv', 'w')
print(IN_name+'-Thelpers.csv')
OUT_Thelpers.write('x,y,peptide_distance,t_birth,t_death\n')

OUT_Bcells=open(IN_name+'-Bcells.csv', 'w')
print(IN_name+'-Bcells.csv')
OUT_Bcells.write('x,y,r,maturity,stress,nab1,t_birth,t_death\n')

IN=open(IN_name)
for line in IN:
    sl=line.strip().split()
    if line[:2] == 't=':
        OUT_no_elements.write(','.join([sl[i] for i in range(1,16,2)])+'\n')
    elif line[:8] == 'Thelper:':
        OUT_Thelpers.write(','.join(sl[1:])+'\n')
    elif line[:6] == 'black:':
        OUT_Bcells.write(','.join(sl[3:])+'\n')


