# Program to read lines and write (output) selected ones
# Python 2.7.3
# if use, thanks to: 
#                   Francisco Paz-Chinchon. 
#                   francisco at dfte.ufrn.br
#                   DFTE-UFRN, Natal, Brazil.
#
import os
import numpy
from sys import exit

# Explanation:
#!/bin/sh
#baja todas las curvas en el directorio seleccionado, desde el ftp
# -np : no guarda los directorios superiores
# -m : todos los archivos del directorio senalado
# -nH : no pone una carpea superior con el nombre del ftp
# ---cut-dirs=4 : la carpeta esta a 4 niveles del pto superior del ftp

#
# Write procedure
def Write(out,result):
  head = "#!/bin/sh \n # chmod a+x wget_file.bash \n # ./wget_file.bash"
  # wget command
  comm = "wget -m -np -nH --cut-dirs=4 --reject '*_slc.fits' ftp://anonymous:@archive.stsci.edu/pub/kepler/lightcurves/"
  with open(out,'w') as f:
    f.writelines(head+'\n')
    for i in result:
      f.writelines(comm+'\n')
  f.close()

# Open list 
filein = "kic_list.txt"
datalin = [line.strip() for line in open(filein, 'r')]	#datalin is an list of strings

# Decide and fill list
temp = []
for ID in datalin:
  if len(ID)==8:
    temp.append("0"+ID[:3]+"/"+"0"+ID)
  elif len(ID)==7:
    temp.append("00"+ID[:2]+"/"+"00"+ID)
  elif len(ID)==6:
    temp.append("000"+ID[:1]+"/"+"000"+ID)
  else:
    print "ID length not contempled", len(ID)
    exit(0)

Write("Wget_SolarThin.bash",temp)

exit(0)

#..................................................
# DEPRECATED:


dao_lin = []
fin = []
Diff = 100
# Open Daospec files and search for the closest line value
path = "/home/fj/gislana"
for (path, dirs, files) in os.walk('.'):
  for file in files:
    if ".daospec" in file:
      # Charge values into 2 lists of floats
      count = 0
      fin_aux = []
      fin_aux = ["#"+file,'#Line'+'\t'+'DAO_line'+'\t'+'DAO_depth']
      
      f =  open(file, 'r')
      temp = [line.strip() for line in f]	#datalin is an list of strings
      cc = 0
      for j in temp:
	if cc > 1:
	  dao_lin.append(float(j.split()[0]))	#spec line
	  dao_dep.append(float(j.split()[2]))	#line depth
	cc = cc+1
      f.close()
      
      #Calculate less difference & store results
      for line in datalin:
	data = float(line)
	aux_dep = "NULL"
	aux_lin = "NULL"
	c2 = 0
	for comp_lin in dao_lin:
	  if Diff(data,comp_lin) < diff: 
	    diff = Diff(data,comp_lin)
	    aux_lin = comp_lin
	    aux_dep = dao_dep[c2]
	    	    
	  c2= c2 + 1 ###
	wr_aux = line+'\t'+str(aux_lin)+'\t'+str(aux_dep)
	fin_aux.append(wr_aux)
	diff = 100	#restore default value to next line iteration
      
      #Write each output file
      Write(file,fin_aux)
      dao_lin = []
      dao_dep = []
      
