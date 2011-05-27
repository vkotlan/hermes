# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 15:48:28 2011

@author: vasek
"""
import re 
vstupniSoubor = 'kotouc_mesh.mesh'

Pole=['mag','temp','elast']
""" definice hranic, ktere maji zustat pro jednotliva pole """
Boundaries=range(len(Pole))
Boundaries[Pole.index('mag')]=[150,155,144,143]
Boundaries[Pole.index('temp')]=[155,151,132,129,137,136,134,135,138,139,142,14,8,126,128,127,101,102,77,78,91,92,93,94,95,96,
           86,85,97,98,105,116,113,112,108,106,117,118,121,122,125,71,35,32,148,147,67,65,64,63,60,59,56,55,52,51,22,24,30,33,
           42,21,69,68,31,34,50,49,48,47,46,45,12,11,74,73]
Boundaries[Pole.index('elast')]=[152,4,154,158,153,159,5,6,3,2,161,157,160,156,1]
""" definice oblasti, ktere maji zustat pro jednotliva pole """
Elements=range(len(Pole))
Elements[Pole.index('mag')]=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
         29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
         61,62,63,64,65,66,67]
Elements[Pole.index('temp')]=[0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,67]
Elements[Pole.index('elast')]=[5,67]


Radky=[]
cisloRadku=0
stopRadek=1000000

soubor=file(vstupniSoubor,'r')
for radek in soubor:  
    cisloRadku=cisloRadku+1
    Radky.append(radek)
    if re.search("curves =",radek):
        stopRadek=cisloRadku
soubor.close() 

print "Puvodni soubor ma {0:d} radku ".format(cisloRadku)

for p in range(0,len(Pole)): 
    vystupniSoubor= 'kotouc_mesh_' + str(Pole[p]) + '.mesh'
    soubor_new=file(vystupniSoubor,'w')
    novychRadku=0
    pocetZmen=0
    Zapis=[]
    print "Probiha redukce site pro {0:s} pole: ".format(Pole[p])
    for r in range(0,len(Radky)): 
        cisloRadku=r+1
        radek=Radky[r]        
        if cisloRadku>stopRadek:
            Zapis.append(radek)
            novychRadku=novychRadku+1
            continue
    
           
        if re.search("  { \d+, \d+, \d+, \d+ },\s",radek):
            # """ normalni elementy """
            for i in Elements[p]:   
                if re.search("  { \d+, \d+, \d+, " + str(i) + " },\s",radek):
                    zapis=radek    
                    break
                else: 
                    zapis = "  {},\n"
        elif re.search("  { \d+, \d+, \d+, \d+ }\s",radek):
            # """ posledni element """
            for i in Elements[p]:   
                if re.search("  { \d+, \d+, \d+, " + str(i) + " }\s",radek):
                    zapis=radek    
                    break
                else: 
                    zapis = "  {}\n"
        elif re.search("  { \d+, \d+, \d+ },\s",radek):
            # """ normalni hranice """
            for i in Boundaries[p]:
                if re.search("  { \d+, \d+, " + str(i) + " },\s",radek):
                    zapis=radek 
                    break
                else: 
                    zapis = ""
        elif re.search("  { \d+, \d+, \d+ }\s",radek):
            # """ posledni hranice """
            for i in Boundaries[p]:
                if re.search("  { \d+, \d+, " + str(i) + " }\s",radek):
                    zapis=radek 
                    break
                else: 
                    zapis = ""
                
            if zapis == "":
                # """ musim prepsat predchozi radek """
                predchRadek=Zapis[novychRadku-1]
                # """ musim odstranit carku
                opravenyRadek=predchRadek[0:-2]+"\n"
                #print predchRadek
                #print opravenyRadek
                Zapis[novychRadku-1]=opravenyRadek
        elif re.search("  { \d+, \d+, -\d+ },\s",radek):
            for i in Boundaries[p]:
                if re.search("  { \d+, \d+, -" + str(i-1) + " },\s",radek):
                    zapis=radek 
                    break
                else: 
                    zapis = ""
        elif re.search("  { \d+, \d+, -\d+ }\s",radek):
            for i in Boundaries[p]:
                if re.search("  { \d+, \d+, -" + str(i-1) + " }\s",radek):
                    zapis=radek 
                    break
                else: 
                    zapis = ""
                
            if zapis == "":
                # """ musim prepsat predchozi radek """
                predchRadek=Zapis[novychRadku-1]
                # """ musim odstranit carku
                opravenyRadek=predchRadek[0:-2]+"\n"
                #print predchRadek
                #print opravenyRadek
                Zapis[novychRadku-1]=opravenyRadek
        else: 
            zapis=radek
                
        if zapis!="":
            Zapis.append(zapis)
            novychRadku=novychRadku+1
            if zapis==(("  {},\n") or ("  {}\n")):
                pocetZmen=pocetZmen+1
        else:
            pocetZmen=pocetZmen+1
    
    for r in range(0,len(Zapis)):
        soubor_new.write(Zapis[r])
   
    soubor_new.close()
    
    print "\tNovych radku: {0:d}".format(novychRadku)
    print "\tSmazano radku: {0:d}".format(cisloRadku-novychRadku)
    print "\tZmeneno radku: {0:d}".format(pocetZmen-(cisloRadku-novychRadku))



    
    

