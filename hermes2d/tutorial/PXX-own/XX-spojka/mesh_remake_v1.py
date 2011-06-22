# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 15:48:28 2011

@author: vasek
"""
import re 
vstupniSoubor = 'mesh_file_curves.mesh'
ProblemNazev = 'spojka'
Pole=['mag','temp','elast']
""" definice hranic, ktere maji zustat pro jednotliva pole """
Boundaries=range(len(Pole))
Boundaries[Pole.index('mag')]=[17,18,19,20]
Boundaries[Pole.index('temp')]=[13,14,15,16]
Boundaries[Pole.index('elast')]=[1,2,3,4,5,6,7,8]
""" definice oblasti, ktere maji zustat pro jednotliva pole """
Elements=range(len(Pole))
Elements[Pole.index('mag')]=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]
Elements[Pole.index('temp')]=[1,2,3,4,5,6,7,8,9,10,11,12,13]
Elements[Pole.index('elast')]=[4]


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
print "Stop radek = {0:d} ".format(stopRadek)

print "Puvodni soubor ma {0:d} radku ".format(cisloRadku)

for p in range(0,len(Pole)): 
    vystupniSoubor= str(ProblemNazev) + '_mesh_' + str(Pole[p]) + '.mesh'
    soubor_new=file(vystupniSoubor,'w')
    novychRadku=0
    pocetZmen=0
    Zapis=[]
    PoleEdgesBody=[]
    pocetEdgu = 0 
    pocetCurves = 0
    print "Probiha redukce site pro {0:s} pole: ".format(Pole[p])
    for r in range(0,len(Radky)): 
        cisloRadku=r+1
        radek=Radky[r]        
        if cisloRadku>=stopRadek:
#            Zapis.append(radek)
#            novychRadku=novychRadku+1
            # prohledat pole edges a smazat curves pokud se body neschodují
            if re.search("  { \d+, \d+, .*},\s",radek):
                for i in range(0,len(PoleEdgesBody)):
                    if re.search(PoleEdgesBody[i] + ",.*},\s",radek):
                        zapis=radek
                        pocetCurves = pocetCurves +1
                        break                                    
                    else:
                        zapis= ""
            elif re.search("  { \d+, \d+, .*}\s",radek):
                for i in range(0,len(PoleEdgesBody)):
                    if re.search(PoleEdgesBody[i] + ",.*}\s",radek):
                        zapis=radek
                        pocetCurves = pocetCurves +1
                        break                                    
                    else:
                        zapis= ""
                        
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
                    pocetEdgu = pocetEdgu+1
                    body = radek.split(",")
                    PoleEdgesBody.append(str(body[0]) + "," + str(body[1]))
                    break
                else: 
                    zapis = ""
        elif re.search("  { \d+, \d+, \d+ }\s",radek):
            # """ posledni hranice """
            for i in Boundaries[p]:
                if re.search("  { \d+, \d+, " + str(i) + " }\s",radek):
                    zapis=radek 
                    pocetEdgu = pocetEdgu+1
                    body = radek.split(",")
                    PoleEdgesBody.append(str(body[0]) + "," + str(body[1]))
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
    print "\tpocet edgu: {0:d}".format(pocetEdgu)
    print "\tpocet curves: {0:d}".format(pocetCurves)
    
    
    print "Vypis Bodu:"
    for i in range(0,len(PoleEdgesBody)): 
        print PoleEdgesBody[i]



    
    

