import numpy as np
import re
for line in open(input("Enter the protein(Any PDB file): ")):
    list = line.split()
    if line.startswith('HELIX'):
        print (line[0:40],line[74:76])
    if line.startswith('SHEET'):
        print (line[0:40],line[74:76])        
print('\n Helix residues are generated') 
for line in open (input("Type the PDB_ID again here for Chain_*: ")):
    list = line.split()
    if list[0] =='ATOM':
         if list[2] =='CA':
           type = list [2]
           if type == 'CA':
            Chain = list[4]
             if Chain == '*':
              residue_no = int(list[4])
              position = list[6:9]
           if residue_no == "J1":
             C1 = line[31:55]
             P1 = re.sub("\s+", ",", C1.strip())
           if residue_no == "J1":
             C2 = line[31:55]
             P2 = re.sub("\s+", ",", C2.strip())
           if residue_no == "J1":
             C3 = line[31:55]
             P3 = re.sub("\s+", ",", C3.strip())
           if residue_no == "J1":
             C4 = line[31:55]
             P4 = re.sub("\s+", ",", C4.strip())
           if residue_no == "J1":
             C5 = line[31:55]
             P5 = re.sub("\s+", ",", C5.strip())
           if residue_no == "J1":
             C6 = line[31:55]
             P6 = re.sub("\s+", ",", C6.strip())
           if residue_no == "J1":
             C7 = line[31:55]
             P7 = re.sub("\s+", ",", C7.strip())
           if residue_no == "J1":
             C8 =  line[31:55]
             P8 = re.sub("\s+", ",", C8.strip())
           if residue_no == "J1":
             C9 = line[31:55]
             P9 = re.sub("\s+", ",", C9.strip())
           if residue_no == "J1":
             C10 = line[31:55]
             P10 = re.sub("\s+", ",", C10.strip())
           if residue_no == "J1":
             C11 = line[31:55]
             P11 = re.sub("\s+", ",", C11.strip())
           if residue_no == "J1":
             C12 = line[31:55]
             P12 = re.sub("\s+", ",", C12.strip())
           if residue_no == "J1":
            C13 = line[31:55]
            P13 = re.sub("\s+", ",", C13.strip())
           if residue_no == "J1":
            C14 = line[31:55]
            P14 = re.sub("\s+", ",", C14.strip())
print('#########################################\n')
print('Angle calculations: \n')
P1 = eval(P1)
P2 = eval(P2)
P3 = eval(P3)
P4 = eval(P4)
P5 = eval(P5)
P6 = eval(P6)
P7 = eval(P7)
P8 = eval(P8)
P9 = eval(P9)
P10 = eval(P10)
P11 = eval(P11)
P12 = eval(P12)
P13 = eval(P13)
P14 = eval(P14)
B1 = np.subtract(P2,P1)
B2 = np.subtract(P3,P2)
B3 = np.subtract(P4,P3)
B4 = np.subtract(P5,P4)
B5 = np.subtract(P6,P5)
B6 = np.subtract(P7,P6)
B7 = np.subtract(P8,P7)
B8 = np.subtract(P9,P8)
B9 = np.subtract(P10,P9)
B10 = np.subtract(P11,P10)
B11 = np.subtract(P12,P11)
B12 = np.subtract(P13,P12)
B13 = np.subtract(P14,P13)
B14 = np.subtract(P3,P2)
B15 = np.subtract(P4,P3)
B16 = np.subtract(P5,P4)
B17 = np.subtract(P6,P5)
B18 = np.subtract(P7,P6)
B19 = np.subtract(P8,P7)
B20 = np.subtract(P9,P8)
B21 = np.subtract(P10,P9)
B22 = np.subtract(P11,P10)
B23 = np.subtract(P12,P11)
B24 = np.subtract(P13,P12)
N1 = np.cross(B1,B2)
N2 = np.cross(B2,B3)
N3 = np.cross(B3,B4)
N4 = np.cross(B4,B5)
N5 = np.cross(B5,B6)
N6 = np.cross(B6,B7)
N7 = np.cross(B7,B8)
N8 = np.cross(B8,B9)
N9 = np.cross(B9,B10)
N10 = np.cross(B10,B11)
N11 = np.cross(B11,B12)
N12 = np.cross(B12,B13)
N13 = np.cross(B14,B15)
N14 = np.cross(B15,B16)
N15 = np.cross(B16,B17)
N16 = np.cross(B17,B18)
N17 = np.cross(B18,B19)
N18 = np.cross(B19,B20)
N19 = np.cross(B20,B21)
N20 = np.cross(B21,B22)
N21 = np.cross(B22,B23)
N22 = np.cross(B23,B24)
B_2 = np.sqrt((B2*B2).sum())
B_4 = np.sqrt((B4*B4).sum())
B_6 = np.sqrt((B6*B6).sum())
B_8 = np.sqrt((B8*B8).sum())
B_10 = np.sqrt((B10*B10).sum())
B_12 = np.sqrt((B12*B12).sum())
B_15 = np.sqrt((B15*B15).sum())
B_17 = np.sqrt((B17*B17).sum())
B_19 = np.sqrt((B19*B19).sum())
B_21 = np.sqrt((B21*B21).sum())
B_23 = np.sqrt((B23*B23).sum())
Unit_B2 = B2/B_2
Unit_B4 = B4/B_4
Unit_B6 = B6/B_6
Unit_B8 = B8/B_8
Unit_B10 = B10/B_10
Unit_B12 = B12/B_12
Unit_B15 = B15/B_15
Unit_B17 = B17/B_17
Unit_B19 = B19/B_19
Unit_B21 = B21/B_21
Unit_B23 = B23/B_23
M1 = np.cross(N1,Unit_B2)
M2 = np.cross(N3,Unit_B4)
M3 = np.cross(N5,Unit_B6)
M4 = np.cross(N7,Unit_B8)
M5 = np.cross(N9,Unit_B10)
M6 = np.cross(N11,Unit_B12)
M7 = np.cross(N13,Unit_B15)
M8 = np.cross(N15,Unit_B17)
M9 = np.cross(N17,Unit_B19)
M10 = np.cross(N19,Unit_B21)
M11 = np.cross(N21,Unit_B23)
X1 = np.dot(N1,N2)
X2 = np.dot(N3,N4)
X3 = np.dot(N5,N6)
X4 = np.dot(N7,N8)
X5 = np.dot(N9,N10)
X6 = np.dot(N11,N12)
X7 = np.dot(N13,N14)
X8 = np.dot(N15,N16)
X9 = np.dot(N17,N18)
X10 = np.dot(N19,N20)
X11 = np.dot(N21,N22)
Y1 = np.dot(M1,N2)
Y2 = np.dot(M2,N4)
Y3 = np.dot(M3,N6)
Y4 = np.dot(M4,N8)
Y5 = np.dot(M5,N10)
Y6 = np.dot(M6,N12)
Y7 = np.dot(M7,N14)
Y8 = np.dot(M8,N16)
Y9 = np.dot(M9,N18)
Y10 = np.dot(M10,N20)
Y11 = np.dot(M11,N22)
dihedral_1 = np.arctan2(Y1,X1)
dihedral_2 = np.arctan2(Y2,X2)
dihedral_3 = np.arctan2(Y3,X3)
dihedral_4 = np.arctan2(Y4,X4)
dihedral_5 = np.arctan2(Y5,X5)
dihedral_6 = np.arctan2(Y6,X6)
dihedral_7 = np.arctan2(Y7,X7)
dihedral_8 = np.arctan2(Y8,X8)
dihedral_9 = np.arctan2(Y9,X9)
dihedral_10 = np.arctan2(Y10,X10)
dihedral_11 = np.arctan2(Y11,X11)
print(("In_radians_1: "),(dihedral_1), ("In_radians_2: "),(dihedral_2))
print(("In_radians_3: "),(dihedral_3), ("In_radians_4: "),(dihedral_4))
print(("In_radians_5: "),(dihedral_5), ("In_radians_6: "),(dihedral_6))
In_degrees_1 = dihedral_1*180/np.pi 
In_degrees_2 = dihedral_2*180/np.pi
In_degrees_3 = dihedral_3*180/np.pi
In_degrees_4 = dihedral_4*180/np.pi
In_degrees_5 = dihedral_5*180/np.pi
In_degrees_6 = dihedral_6*180/np.pi
In_degrees_7 = dihedral_7*180/np.pi 
In_degrees_8 = dihedral_8*180/np.pi
In_degrees_9 = dihedral_9*180/np.pi
In_degrees_10 = dihedral_10*180/np.pi
In_degrees_11 = dihedral_11*180/np.pi
print('#########################################\n')
print('H_Angle calculations: \n')
print((">> Angle between(H1 - H2)in_degrees_1: "),round(In_degrees_1,2))
print((">> Angle between(H2 - H3)in_degrees_2: "),round(In_degrees_2,2))
print((">> Angle between(H3 - H4)in_degrees_3: "),round(In_degrees_3,2))
print((">> Angle between(H4 - H5)in_degrees_4: "),round(In_degrees_4,2))
print((">> Angle between(H5 - H6)in_degrees_5: "),round(In_degrees_5,2))
print((">> Angle between(H6 - H7)in_degrees_6: "),round(In_degrees_6,2))
print('#########################################\n')
print('L_Angle calculations: \n')
print((">> Angle between(L1 - L2)in_degrees_7: "),round(In_degrees_7,2))
print((">> Angle between(L2 - L3)in_degrees_8: "),round(In_degrees_8,2))
print((">> Angle between(L3 - L4)in_degrees_9: "),round(In_degrees_9,2))
print((">> Angle between(L4 - L5)in_degrees_10: "),round(In_degrees_10,2))
print((">> Angle between(L5 - L6)in_degrees_11: "),round(In_degrees_11,2))
with open("Output.txt",'w+') as f:
 f.write(" H_Angles {}\n\n")
 f.write(" 1.Angle between(H1 - H2)in_degrees_1: {}\n\n".format(In_degrees_1))
 f.write(" 2.Angle between(H2 - H3)in_degrees_2: {}\n\n".format(In_degrees_2))
 f.write(" 3.Angle between(H3 - H4)in_degrees_3: {}\n\n".format(In_degrees_3))
 f.write(" 4.Angle between(H4 - H5)in_degrees_4: {}\n\n".format(In_degrees_4))
 f.write(" 5.Angle between(H5 - H6)in_degrees_5: {}\n\n".format(In_degrees_5))
 f.write(" 6.Angle between(H6 - H7)in_degrees_6: {}\n\n".format(In_degrees_6))
 f.write(" L_Angles {}\n\n")
 f.write(" 7.Angle between(L1 - L2)in_degrees_7: {}\n\n".format(In_degrees_7))
 f.write(" 8.Angle between(L2 - L3)in_degrees_8: {}\n\n".format(In_degrees_8))
 f.write(" 9.Angle between(L3 - L4)in_degrees_9: {}\n\n".format(In_degrees_9))
 f.write(" 10.Angle between(L4 - L5)in_degrees_10: {}\n\n".format(In_degrees_10))
 f.write(" 11.Angle between(L5 - L6)in_degrees_11: {}\n\n".format(In_degrees_11))
 f.close()
