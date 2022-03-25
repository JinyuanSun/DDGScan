#! /usr/bin/env python
import numpy as np

class Atom:
    #  1 -  6        Record name   "ATOM  "
    #  7 - 11        Integer       serial       Atom  serial number.
    # 13 - 16        Atom          name         Atom name.
    # 17             Character     altLoc       Alternate location indicator.
    # 18 - 20        Residue name  resName      Residue name.
    # 22             Character     chainID      Chain identifier.
    # 23 - 26        Integer       resSeq       Residue sequence number.
    # 27             AChar         iCode        Code for insertion of residues.
    # 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    # 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    # 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    # 55 - 60        Real(6.2)     occupancy    Occupancy.
    # 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    # 77 - 78        LString(2)    element      Element symbol, right-justified.
    # 79 - 80        LString(2)    charge       Charge  on the atom.
    def __init__(self, record):
        self.record = self.checkrecord(record)
        self.serial = int(self.record[6:11])
        self.name = self.record[12:16]
        self.altLoc = self.record[16]
        self.resName = self.record[17:20]
        self.chainID = self.record[21]
        self.resSeq = int(self.record[22:26])
        self.iCode = self.record[26]
        self.x = self.record[30:38]
        self.y = self.record[38:46]
        self.z = self.record[46:54]
        self.occupancy = self.record[54:60]
        self.tempFactor = self.record[60:66]
        self.element = self.record[76:78]
        self.charge = self.record[78:80]

    def get_coord(self):
        return np.array([self.x, self.y, self.z], dtype=float)

    def checkrecord(self, record):
        if len(record) == 80:
            return record
        if len(record) > 80:
            return record[:80]
        if len(record) < 80:
            pad = " " * 80
            return record + pad[len(record):]


class AminoAcid:
    def __init__(self):
        self.resnum: int
        self.resname: str
        self.atomdict = {}


    def readline(self, line):
        self.resnum = line[17:20]


if __name__ == '__main__':
    record = "ATOM     10  CA  THR A   2      28.716  29.848  34.951  1.00 13.11           C              "
    CA2 = Atom(record)
    print(CA2.get_coord())
