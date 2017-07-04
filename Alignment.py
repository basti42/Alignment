
class AlignmentString:
    """
    Beinhaltet moegliche Alignments, inklusive der Gaps und
    matches, nicht matches
    """
    def __init__(self, AlignStringA, AlignStringB):
        self.AlignStringA = AlignStringA
        self.AlignStringB = AlignStringB



class MatrixElement:
    """
    Elemente aus welchen die Matrix aufgebaut werden
    DIAG, TOP, LEFT boolean welche zum Traceback gebraucht werden
    coordinates speichert die Position in der ElementMatrix
    """
    def __init__(self, value, coordinate):
        self.value = value
        self.coordinates = coordinate
        self.DIAG = False
        self.TOP = False
        self.LEFT = False
        self.neighbours = []

    def toString(self):
        """
        return string representation fuer ein MatrixElement
        in der Form: [VALUE, TOP, DIAG, LEFT], Booleans als 0/1 dargestellt
        """
        return "["+str(self.value)+" | " +str(int(self.TOP))+"/" +str(int(self.DIAG))+"/" +str(int(self.LEFT)) + " | " + str(self.coordinates) +"]"



class Alignment:
    """
    Klassen fuer ein Objekt von Alignment, welches sowohl das
    lokale als auch das globale Alignment berechnen und zurueckgeben kann
    Benutzt dynamische Programmierung in einer ElementMatrix
    """
    def __init__(self, StringA, StringB):
        self.StringA = StringA.lower()
        self.StringB = StringB.lower()
        self.Table = [[0 for x in range(len(StringA))] for y in range(len(StringB))]
        self.currMax = 0
        self.indMaxScore = (len(StringA), len(StringB))
        self.match = 0
        self.mismatch = -1
        self.indel = -1
        self.possAlignments = []
        self.buildMatrix()


    def buildMatrix(self):
        """
        Erstellt die Matrix als Tabelle als Grundlage des Alignments
        Values und Pointer fuer den Traceback werden hier berechnet und Erstellt
        # StringA --> spalten
        # StringB --> zeilen
        """
        # erste Zeile und Spalte initialisiert, dann Rest der Matrix auffuellen
        for i in range(len(self.StringA)):
            self.Table[0][i] = MatrixElement(i*self.indel, (0, i))
            if i > 0: self.Table[0][i].neighbours.append((0,i-1)) #nur mgl wenn i > 0 ist
        for j in range(1,len(self.StringB)):
            self.Table[j][0] = MatrixElement(j*self.indel, (j, 0))
            if i > 0: self.Table[j][0].neighbours.append((j-1, 0)) #nur mgl wenn j > 0 ist
        # Rest der Matrix auffuellen
        for row in range(1, len(self.StringB)):
            for col in range(1, len(self.StringA)):
                c = self.match if self.StringA[col-1]==self.StringB[row-1] else self.mismatch
                newval = max(self.Table[row-1][col-1].value + c, self.Table[row-1][col].value + self.indel, self.Table[row][col-1].value + self.indel)
                self.Table[row][col] = MatrixElement(newval, (row, col))
                # wenn newval > currMax, update currMax und die zugehoerigen Indices
                if newval > self.currMax:
                    self.currMax = newval
                    self.indMaxScore = (row, col)
                # Pointer des neuen Elementes setzen
                # zusaetzlich die Liste der Nachbarn mit den Koordinaten
                # der verknuepften Elemente auffuellen, equiv zu Pointer
                # TODO evtl kann Pointer auch ganz weggelassen werden
                if newval == self.Table[row-1][col-1].value + c: # DIAG
                    self.Table[row][col].DIAG = True
                    self.Table[row][col].neighbours.append((row-1, col-1))
                if newval == self.Table[row-1][col].value + self.indel: # TOP
                    self.Table[row][col].TOP = True
                    self.Table[row][col].neighbours.append((row-1, col))
                if newval == self.Table[row][col-1].value + self.indel: #LEFT
                    self.Table[row][col].LEFT = True
                    self.Table[row][col].neighbours.append((row, col-1))


    def displayTable(self):
        """
        Erstellt eine String representation der Matrix mit Hilfe der toString Methode
        aus der MatrixElement Klasse, fuer STDOUT
        """
        for row in self.Table:
            line = ""
            for col in row:
                line += str(col.toString()) + " "
            print(line)

    def globalAlign(self):
        """
        erstellt ein globales Alignment aufgrund der ElementMatrix Table
        """
        # Koordinaten des Alignmentbeginns (index (-1, -1)) in vars entpacken
        row, col = self.Table[-1][-1].coordinates
        for neighbour in self.Table[row][col].neighbours:
            # type(neighbour) = tuple --> entpacken, indizieren
            self.__alignStep__(neighbour[0], neighbour[1], "", "")
        print(row, col)


    def __alignStep__(self, i, j, rowString, colString):
        """
        HELPER Methode zum rekursiven Aufbau des Alignments
        @params row, col : integer, indices zur Bestimmung des aktuellen Elements
        @params rowString, colString : string, aktuelle Alignmentstrings
        """
        # auf ENDE der rückführung checken
        if (i == 0 and j == 0):
            self.possAlignments.append(AlignmentString(rowString, colString))
            #self.alignB.append(colString)

        else:
            # nur wenn Zeile noch nicht fertig / 0 ist
            if i > 0:
                # wenn top, dann gap zum colString und rekursion
                if self.Table[i][j].TOP:
                    rowString = self.StringA[i-1] + rowString
                    colString = "-" + colString
                    self.__alignStep__(i-1, j, rowString, colString)

            # nur wenn Spalte noch nicht fertig / 0 ist
            if j > 0:
                # wenn Left, dann gap zum rowString und rekursion
                if self.Table[i][j].LEFT:
                    rowString = "-" + rowString
                    colString = self.StringB[j-1] + colString
                    self.__alignStep__(i,j-1, rowString, colString)

            # diag nur wenn diag noch mgl ist.
            if i > 0 and j > 0:
                # wenn diag dann alignen und rekursion
                if self.Table[i][j].DIAG:
                    rowString = self.StringA[i-1] + rowString
                    colString = self.StringB[j-1] + colString
                    self.__alignStep__(i-1, j-1, rowString, colString)


# --------------------------- Testing -------------------------------------------
align = Alignment("TAGAG", "CGAGA")
align.displayTable()
align.globalAlign()



"""
TODO
Traceback
set currMax und passende Position properly
global und local Alignment Methoden
"""
