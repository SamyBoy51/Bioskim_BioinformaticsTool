from collections import Counter

#Diccionarios que le robé a dani pq no quiero copiar esta mierda (y lo modifique pq estaba una mierda tbmn)
DiccionarioUwU = {
    'AUA': ('Ile','I'), 'AUC': ('Ile','I'), 'AUU': ('Ile','I'), 'AUG': ('Met','M'), 'ACA': ('Thr','T'), 'ACC': ('Thr','T'), 'ACG': ('Thr','T'), 'ACU': ('Thr','T'),
    'AAC': ('Asn','N'), 'AAU': ('Asn','N'), 'AAA': ('Lys','K'), 'AAG': ('Lys','K'), 'AGC': ('Ser','S'), 'AGU': ('Ser','S'), 'AGA': ('Arg','R'), 'AGG': ('Arg','R'),
    'CUA': ('Leu','L'), 'CUC': ('Leu','L'), 'CUG': ('Leu','L'), 'CUU': ('Leu','L'), 'CCA': ('Pro','P'), 'CCC': ('Pro','P'), 'CCG': ('Pro','P'), 'CCU': ('Pro','P'),
    'CAC': ('His','H'), 'CAU': ('His','H'), 'CAA': ('Gln','Q'), 'CAG': ('Gln','Q'), 'CGA': ('Arg','R'), 'CGC': ('Arg','R'), 'CGG': ('Arg','R'), 'CGU': ('Arg','R'),
    'GUA': ('Val','V'), 'GUC': ('Val','V'), 'GUG': ('Val','V'), 'GUU': ('Val','V'), 'GCA': ('Ala','A'), 'GCC': ('Ala','A'), 'GCG': ('Ala','A'), 'GCU': ('Ala','A'),
    'GAC': ('Asp','D'), 'GAU': ('Asp','D'), 'GAA': ('Glu','E'), 'GAG': ('Glu','E'), 'GGA': ('Gly','G'), 'GGC': ('Gly','G'), 'GGG': ('Gly','G'), 'GGU': ('Gly','G'),
    'UCA': ('Ser','S'), 'UCC': ('Ser','S'), 'UCG': ('Ser','S'), 'UCU': ('Ser','S'), 'UUC': ('Phe','F'), 'UUU': ('Phe','F'), 'UUA': ('Leu','L'), 'UUG': ('Leu','L'),
    'UAC': ('Tyr','Y'), 'UAU': ('Tyr','Y'), 'UGC': ('Cys','C'), 'UGU': ('Cys','C'), 'UGG': ('Trp','W'),
    'UAA': ('STOP','*'), 'UAG': ('STOP','*'), 'UGA': ('STOP','*')
}


def readFile(file):
    try:
        #deletion part of the non important elements
        with open(file, 'r') as data: #data is the name of the file variable for the reading process
            lineas = data.readlines()

            seq = "".join(lineas).replace("\n", "").upper()
            if seq.startswith(">"):
                seq = "".join(lineas[1:]).replace("\n", "").upper()
            else:
                seq = "".join(lineas).replace("\n", "").upper()
            qtyGuanina = seq.count("G")
            qtyCitosina = seq.count("C")
            total = len(seq)
            if total ==0: return 0
            GCpercentage = ((qtyCitosina + qtyGuanina)/total) *100
            return round(GCpercentage, 2), total

    except Exception as e:
        return f"Error: {str(e)}"


def analysis(file, rango):
    prot1 = []
    prot3 = []
    try:
        inicio,fin = rango.split('-')
        inicio, fin = int(inicio), int(fin)

        #nonCoding strand
        ncStrand = []
        #mArn
        mArn = []
        with open(file, 'r') as data: #data is the name of the file variable for the reading process
            lineas = data.readlines()
            seq = "".join(lineas[1:]).replace("\n", "").upper()
            strand = seq[inicio-1: fin]

            if not strand:
                return "Rango fuera de los limites de la secuencia."
            #Hier as kontinue die oprazionen.
            #1. Cadena 3' a 5'.
            for x in strand:
                if x == 'A':
                    ncStrand.append('T')
                elif x =='G':
                    ncStrand.append("C")
                elif x=='C':
                    ncStrand.append("G")
                elif x=='T':
                    ncStrand.append("A")

            #2. mArn
            for y in ncStrand:
                if y =="C":
                    mArn.append('G')
                elif y =="G":
                    mArn.append('C')
                elif y =="T":
                    mArn.append('A')
                elif y =="A":
                    mArn.append('U')



            res1 = f"Secuencia exraída desde {inicio} a {fin}: \n {strand}"
            res2 = mArn
            res3 = ncStrand
            res4 = "".join(strand)
            #Remember if u place return a,b. youll have an array of the different str values. such that outppu = [a,b].

            # 3. mArn to Protein:
            StrmArn = "".join(mArn)
            for i in range(0,len(StrmArn)// 3 * 3,3):
                codon = StrmArn[i:i+3]
                e3,e1 = DiccionarioUwU.get(codon,('???','?'))

                if e3 == 'STOP':
                    break
                prot3.append(e3)
                prot1.append(e1)

            totAA = len(prot1)
            cuenta = Counter(prot1)
            stats = {
                aa: {
                    "cantidad": cant,
                    "porcentaje": round((cant/ totAA)*100,2),
                    "": "\n"
                } for aa, cant in cuenta.items()
            }if totAA > 0 else{}


            return res1, res2, res3, "".join(prot1), stats, res4
    except ValueError:
        return "Error: utilice el formato de rango indicado."
    except FileNotFoundError:
        return "Error: no se encontró el archivo de ADN."
    except Exception as e:
        return f"Error inesperado: {str(e)}"

def Transcripcion(ADN):
    prot1 = []
    mArn = ADN.replace("T","U")
    StrmArn = "".join(mArn)
    for i in range(0, len(StrmArn) // 3 * 3, 3):
        codon = mArn[i:i + 3]
        e3, e1 = DiccionarioUwU.get(codon, ('???', '?'))

        if e3 == 'STOP':
            break
        prot1.append(e1)
        res1 = "".join(prot1)
    return res1