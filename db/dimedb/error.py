inchi = "asdasd"

from subprocess import check_output

def get_smiles(inchi):
    val = check_output('openbabel.obabel -iinchi -:"%s" -osmi' %inchi, shell=True).decode().strip()
    if val != "":
        return val
    else:
        return None

get_smiles(inchi)