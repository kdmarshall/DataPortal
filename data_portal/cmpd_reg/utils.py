from matplotlib.colors import ColorConverter
from django.core.files.base import ContentFile
from django.db import connection

from rdkit import Chem 
from rdkit.Chem import Draw
from models import (Compound,
					Property,
					Fingerprint)
from django_rdkit.models import (MOL_TO_CTAB,
								 MOL_INCHI,
								 MOL_INCHIKEY)
import StringIO
import psycopg2

"""
This module is for useful cheminformatics related functions
"""

def bulk_load():
	pass

def get_or_create_compound(smiles):
	"""
	Will return the compound object if found
	in the database or will create a new one
	and return it. Also returned a boolean 
	indicating whether or not it was created.
	If more than one Compound object returned,
	returns None since is an error.
	"""
	compounds = Compound.objects.filter(molecule__exact=smiles)
	if not compounds or len(compounds) == 0:
		mol = Chem.MolFromSmiles(smiles)
		compound = Compound(smiles=smiles,
							molecule=mol,
							inchi=MOL_INCHI(mol),
							inchi_key=MOL_INCHIKEY(mol),
							ctab=Chem.MolToMolBlock(mol))
		compound.save()
		return compound, True
	elif len(compounds) == 1:
		return compounds[0], False
	else:
		return None, False

def desalt_neutralize(smiles, return_mol=False):
	"""
	Accepts SMILES and returns a neutralized 
	SMILES string or mol object.
	"""
	pass

def create_compound_image(compound,smiles=False, size=(300,300), png_path=None):
	"""
	Returns a PIL Image object of compound.
	Writes to PNG file if path is supplied.
	"""
	if smiles:
		mol = Chem.MolFromSmiles(compound)
		# if mol == None:
		# 	return None
	else:
		mol = compound
	img = Draw.MolToImage(mol, size=size, wedgeBonds=True, kekulize=True)
	if not img:
		return None
	if png_path:
		img.save(png_path)
	return img

def create_image_tag(pil):
	output = StringIO.StringIO()
	pil.save(output, "PNG")
	contents = output.getvalue().encode("base64")
	output.close()
	img_tag = '<img src="data:image/png;base64,{}" />'.format(contents)
	return img_tag

def create_substruct_image(mol, substruct_smarts, size=(300,300), highlight_color='aqua', png_path=None):
	"""
	Returns a PIL Image object of compound with
	matching substructure highlighted. Writes to PNG 
	file if path is supplied.
	"""
	pattern = Chem.MolFromSmarts(substruct_smarts)
	matching = mol.GetSubstructMatch(pattern)
	try:
		color = ColorConverter().to_rgb(highlight_color)
	except Exception:
		print "Cannot find color {}. Setting to default aqua.".format(highlight_color)
		color = ColorConverter().to_rgb('aqua')
	color = tuple(int(i) for i in color)
	img = Draw.MolToImage(mol, highlightAtoms=matching, highlightColor=color, size=size, wedgeBonds=True, kekulize=True)
	if png_path:
		img.save(png_path)
	return img


def write_to_sdf(mol_list, file_path):
	pass


def set_chiral_flag():
    SET_FLAG = 'set rdkit.do_chiral_sss=true'
    try:
        cursor = connection.cursor()
        cursor.execute(SET_FLAG)
        print "**** CHIRAL FLAG SET ****"

    except psycopg2.DatabaseError, e:
        print 'Error %s' % e

def remove_chiral_flag():
    REMOVE_FLAG = 'set rdkit.do_chiral_sss=false'
    try:
        cursor = connection.cursor()
        cursor.execute(REMOVE_FLAG)
        print "**** CHIRAL FLAG REMOVED ****"

    except psycopg2.DatabaseError, e:
        print 'Error %s' % e

def substructure_search(smiles_substructure):
	compounds = Compound.objects.filter(molecule__hassubstruct=smiles_substructure)
	if len(compounds) == 0:
		return None
	else:
		return compounds