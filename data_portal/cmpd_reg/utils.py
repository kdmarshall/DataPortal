from matplotlib.colors import ColorConverter
from django.core.files.base import ContentFile
from django.db import connection

from rdkit import Chem 
from rdkit.Chem import Draw
from models import (Compound,
					Property,
					Fingerprint)
from django_rdkit.models import *
from django_rdkit.config import config
import StringIO
import psycopg2

"""
This module is for useful cheminformatics related functions
"""

def get_or_create_compound(smiles):
	"""
	Will return the compound object if found
	in the database or will create a new one
	and return it. Also returned a boolean 
	indicating whether or not it was created.
	If more than one Compound object returned,
	returns None since is an error.
	"""
	config.do_chiral_sss = True
	compounds = Compound.objects.filter(molecule__exact=smiles)
	if not compounds or len(compounds) == 0:
		mol = Chem.MolFromSmiles(smiles)

		property = Property(amw=Chem.rdMolDescriptors.CalcExactMolWt(mol),
							hba=Chem.rdMolDescriptors.CalcNumHBA(mol),
							hbd=Chem.rdMolDescriptors.CalcNumHBD(mol),
							tpsa=Chem.rdMolDescriptors.CalcTPSA(mol))
		property.save()

		bfp = Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, 512)
		fingerprint = Fingerprint(bfp=bfp)
		fingerprint.save()

		inchi = Chem.MolToInchi(mol)
		inchi_key = Chem.InchiToInchiKey(inchi)
		compound = Compound(smiles=smiles,
							molecule=mol,
							inchi=inchi,
							inchi_key=inchi_key,
							ctab=Chem.MolToMolBlock(mol),
							fingerprint=fingerprint,
							property=property)
		compound.save()
		return compound, True
	elif len(compounds) == 1:
		return compounds[0], False
	else:
		raise MultipleMoleculeMatchException("Exact match search returned multiple molecules")

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


def write_to_sdf(mol_id_list):
	"""
	Pass the function a list of mol ids
	and the function returns an open StringIO object.
	"""
	mol_list = [Compound.objects.get(pk=int(id)) for id in mol_id_list]
	sio = StringIO.StringIO()
	for mol in mol_list:
		sio.write(mol.ctab+'\n')
		sio.write('>  <{property}>\n{value}\n\n'.format(property='Compound_ID',value=str(mol.id)))
		sio.write('>  <{property}>\n{value}\n\n'.format(property='AMW',value=str(mol.property.amw)))
		sio.write('>  <{property}>\n{value}\n\n'.format(property='HBA',value=str(mol.property.hba)))
		sio.write('>  <{property}>\n{value}\n\n'.format(property='HBD',value=str(mol.property.hbd)))
		sio.write('>  <{property}>\n{value}\n\n'.format(property='TPSA',value=str(mol.property.tpsa)))
		sio.write('$$$$\n')
	return sio


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

def substructure_search(substructure, is_smarts=False):
	config.do_chiral_sss = True
	compounds = Compound.objects.filter(molecule__hassubstruct=substructure)
	if len(compounds) == 0:
		return None
	else:
		return compounds

def similarity_search(smiles, method='tanimoto', threshold=0.5):
	"""
	Perform similarity search. Method default is tanimoto with
	a threshold of 0.5
	"""
	config.do_chiral_sss = True
	config.tanimoto_threshold = threshold
	method_options = ('tanimoto', 'dice')
	if method not in method_options:
		raise InvalidSimilarityMethodException("Method must be either tanimoto or dice")
	if method == 'tanimoto':
		compounds = Compound.objects.filter(fingerprint__bfp__tanimoto=MORGANBV_FP(Value(smiles)))
	else:
		compounds = Compound.objects.filter(fingerprint__bfp__dice=MORGANBV_FP(Value(smiles)))
	if len(compounds) == 0:
		return None
	else:
		return compounds


class MultipleMoleculeMatchException(Exception):
	"""
	Raised when an exact match search returns more
	than one molecule.
	"""
	pass

class InvalidSimilarityMethodException(Exception):
	"""
	Raised when an incorrect similarity method is supplied.
	"""
	pass