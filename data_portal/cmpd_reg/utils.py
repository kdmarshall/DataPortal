from matplotlib.colors import ColorConverter
from django.core.files.base import ContentFile
from django.db import connection

from rdkit import Chem 
from rdkit.Chem import Draw
from rdkit.Chem import SaltRemover
from rdkit.Chem import AllChem
from django_rdkit.models import *
from django_rdkit.config import config

from models import (Compound,
					Property,
					Fingerprint)

import StringIO
import psycopg2

"""
This module is for useful cheminformatics related functions
"""

_reactions=None

def get_or_create_compound(smiles, clean=False, inchi_standard=False):
	"""
	Will return the compound object if found
	in the database or will create a new one
	and return it. Also returned a boolean 
	indicating whether or not it was created.
	If more than one Compound object returned,
	an exception is raised.
	If clean=True, input SMILES will be desalted and
	neutralized.
	If inchi_standard=True, INCHI will be used to determine
	if a compounds already exists in the databases. Since
	the standard INCHI canonicalized tautomers, you will see
	collisions between tautomers such as E/Z imines.
	"""
	config.do_chiral_sss = True
	if clean:
		smiles = desalt_neutralize(smiles, return_smiles=True)
	compounds = Compound.objects.filter(molecule__exact=smiles)
	if not compounds or len(compounds) == 0:
		mol = Chem.MolFromSmiles(smiles)
		inchi = Chem.MolToInchi(mol)
		if inchi_standard:
			inchi_check = Compound.objects.filter(inchi=inchi)
			if len(inchi_check) == 0:
				pass
			elif len(inchi_check) == 1:
				return inchi_check[0], False
			else:
				matched_smiles = " AND ".join([inchi_obj.smiles for inchi_obj in inchi_check])
				raise MoleculeMatchException("Inchi from SMILES {} matched an Inchi existing in database with SMILES {}".format(smiles, matched_smiles))
		property = Property(amw=Chem.rdMolDescriptors.CalcExactMolWt(mol),
							hba=Chem.rdMolDescriptors.CalcNumHBA(mol),
							hbd=Chem.rdMolDescriptors.CalcNumHBD(mol),
							tpsa=Chem.rdMolDescriptors.CalcTPSA(mol))
		property.save()

		bfp = Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, 512)
		fingerprint = Fingerprint(bfp=bfp)
		fingerprint.save()

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
		if inchi_standard:
			mol = Chem.MolFromSmiles(smiles)
			inchi = Chem.MolToInchi(mol)
			inchi_check = Compound.objects.filter(inchi=inchi)
			if len(inchi_check) != 1:
				raise MoleculeMatchException("Inchi from SMILES {} did not match one Inchi existing in database".format(smiles))
		return compounds[0], False
	else:
		raise MoleculeMatchException("Exact match search returned multiple molecules")

def desalt_neutralize(smiles, return_smiles=False, remove_only=None):
	"""
	Accepts SMILES and returns a neutralized and desalted
	SMILES string or Mol object. Keyword remove_only
	accepts a string in the format of SMARTS optional list of salt elements
	that should only be removed. For example, remove_only="[Cl,Br]"
	"""
	# TODO spend more time understanding the desalt and neutral process. Make sure
	# it is not changing the molecule too drastically
	desalted = desalt(smiles, remove_only=remove_only)
	desalted_neutralized = neutralize(Chem.MolToSmiles(desalted, isomericSmiles=True))
	if return_smiles:
		return desalted_neutralized[0]
	else:
		return Chem.MolFromSmiles(desalted_neutralized[0])
	

def desalt(smiles, remove_only=None):
	"""
	Accepts SMILES and returns a desalted Mol object.
	Keyword remove_only accepts a string in the format of SMARTS optional list of salt elements
	that should only be removed. For example, remove_only="[Cl,Br]"
	"""
	remover = SaltRemover.SaltRemover(defnData=remove_only) if remove_only else SaltRemover.SaltRemover()
	mol = Chem.MolFromSmiles(smiles)
	resolved = remover.StripMol(mol, dontRemoveEverything=True)
	return resolved

def create_compound_image(compound,smiles=False, size=(300,300), png_path=None):
	"""
	Returns a PIL Image object of compound.
	Writes to PNG file if path is supplied.
	"""
	#TODO instead of smiles keyword flag, just check if compound is instance
	# of basestring or rdkit MOL object
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

def substructure_search(substructure, is_smarts=False, inverted=False):
	"""
	Perform substructure search. If is_smarts is true, uses smarts pattern
	to perform query. If inverted is true, checks to see if compound in
	database is a substructure of the query molecule.
	"""
	config.do_chiral_sss = True
	query = QMOL(Value(substructure)) if is_smarts else substructure
	if inverted:
		compounds = Compound.objects.filter(molecule__issubstruct=query)
	else:
		compounds = Compound.objects.filter(molecule__hassubstruct=query)
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

def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

def neutralize(smiles, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return (Chem.MolToSmiles(mol,True), True)
    else:
        return (smiles, False)

class MoleculeMatchException(Exception):
	"""
	Raised when an exact match search returns an unexpected result.
	"""
	pass

class InvalidSimilarityMethodException(Exception):
	"""
	Raised when an incorrect similarity method is supplied.
	"""
	pass

def test_desalt_neutralize():
	compounds = Compound.objects.all()
	for compound in compounds:
		if '-' in compound.smiles or '+' in compound.smiles or '.' in compound.smiles:
			print "{} -> {}".format(compound.smiles, desalt_neutralize(compound.smiles, return_smiles=True))