from django.shortcuts import render
from django.http import HttpResponse
from django.core.files.base import ContentFile
from django_rdkit.config import config

from utils import create_compound_image, get_or_create_compound, create_image_tag, substructure_search, create_substruct_image

import StringIO

def index(request):
	return render(request, 'cmpd_reg/index.html', {})


def creg_index(request):
	config.do_chiral_sss = True
	if request.method == 'POST':
		smiles = request.POST.get("smiles", None)
		try:
			result_compound, is_new = get_or_create_compound(smiles)
		except:
			return render(request, 'cmpd_reg/cmpd-reg.html', {
				'smiles_error': "Error! Could not convert SMILES '{}' to molecule object".format(smiles),
			})
		searched_image = create_compound_image(smiles, smiles=True)
		if searched_image:
			img_tag = create_image_tag(searched_image)
			return render(request, 'cmpd_reg/cmpd-reg.html', {
				'searched_smiles': smiles,
				'searched_compound_image': img_tag,
				'returned_compound': str(result_compound.id),
				'is_new': str(is_new),
			})
		else:
			return render(request, 'cmpd_reg/cmpd-reg.html', {
				'smiles_error': "Error! Could not convert SMILES '{}' to image".format(smiles),
			})

	return render(request, 'cmpd_reg/cmpd-reg.html', {})


def substructure(request):
	config.do_chiral_sss = True
	if request.method == 'POST':
		smiles = request.POST.get("smiles", None)
		searched_image = create_compound_image(smiles, smiles=True)
		if searched_image:
			searched_img_tag = create_image_tag(searched_image)
			matched_substructures = substructure_search(smiles)
			if not matched_substructures:
				return render(request, 'cmpd_reg/substructure.html', {
					'searched_smiles': smiles,
					'searched_compound_image': searched_img_tag,
					'no_matches': True,
				})
			matched_output = []
			for substruct in matched_substructures:
				cmpd_obj = {}
				cmpd_obj['id'] = str(substruct.id)
				cmpd_obj['smiles'] = substruct.smiles
				pil_image = create_substruct_image(substruct.molecule, smiles)
				cmpd_obj['image'] = create_image_tag(pil_image)
				matched_output.append(cmpd_obj)
			return render(request, 'cmpd_reg/substructure.html', {
				'searched_smiles': smiles,
				'searched_compound_image': searched_img_tag,
				'matched_substructures': matched_output,
			})
		else:
			return render(request, 'cmpd_reg/substructure.html', {
				'smiles_error': "Error! Could not convert SMILES '{}' to image".format(smiles),
			})
	return render(request, 'cmpd_reg/substructure.html', {})

def similarity(request):
	return render(request, 'cmpd_reg/similarity.html', {})

def bulk_loader(request):
	return render(request, 'cmpd_reg/bulk_loader.html', {})
