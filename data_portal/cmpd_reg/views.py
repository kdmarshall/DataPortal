from django.shortcuts import render
from django.http import HttpResponse
from django.core.files.base import ContentFile
from django.core.serializers.json import DjangoJSONEncoder
from django.utils.safestring import mark_safe

from utils import (create_compound_image,
				   get_or_create_compound,
				   create_image_tag,
				   substructure_search,
				   create_substruct_image,
				   similarity_search,
				   write_to_sdf)

import StringIO
import csv
import json
import time

def index(request):
	return render(request, 'cmpd_reg/index.html', {})


def creg_index(request):
	if request.method == 'POST':
		smiles = request.POST.get("smiles", None)
		inchi = request.POST.get("inchi", False)
		inchi_check = True if inchi == 'on' else False
		if not smiles:
			return render(request, 'cmpd_reg/cmpd-reg.html', {})
		try:
			result_compound, is_new = get_or_create_compound(smiles, clean=True, inchi_standard=inchi_check)
		except Exception, e:
			print str(e)
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
	if request.method == 'POST':
		smiles = request.POST.get("smiles", None)
		smarts = request.POST.get("smarts", None)
		search_type = request.POST.get("search_type", None)
		search_type_options = {
			'Normal': False,
			'Inverted': True
		}
		if smarts:
			print "Using SMARTS pattern: {}".format(smarts)
			matched_substructures = substructure_search(smarts, is_smarts=True, inverted=search_type_options[search_type])
			if not matched_substructures:
				return render(request, 'cmpd_reg/substructure.html', {
					'searched_smiles': smarts,
					'no_matches': True,
				})
			matched_output = []
			id_list = []
			for substruct in matched_substructures:
				cmpd_obj = {}
				cmpd_obj['id'] = str(substruct.id)
				id_list.append(str(substruct.id))
				cmpd_obj['smiles'] = substruct.smiles
				cmpd_obj['inchi'] = substruct.inchi
				cmpd_obj['property'] = substruct.property
				pil_image = create_substruct_image(substruct.molecule, smarts)
				cmpd_obj['image'] = create_image_tag(pil_image)
				matched_output.append(cmpd_obj)
			return render(request, 'cmpd_reg/substructure.html', {
				'searched_smiles': smarts,
				'matched_substructures': matched_output,
				'id_list': mark_safe(json.dumps(id_list, cls=DjangoJSONEncoder)),
			})
		searched_image = create_compound_image(smiles, smiles=True)
		if searched_image:
			searched_img_tag = create_image_tag(searched_image)
			matched_substructures = substructure_search(smiles, inverted=search_type_options[search_type])
			if not matched_substructures:
				return render(request, 'cmpd_reg/substructure.html', {
					'searched_smiles': smiles,
					'searched_compound_image': searched_img_tag,
					'no_matches': True,
				})
			matched_output = []
			id_list = []
			for substruct in matched_substructures:
				cmpd_obj = {}
				cmpd_obj['id'] = str(substruct.id)
				id_list.append(str(substruct.id))
				cmpd_obj['smiles'] = substruct.smiles
				cmpd_obj['inchi'] = substruct.inchi
				cmpd_obj['property'] = substruct.property
				pil_image = create_substruct_image(substruct.molecule, smiles)
				cmpd_obj['image'] = create_image_tag(pil_image)
				matched_output.append(cmpd_obj)
			return render(request, 'cmpd_reg/substructure.html', {
				'searched_smiles': smiles,
				'searched_compound_image': searched_img_tag,
				'matched_substructures': matched_output,
				'id_list': mark_safe(json.dumps(id_list, cls=DjangoJSONEncoder)),
			})
		else:
			return render(request, 'cmpd_reg/substructure.html', {
				'smiles_error': "Error! Could not convert SMILES '{}' to image".format(smiles),
			})
	return render(request, 'cmpd_reg/substructure.html', {})

def similarity(request):
	if request.method == 'POST':
		smiles = request.POST.get("smiles", None)
		method = request.POST.get('similarity_method', None)
		threshold = request.POST.get('threshold', None)
		searched_image = create_compound_image(smiles, smiles=True)
		if searched_image:
			searched_img_tag = create_image_tag(searched_image)
			returned_compounds = similarity_search(smiles, method=method.lower(), threshold=float(threshold))
			if not returned_compounds:
				return render(request, 'cmpd_reg/similarity.html', {
					'searched_smiles': smiles,
					'searched_compound_image': searched_img_tag,
					'no_matches': True,
				})
			matched_output = []
			id_list = []
			for structure in returned_compounds:
				cmpd_obj = {}
				cmpd_obj['id'] = str(structure.id)
				id_list.append(str(structure.id))
				cmpd_obj['smiles'] = structure.smiles
				cmpd_obj['inchi'] = structure.inchi
				cmpd_obj['property'] = structure.property
				pil_image = create_compound_image(structure.molecule)
				cmpd_obj['image'] = create_image_tag(pil_image)
				matched_output.append(cmpd_obj)
			return render(request, 'cmpd_reg/similarity.html', {
				'searched_smiles': smiles,
				'searched_compound_image': searched_img_tag,
				'matched_structures': matched_output,
				'id_list': mark_safe(json.dumps(id_list, cls=DjangoJSONEncoder)),
			})
		else:
			return render(request, 'cmpd_reg/similarity.html', {
				'smiles_error': "Error! Could not convert SMILES '{}' to image".format(smiles),
			})

	return render(request, 'cmpd_reg/similarity.html', {})

def bulk_loader(request):
	if request.method == "POST":
		inchi = request.POST.get("inchi", False)
		inchi_check = True if inchi == 'on' else False
		output_stats = {'new':0,'existing':0,'errors':{'count':0,'report':[]}}
		reader = csv.DictReader(request.FILES['file'])
		row_count = 1
		for row in reader:
			row_count += 1
			if 'SMILES' not in row.keys():
				return render(request, 'cmpd_reg/bulk_loader.html', {
					'load_error': "Error! CSV must contain a column header named 'SMILES'"
				})
			try:
				_, is_new = get_or_create_compound(row['SMILES'], clean=True, inchi_standard=inchi_check)
				if is_new:
					output_stats['new'] += 1
				else:
					output_stats['existing'] += 1
			except Exception, e:
				print str(e)
				output_stats['errors']['count'] += 1
				output_stats['errors']['report'].append({'smiles':row['SMILES'],'msg':str(e),'row':str(row_count)})
		return render(request, 'cmpd_reg/bulk_loader.html', {
			'load_results': output_stats
		})

	return render(request, 'cmpd_reg/bulk_loader.html', {})

def download_structures(request):
	file_ids = request.GET.getlist('ids[]', None)
	s = write_to_sdf(file_ids)
	response = HttpResponse(s.getvalue(), content_type='text/plain')
	response['Content-Disposition'] = 'attachment; filename="%s.sdf"' % "search_results_{}".format(time.strftime("%m%d%Y"))
	return response
