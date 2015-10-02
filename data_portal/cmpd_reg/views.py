from django.shortcuts import render
from django.http import HttpResponse
from django.core.files.base import ContentFile

from utils import create_compound_image, get_or_create_compound

import StringIO

def creg_index(request):
	if request.method == 'POST':
		smiles = request.POST.get("smiles", None)
		result_compound, is_new = get_or_create_compound(smiles)
		searched_image = create_compound_image(smiles, smiles=True)
		if searched_image:
			output = StringIO.StringIO()
			searched_image.save(output, "PNG")
			contents = output.getvalue().encode("base64")
			output.close()
			img_tag = '<img src="data:image/png;base64,' + contents + '" />'
			return render(request, 'cmpd_reg/cmpd-reg.html', {
				'searched_smiles': smiles,
				'searched_compound_image': img_tag,
				'new_compound': str(is_new),
			})
		else:
			return render(request, 'cmpd_reg/cmpd-reg.html', {
				'smiles_error': 'Error! Could not convert smiles {} to image'.format(smiles),
			})

	return render(request, 'cmpd_reg/cmpd-reg.html', {})
