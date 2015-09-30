from __future__ import unicode_literals

from django_rdkit import models as rdkit_models

class MoleculeModel(rdkit_models.Model):

    molecule = models.MolField()


class SmilesModel(rdkit_models.Model):

    smiles = models.CharField(max_length=2048, blank=True, null=False)
    molecule = models.MolField(null=True)

