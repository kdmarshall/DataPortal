# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django_rdkit.operations import GiSTIndex

class Migration(migrations.Migration):

    dependencies = [
        ('cmpd_reg', '0001_initial'),
    ]

    operations = [
	GiSTIndex('MoleculeModel','molecule')
    ]
