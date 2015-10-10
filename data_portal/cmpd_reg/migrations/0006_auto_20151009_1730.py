# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django_rdkit.operations import GiSTIndex

class Migration(migrations.Migration):

    dependencies = [
        ('cmpd_reg', '0005_remove_compound_image'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='fingerprint',
            name='sfp',
        ),
        migrations.RemoveField(
            model_name='property',
            name='logp',
        ),
	
    ]
