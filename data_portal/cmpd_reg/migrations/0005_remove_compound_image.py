# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django_rdkit.operations import GiSTIndex

class Migration(migrations.Migration):

    dependencies = [
        ('cmpd_reg', '0004_auto_20151002_2137'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='compound',
            name='image',
        ),
	GiSTIndex('Compound','molecule')
    ]
