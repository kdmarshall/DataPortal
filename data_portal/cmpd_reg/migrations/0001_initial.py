# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django_rdkit.models.fields


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='MoleculeModel',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('molecule', django_rdkit.models.fields.MolField()),
            ],
        ),
        migrations.CreateModel(
            name='SmilesModel',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('smiles', models.CharField(max_length=2048, blank=True)),
                ('molecule', django_rdkit.models.fields.MolField(null=True)),
            ],
        ),
    ]
