# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import django.db.models.deletion
import django_rdkit.models.fields


class Migration(migrations.Migration):

    dependencies = [
        ('cmpd_reg', '0002_create_compound_molecule_index'),
    ]

    operations = [
        migrations.CreateModel(
            name='Compound',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('smiles', models.CharField(max_length=2048)),
                ('molecule', django_rdkit.models.fields.MolField()),
                ('inchi', models.CharField(max_length=2048, null=True, blank=True)),
                ('inchi_key', models.CharField(max_length=2048, null=True, blank=True)),
                ('ctab', models.TextField(null=True, blank=True)),
                ('datetime_loaded', models.DateTimeField(auto_now_add=True, verbose_name='Datetime When Compound Registered')),
                ('image', models.FileField(upload_to='cmpd_images/', verbose_name='Compound Image')),
            ],
            options={
                'db_table': 'compound',
            },
        ),
        migrations.CreateModel(
            name='Fingerprint',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('bfp', django_rdkit.models.fields.BfpField(null=True, verbose_name='Binary Fingerprint')),
                ('sfp', django_rdkit.models.fields.SfpField(null=True, verbose_name='Sparse Vector Fingerprint')),
            ],
            options={
                'db_table': 'compound_fingerprint',
                'verbose_name': 'Compound Fingerprint',
                'verbose_name_plural': 'Compound Fingerprints',
            },
        ),
        migrations.CreateModel(
            name='Property',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('amw', models.DecimalField(null=True, verbose_name='Average Molecular Weight', max_digits=10, decimal_places=3, blank=True)),
                ('hbd', models.PositiveSmallIntegerField(null=True, verbose_name='Hydrogen Bond Donors', blank=True)),
                ('hba', models.PositiveSmallIntegerField(null=True, verbose_name='Hydrogen Bond Acceptors', blank=True)),
                ('logp', models.DecimalField(null=True, verbose_name='LogP', max_digits=10, decimal_places=3, blank=True)),
                ('tpsa', models.DecimalField(null=True, verbose_name='Total Polar Surface Area', max_digits=10, decimal_places=3, blank=True)),
            ],
            options={
                'db_table': 'compound_property',
                'verbose_name': 'Compound Property',
                'verbose_name_plural': 'Compound Properties',
            },
        ),
        migrations.DeleteModel(
            name='MoleculeModel',
        ),
        migrations.DeleteModel(
            name='SmilesModel',
        ),
        migrations.AddField(
            model_name='compound',
            name='fingerprint',
            field=models.OneToOneField(null=True, on_delete=django.db.models.deletion.SET_NULL, blank=True, to='cmpd_reg.Fingerprint', verbose_name='Compound Fingerprint'),
        ),
        migrations.AddField(
            model_name='compound',
            name='property',
            field=models.OneToOneField(null=True, on_delete=django.db.models.deletion.SET_NULL, blank=True, to='cmpd_reg.Property', verbose_name='Compound Property'),
        ),
    ]
