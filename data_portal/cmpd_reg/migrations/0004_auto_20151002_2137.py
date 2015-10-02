# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cmpd_reg', '0003_auto_20151001_0159'),
    ]

    operations = [
        migrations.AlterField(
            model_name='compound',
            name='image',
            field=models.ImageField(upload_to='cmpd_images/', null=True, verbose_name='Compound Image'),
        ),
    ]
