from django.contrib import admin
from models import (Compound,
					Property,
					Fingerprint)

class CompoundAdmin(admin.ModelAdmin):
	pass

admin.site.register(Compound, CompoundAdmin)
