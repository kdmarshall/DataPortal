from django.contrib import admin
from models import (Compound,
					Property,
					Fingerprint)

class CompoundAdmin(admin.ModelAdmin):
	pass

class PropertyAdmin(admin.ModelAdmin):
	pass

class FingerprintAdmin(admin.ModelAdmin):
	pass

admin.site.register(Compound, CompoundAdmin)
admin.site.register(Property, PropertyAdmin)
admin.site.register(Fingerprint, FingerprintAdmin)
