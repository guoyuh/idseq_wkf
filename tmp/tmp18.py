
import pandas as pd

import re


a="k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Dipodascaceae|g__Yarrowia%61"
b="k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Eremothecium|s__Eremothecium_gossypii%2"
c="k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Saccharomycetaceae|g__Saccharomyces|s__Saccharomyces_paradoxus%3"
genus_re = re.search(r'.*?\|s__(\w+)%(\d+)$', b)
genus_re2 = re.search('\|s__(\w+)%(\d+)', b)

genus_re3 = re.search(r".*?\|g__(\w+)\|s__(\w+)%(\d+)",b)

print(genus_re.group())
print(genus_re.group(1))
print(genus_re.group(2))

print("=================================================")

print(genus_re2.group())
print(genus_re2.group(1))
print(genus_re2.group(2),type(genus_re2.group(2)))
print("------------------------------------------------")

print(genus_re3.group())
print(genus_re3.group(1))
print(genus_re3.group(2))
print(genus_re3.group(3))
print("++++++++++++++++++++++++++++++++++++++++++++++++++++")
fungi_re = re.match('k__Eukaryota\|k__Fungi\|.*?\|g__(\w+)\|s__(\w+)%(\d+)', c)
print(fungi_re.group())
print(fungi_re.group(1))
print(fungi_re.group(2))
print(fungi_re.group(3))

print("***************************************************")
d="k__Eukaryota|k__Fungi|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Debaryomycetaceae|g__Candida|s__Candida_parapsilosis%133"
e="k__Eukaryota|p__Apicomplexa|c__Conoidasida|o__Eucoccidiorida|f__Cryptosporidiidae|g__Cryptosporidium|s__Cryptosporidium_hominis%101"
fungi_re = re.match('k__Eukaryota\|k__Fungi\|.*?\|g__(\w+)\|s__(\w+)%(\d+)', d)
parasites_re = re.match('k__Eukaryota\|.*?\|([kpcofgs])__(\w+)%(\d+)', d)

if fungi_re:
    print("k__Fungi")

if parasites_re:
    print("parasites")

fungi_re2 = re.match('k__Eukaryota\|k__Fungi\|.*?\|g__(\w+)\|s__(\w+)%(\d+)', e)
parasites_re2 = re.match('k__Eukaryota\|.*?\|([kpcofgs])__(\w+)%(\d+)', e)

if (not parasites_re2 is None) & (fungi_re2 is None):
    print("only parasites")