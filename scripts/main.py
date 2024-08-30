from cyvcf2 import VCF

i = 0
for variant in VCF("./data/raw/600K_immigrants.vcf"):
    print(variant.gt_bases)
    print(variant.ID)
    i += 1
    if i == 10:
        break

# extract the individuals and variants to a dataframe to then train a random
# forest:
