#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

from Bio import SeqIO
import pandas as pd
import re
infile = '/mnt/data/kraken2_mask/library/bacteria/library.fna'


# with open (infile,'r') as fh:
#     record = SeqIO.parse(fh,'fasta')

#     for seq in record:
#         #print(seq.id,'|||||||',seq.name)
#         if seq.id.split("|")[1] == "28131":
#             print(seq.id)
#             print(len(seq.seq))





a = "/mnt/home/huanggy/project/20220408/result/21JS944006/align/21JS944006.refseq.cove.xls"

# tb = pd.read_table(a,sep="\t",header=0)
#
# print(tb.groupby("Genome")["Covered Bases","Covered Fraction","Variance","Length" ,"my_refseq_Read_Count","RPKM"].mean())
# sub_tb = tb.groupby("Genome")["Covered Bases","Covered Fraction","Variance","Length" ,"my_refseq_Read_Count","RPKM"].mean()
# sub_tb.to_csv("/mnt/home/huanggy/project/20220408/result/21JS944006/align/21JS944006.tmp",sep="\t")

a = {'Schizosaccharomyces pombe': 1, 'Pichia kudriavzevii': 1, 'Torulaspora delbrueckii': 1, 'Zygosaccharomyces rouxii': 1, 'Kockovaella imperatae': 1, 'Brettanomyces bruxellensis': 1, 'Aspergillus tubingensis': 1, 'Penicillium griseofulvum': 1, 'Penicillium roqueforti': 1, 'Pleurotus ostreatus': 1, 'Ascochyta rabiei': 1, 'Colletotrichum truncatum': 1, '[Candida] glabrata': 1, 'Candida parapsilosis': 1, 'Diutina rugosa': 1, 'Alternaria alternata': 1, 'Brettanomyces nanus': 1, 'Saccharomyces paradoxus': 1, 'Penicillium expansum': 1, 'Kluyveromyces lactis': 1, 'Purpureocillium lilacinum': 1, 'Nosema ceranae': 1, 'Aspergillus thermomutatus': 1, 'Aspergillus candidus': 1, 'Zygotorulaspora mrakii': 1, 'Fusarium subglutinans': 1, 'Lasiodiplodia theobromae': 1, 'Eremothecium sinecaudum': 1, '[Candida] haemuloni': 1, 'Scheffersomyces spartinae': 1, 'Wickerhamiella sorbophila': 1, 'Torulaspora globosa': 1, 'Suillus bovinus': 1, 'Suillus paluster': 1, 'Suillus subalutaceus': 1, 'Suillus subaureus': 1, 'Fusarium venenatum': 1, 'Debaryomyces fabryi': 1, 'Trichoderma citrinoviride': 1, 'Tilletiopsis washingtonensis': 1, 'Penicillium solitum': 1, 'Kazachstania barnettii': 1, 'Linderina pennispora': 1, 'Aspergillus caelatus': 1, 'Lobosporangium transversale': 1, 'Arthroderma uncinatum': 1, 'Aspergillus viridinutans': 1, 'Malassezia restricta': 1, 'Malassezia pachydermatis': 1, 'Diaporthe citri': 1, 'Moesziomyces antarcticus': 1, 'Botrytis porri': 1, 'Aspergillus udagawae': 1, 'Exophiala spinifera': 1, 'Drechmeria coniospora': 1, 'Apiotrichum porosum': 1, 'Aspergillus bombycis': 1, 'Letharia columbiana': 1, 'Ramularia collo-cygni': 1, 'Suillusplorans': 1, 'Alternaria atra': 1, 'Ustilago hordei': 1, 'Talaromyces rugulosus': 1, 'Cercospora beticola': 1, 'Aspergillus pseudotamarii': 1, 'Botrytis byssoidea': 1, 'Sparassis crispa': 1, 'Pyricularia grisea': 1, 'Mollisia scopiformis': 1, 'Alternaria arborescens': 1, 'Marasmius oreades': 1, 'Aspergillus chevalieri': 1, 'Fusarium mangiferae': 1, 'Aspergillus alliaceus': 1, 'Exophiala mesophila': 1, 'Cryptococcus neoformans var. neoformans JEC21': 1, 'Exophiala oligosperma': 1, 'Acaromyces ingoldii': 1, 'Coccidioides posadasii C735 delta SOWgp': 1, 'Aspergillus nidulans FGSC A4': 1, 'Fusarium graminearum PH-1': 1, 'Fusarium coffeatum': 1, 'Cryptococcus neoformans var. grubii H99': 1, 'Diplodia corticola': 1, 'Candida albicans SC5314': 1, 'Ustilago maydis 521': 1, 'Coprinopsis cinerea okayama7#130': 1, 'Pyricularia oryzae 70-15': 1, 'Neohortaea acidophila': 1, 'Coccidioides immitis RS': 1, 'Verruconis gallopava': 1, 'Fonsecaea monophora': 1, 'Paecilomyces variotii': 1, 'Sporisorium graminicola': 1, 'Cryptococcus neoformans var. neoformans B-3501A': 1, 'Yarrowia lipolytica CLIB122': 1, 'Debaryomyces hansenii CBS767': 1, 'Eremothecium gossypii ATCC 10895': 1, 'Encephalitozoon cuniculi GB-M1': 1, 'Aspergillus lentulus': 1, 'Meyerozyma guilliermondii ATCC 6260': 1, 'Candida tropicalis MYA-3404': 1, 'Macroventuria anomochaeta': 1, 'Chaetomium globosum CBS 148.51': 1, 'Clavispora lusitaniae ATCC 42720': 1, 'Westerdykella ornata': 1, 'Parastagonospora nodorum SN15': 1, 'Scheffersomyces stipitis CBS 6054': 1, 'Daldinia childiae': 1, 'Aspergillus fumigatus Af293': 1, 'Aspergillus fischeri NRRL 181': 1, 'Botrytis cinerea B05.10': 1,'Aspergillus flavus NRRL3357': 1, 'Fusarium verticillioides 7600': 1, 'Zymoseptoria tritici IPO323': 1, 'Uncinocarpus reesii 1704': 1, 'Histoplasma capsulatum NAm1': 1, 'Aspergillus terreus NIH2624': 1, 'Pseudogymnoascus verrucosus': 1, 'Aspergillus clavatus NRRL 1': 1, 'Exophiala xenobiotica': 1, 'Neurospora crassa OR74A': 1, 'Cryptococcus gattii WM276': 1, 'Lodderomyces elongisporus NRRL YB-4239': 1, 'Pseudocercospora fijiensis CIRAD86': 1, 'Trematosphaeria pertusa': 1, 'Trichoderma gamsii': 1, 'Schizosaccharomyces japonicus yFS275': 1, 'Trichoderma virens Gv29-8': 1, 'Puccinia graminis f. sp. tritici CRL 75-36-700-3': 1, '[Candida] pseudohaemulonii': 1, 'Aspergillus niger CBS 513.88': 1, 'Malassezia globosa CBS 7966': 1, 'Pyrenophora tritici-repentis Pt-1C-BFP': 1, 'Fusarium oxysporum f. sp. lycopersici 4287': 1, 'Trichoderma reesei QM6a': 1, 'Vanderwaltozyma polyspora DSM 70294': 1, 'Talaromyces stipitatus ATCC 10500': 1, 'Talaromyces marneffei ATCC 18224': 1, 'Trichoderma atroviride IMI 206040': 1, 'Rhizoctonia solani': 1, 'Ogataea polymorpha': 1, 'Pseudovirgaria hyperparasitica': 1, 'Schizosaccharomyces octosporus yFS286': 1, 'Laccaria bicolor S238N-H82':1, '[Candida] auris': 1, 'Verticillium dahliae VdLs.17': 1, 'Penicillium rubens Wisconsin 54-1255': 1, 'Paracoccidioides lutzii Pb01': 1, 'Paracoccidioides brasiliensis Pb18': 1, 'Aspergillus oryzae RIB40': 1, 'Neurospora tetrasperma FGSC 2508': 1, 'Podospora anserina S mat+': 1, 'Verticillium alfalfae VaMs.102': 1, 'Nannizzia gypsea CBS 118893': 1, 'Microsporum canis CBS 113480': 1, 'Saccharomyces cerevisiae S288C': 1, 'Lachancea thermotolerans CBS 6340': 1, 'Blastomyces gilchristii SLH14081': 1, 'Trichophyton rubrum CBS 118892': 1, 'Letharia lupina': 1, 'Scedosporium apiospermum': 1, 'Cladophialophora immunda': 1, 'Thermothelomyces thermophilus ATCC 42464': 1, 'Candida dubliniensis CD36': 1, 'Mytilinidion resinicola': 1, 'Thermothielavioides terrestris NRRL 8126': 1, 'Tremella mesenterica DSM 1558': 1, 'Serpula lacrymans var. lacrymans S7.9': 1, 'Schizophyllum commune H4-8': 1, 'Rhodotorula graminis WP1': 1, 'Yamadazyma tenuis ATCC 10573': 1, 'Agaricus bisporus var. burnettii JB137-S8': 1, 'Fibroporia radiculosa': 1, 'Spathaspora passalidarum NRRL Y-27907': 1, 'Komagataella phaffii GS115': 1, 'Gaeumannomyces tritici R3-111a-1': 1, 'Colletotrichum graminicola M1.001': 1, 'Spizellomyces punctatus DAOM BR117': 1, 'Phanerochaete carnosa HHB-10118-sp': 1, 'Schizosaccharomyces cryophilus OY26': 1, 'Beauveria bassiana ARSEF 2860': 1, 'Metarhizium acridum CQMa 102': 1, 'Metarhizium robertsii ARSEF 23': 1, 'Grosmannia clavigera kw1407': 1, 'Pseudogymnoascus destructans': 1, 'Tuber melanosporum Mel28': 1, 'Fusarium oxysporum NRRL 32931': 1, 'Fusarium vanettenii 77-13-4': 1, 'Cryphonectria parasitica EP155': 1, 'Trichophyton verrucosum HKI 0517': 1, 'Trichophyton benhamiae CBS 112371': 1, 'Bipolaris maydis ATCC 48331': 1, 'Sclerotinia sclerotiorum 1980 UF-70': 1, 'Bipolaris sorokiniana ND90Pr': 1, 'Gloeophyllum trabeum ATCC 11539': 1, 'Postia placenta MAD-698-R-SB12': 1, 'Wallemia mellicola CBS 633.66': 1, 'Exserohilum turcica Et28A': 1, 'Lindgomyces ingoldianus': 1, 'Wickerhamomyces anomalus NRRL Y-366-8': 1, 'Batrachochytrium dendrobatidis JAM81': 1, 'Colletotrichum fructicola': 1, 'Colletotrichum siamense': 1, 'Aspergillus aculeatus ATCC 16872': 1, 'Sphaerulina musiva SO2202': 1, 'Fomitiporia mediterranea MF3/22': 1, 'Saitoella complicata NRRL Y-17804': 1, 'Baudoinia panamericana UAMH 10762': 1, 'Trametes versicolor FP-101664 SS1': 1, 'Stereum hirsutum FP-91666 SS1': 1, 'Dichomitus squalens LYAD-421 SS1': 1, 'Punctularia strigosozonata HHB-11173 SS5': 1, 'Coniophora puteana RWD-64-598 SS2': 1, 'Rhizophagus irregularis DAOM 181602=DAOM 197198': 1, 'Heterobasidion irregulare TC 32-1': 1, 'Melampsora larici-populina 98AG31': 1, 'Orbilia oligospora ATCC 24927': 1, 'Chaetomium thermophilum var. thermophilum DSM 1495': 1, 'Colletotrichum higginsianum IMI 349063': 1, 'Pichia membranifaciens NRRL Y-2026': 1, 'Phycomyces blakesleeanus NRRL 1555(-)': 1, 'Mixia osmundae IAM 14324': 1, 'Sordaria macrospora k-hell': 1, 'Sugiyamaella lignohabitans': 1, 'Fonsecaea nubica': 1, 'Amorphotheca resinae ATCC 22711': 1, 'Exophiala dermatitidis NIH/UT8656': 1, 'Metschnikowia bicuspidata var. bicuspidata NRRL YB-4993': 1, 'Ogataea angusta': 1, 'Ogataea parapolymorpha DL-1': 1, 'Encephalitozoon intestinalis ATCC 50506': 1, 'Cutaneotrichosporon oleaginosum': 1, 'Nematocida parisii ERTm1': 1, 'Encephalitozoon hellem ATCC 50504': 1, 'Bipolaris zeicola 26-R-13': 1, 'Bipolaris oryzae ATCC 44560': 1, 'Bipolaris victoriae FI3': 1, 'Eremothecium cymbalariae DBVPG#7215': 1, 'Agaricus bisporus var. bisporus H97': 1, 'Vavraia culicis subsp. floridensis': 1, 'Cordyceps militaris CM01': 1, 'Trichoderma harzianum CBS 226.95': 1, 'Cyberlindnera jadinii NRRL Y-1542': 1, 'Hyphopichia burtonii NRRL Y-1933': 1, 'Babjeviella inositovora NRRL Y-12698': 1, 'Suhomyces tanzawaensis NRRL Y-17324': 1, 'Leptosphaeria maculans JN3': 1, 'Vittaforma corneae ATCC 50505': 1, 'Kluyveromyces marxianus DMKU3-1042': 1, 'Fusarium pseudograminearum CS3096': 1, 'Aspergillus versicolor CBS 583.65': 1, 'Aspergillus sydowii CBS 593.65': 1, 'Tilletiaria anomala UBC 951': 1, 'Wickerhamomyces ciferrii': 1, 'Trichoderma asperellum CBS 433.97': 1, 'Aureobasidium pullulans EXF-150': 1, 'Aureobasidium melanogenum CBS 110374': 1, 'Aureobasidium namibiae CBS 147.97': 1, 'Aureobasidium subglaciale EXF-2481': 1, 'Verticillium nonalfalfae': 1, 'Naumovozyma castellii CBS 4309': 1, 'Aspergillus luchuensis': 1, 'Pneumocystis murina B123': 1, 'Naumovozyma dairenensis CBS 421': 1, 'Tetrapisispora blattae CBS 6284': 1, 'Tetrapisispora phaffii CBS 4417': 1, 'Kazachstania africana CBS 2517': 1, 'Kazachstania naganishii CBS 8797': 1, "Marssonina brunnea f. sp. 'multigermtubi' MB_m1": 1, 'Aspergillus wentii DTO 134E9': 1, 'Penicilliopsis zonata CBS 506.65': 1, 'Zasmidium cellare ATCC 36951': 1, 'Saccharomyces eubayanus': 1, 'Metarhizium album ARSEF 1941': 1, 'Cordyceps fumosorosea ARSEF 2679': 1, 'Fusarium odoratissimum NRRL 54006': 1, 'Phialemoniopsis curvata': 1, 'Geosmithia morbida': 1, 'Colletotrichum karsti': 1, 'Hyaloscypha bicolor E': 1, 'Glarea lozoyensis ATCC 20868': 1, 'Rhodotorula toruloides NP11': 1, 'Candida orthopsilosis Co 90-125': 1, 'Pseudomassariella vexata': 1, 'Didymella exigua CBS 183.55': 1, 'Ustilaginoidea virens': 1, 'Aspergillus glaucus CBS 516.65': 1, 'Coniosporium apollinis CBS 100218': 1, 'Cucurbitaria berberidis CBS 394.84': 1, 'Penicillium digitatum Pd1': 1, 'Aplosporella prunicola CBS 121167': 1, 'Encephalitozoon romaleae SJ-2008': 1, 'Capronia coronata CBS 617.96': 1, 'Capronia epimyces CBS 606.96': 1, 'Cladophialophora psammophila CBS 110553': 1, 'Cladophialophora yegresii CBS 114405': 1, 'Exophiala aquamarina CBS 119918': 1, 'Trichosporon asahii var. asahii CBS 2479': 1, 'Alternaria burnsii': 1, 'Talaromyces amestolkiae': 1, 'Colletotrichum orchidophilum': 1, 'Colletotrichum scovillei': 1, 'Colletotrichum aenigma': 1, 'Aspergillus tanneri': 1, 'Aspergillus puulaauensis': 1, 'Cyphellophora europaea CBS 101466': 1, 'Fusarium proliferatum ET1': 1, 'Pestalotiopsis fici W106-1': 1, 'Malassezia sympodialis ATCC 42132': 1, '[Candida] duobushaemulonis': 1, 'Lachancea lanzarotensis': 1, 'Endocarpon pusillum Z07020': 1, 'Metarhizium brunneum ARSEF 3297': 1, 'Anthracocystis flocculosa PF-1': 1, 'Cladophialophora carrionii CBS 160.54': 1, 'Fusarium fujikuroi IMI 58289': 1, 'Meira miltonrushii': 1, 'Phaeoacremonium minimum UCRPA7': 1, 'Cryptococcus wingfieldii CBS 7118': 1, 'Cryptococcus amylolentus CBS 6039': 1, 'Kwoniella pini CBS 10737': 1, 'Kwoniella bestiolae CBS 10118': 1, 'Kwoniella dejecticola CBS 10117': 1, 'Kwoniella mangroviensis CBS 8507': 1, 'Wallemia ichthyophaga EXF-994': 1, 'Pseudozyma hubeiensis SY62': 1, 'Sodiomyces alkalinus F11': 1, 'Laetiporus sulphureus 93-53': 1, 'Dissoconium aciculare CBS 342.82': 1, 'Lachnellula hyalina': 1, 'Xylona heveae TC161': 1, 'Rhizopus microsporus ATCC 52813': 1, 'Aspergillus welwitschiae': 1, 'Ascoidea rubescens DSM 1968': 1, 'Ordospora colligata OC4': 1, 'Kalmanozyma brasiliensis GHG001': 1, 'Fonsecaea erecta': 1, 'Pochonia chlamydosporia 170': 1, 'Kuraishia capsulata CBS 1993': 1, 'Aspergillus ruber CBS 135680': 1, 'Eremomyces bilateralis CBS 781.70': 1, 'Dothidotthia symphoricarpi CBS 119687': 1, 'Aspergillus campestris IBT 28561': 1, 'Aspergillus steynii IBT 23096': 1, 'Aspergillus novofumigatusIBT 16806': 1, 'Aspergillus ochraceoroseus IBT 24754': 1, 'Sporothrix schenckii 1099-18': 1, 'Sporothrix brasiliensis 5110': 1, 'Rasamsonia emersonii CBS 393.64': 1, 'Pneumocystis jirovecii RU7': 1, 'Pneumocystis carinii B80': 1, 'Talaromyces atroroseus': 1, 'Fonsecaea pedrosoi CBS 271.37': 1, 'Rhinocladiella mackenziei CBS 650.93': 1, 'Cladophialophora bantiana CBS 173.52': 1, 'Fonsecaea multimorphosa CBS 102226': 1, 'Aspergillus neoniger CBS 115656': 1, 'Aspergillus vadensis CBS 113365': 1, 'Aspergillus japonicus CBS 114.51': 1, 'Aspergillus piperis CBS 112811': 1, 'Aspergillus eucalypticola CBS 122712': 1, 'Aspergillus uvarum CBS 121591': 1, 'Aspergillus ibericus CBS 121593': 1, 'Aspergillus costaricaensis CBS 115574': 1, 'Aspergillus fijiensis CBS 313.89': 1, 'Aspergillus heteromorphus CBS 117.55': 1, 'Aspergillus aculeatinus CBS 121060': 1, 'Aaosphaeria arxii CBS 175.79': 1, 'Aspergillus niger CBS 101883': 1, 'Aspergillus brunneoviolaceus CBS 621.78': 1, 'Aspergillus sclerotioniger CBS 115572': 1, 'Guyanagaster necrorhizus MCA 3950': 1, 'Aspergillus homomorphus CBS 101889': 1, 'Aspergillus saccharolyticus JOP 1030-1': 1, 'Paraphaeosphaeria sporulosa': 1, 'Botrytis sinoallii': 1, 'Mitosporidium daphniae': 1, 'Aspergillus pseudonomiae': 1, 'Aspergillus nomiae NRRL 13137': 1, 'Aspergillus pseudoviridinutans': 1, 'Ceraceosorus guamensis': 1, 'Fusarium tjaetaba': 1, 'Jaminaea rosea': 1, 'Pyricularia pennisetigena': 1, 'Phialophora attinorum': 1, 'Pseudomicrostroma glucosiphilum': 1, 'Kwoniella shandongensis': 1, 'Cantharellus anzutake': 1, 'Synchytrium microbalum': 1, 'Aspergillus mulundensis': 1, 'Penicillium arizonense': 1, 'Dacryopinax primogenitus': 1, 'Suillus clintonianus': 1, 'Suillus discolor': 1, 'Suillus fuscotomentosus': 1, 'Ogataea haglerorum': 1, 'Botrytis fragariae': 1, 'Mycena indigotica': 1, 'Botrytis deweyae': 1, 'Saprochaete ingens': 1, 'Venustampulla echinocandica': 1, 'Protomyces lactucae-debilis': 1}


# for k in a.keys():
#     if k.find("Aspergillus niger") != -1 :
#         print("======================")


def is_fungi(fungi_assemble_refseq="/mnt/data/kraken2_mask/library/fungi/assembly_summary.txt"):
    fungi_organism_name = {}
    with open(fungi_assemble_refseq) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue

            organism_name = line.strip().split("\t")[7]
            if organism_name not in fungi_organism_name:
                fungi_organism_name[organism_name] = 1
            else:
                pass
    return fungi_organism_name

def get_category_counts(before_filter_mpa_report):
    """
    :param before_filter_mpa_report: kreport2mpa ---> before_filter_mpa_report
    :return:
    """
    fungi_amount = 0
    bacteria_amount = 0
    parasites_amount = 0
    viruses_amount = 0
    fungi_organism_name = is_fungi(fungi_assemble_refseq="/mnt/data/kraken2_mask/library/fungi/assembly_summary.txt")

    with open(before_filter_mpa_report, "r") as fh:
        for line in fh:
            line = line.strip()
            # fungi_re = re.match('k__Eukaryota\|k__Fungi\t(\d+)', line)
            # Eukaryota_re = re.match('k__Eukaryota\|.*?\|g__(\w+)\|s__(\w+)\t(\d+)', line)

            # 真核k__Eukaryota
            Eukaryota_re = re.match('k__Eukaryota\t(\d+)', line)
            if Eukaryota_re:
                # 真核生物
                Eukaryota_amount = int(Eukaryota_re.group(1))
                print("==============================>", Eukaryota_amount)

            species_re = re.search(r'\|s__(\w+)\t(\d+)', line)
            if species_re:
                # 真菌
                species_name = species_re.group(1).replace("_", " ", 2)
                if species_name == "Aspergillus niger":
                    print("Aspergillus niger---------------------------->")
                # if species_name in fungi_organism_name:
                #     fungi_amount +=1
                #     print (fungi_amount)
                for k in fungi_organism_name.keys():
                    if k.find(species_name) != -1:
                        fungi_amount += 1
                        print(fungi_amount)

            # 细菌 k__Bacteria
            bacteria_re = re.match('k__Bacteria\t(\d+)', line)
            if bacteria_re:
                bacteria_amount = int(bacteria_re.group(1))
                print('bacteria_amount', bacteria_amount)

            viruses_re = re.match('k__Viruses\t(\d+)', line)
            if viruses_re:
                viruses_amount = int(viruses_re.group(1))
                print('viruses_amount', viruses_amount)

    parasites_amount = Eukaryota_amount - fungi_amount
    print("Eukaryota_amount:++++++++++++++",Eukaryota_amount)
    print("fungi_amount:", fungi_amount, "bacteria_amount:", bacteria_amount, "parasites_amount:", parasites_amount,
          "viruses_amount:", viruses_amount)
    return fungi_amount, bacteria_amount, parasites_amount, viruses_amount

fungi_amount, bacteria_amount, parasites_amount, viruses_amount = get_category_counts(before_filter_mpa_report="/mnt/home/huanggy/project/20220408/result/21JS944006/kraken2/21JS944006.kreport.mpa.txt")