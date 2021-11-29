# -*- coding: utf-8 -*
# @Time    : 2021/03/06 16:31
# @Author  : yellow_huang
import json,time, requests,sys,os,re
"""
Ref : https://api.oncokb.org/oncokb-website/api#annotate-mutations-by-protein-change
"""
def getbyProteinChange(genename,protein_change,API_Token='90ec77a3-cfc6-44b7-9cd8-abd952a3749b'):
    # cmd = (
    #     f'curl -X GET "https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol={genename}&alteration={protein_change}" '
    #     f'-H "accept: application/json" '
    #     f'-H "Authorization: Bearer {API_Token}" '
    # )
    baseurl = f'https://oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol={genename}&alteration={protein_change}'
    headers = {'Authorization': 'Bearer 90ec77a3-cfc6-44b7-9cd8-abd952a3749b', 'accept': 'application/json'}
    r = requests.get(url=baseurl,timeout=30,headers = headers)
    if r.status_code != 200:
        print("connection err!!! return {}".format(r.status_code), file=sys.stderr)
        print(r.url)
        exit(1)
    return  r.content
    #out = subprocess.check_output(cmd,shell=True)
    #return  out

def myreadline(data):
    with open (data,'r') as fh:
        for line in fh:
            yield line.strip()


def anno_cosmic(gene,protein_change,cosmic_database):

    if re.search(r'\*',protein_change):
        protein_change = protein_change.replace('*',"\*")

    ##使用精确搜索，不要使用模糊匹配 grep -w
    ## grep 'TP53' CosmicCodingMuts.vcf | grep -E -w 'E287\*'
    cmd = (
        f'grep -E  {gene} {cosmic_database} | grep -E -w {protein_change} '
    )
    res = os.popen(cmd)
    res_list = res.readlines()
    anno_cosmic_list = []
    if res_list:
        for item in res_list:
            print(item.split('\t'))
            chrom, pos, cosv, ref, alt, *_, info = item.split('\t')
            matchgene = re.search(r'GENE=(.*?);', info).groups()[0]  ###是选择基因还是选择转录本？？？？？？？？？？？？？
            match_cosm = re.search(r'LEGACY_ID=(COSM\d+);', info).groups()[0]
            match_cds = re.search(r'CDS=(c\..*?);', info).groups()[0]
            match_AA = re.search(r'AA=(p\..*?);', info).groups()[0] #HGVS_SHOT
            HGVSC = re.search(r'GVSC=.*?:(c\..*?);', info).groups()[0]
            HGVSP = re.search(r'HGVSP=.*?:(p\..*?);', info).groups()[0]
            HGVSG = re.search(r'HGVSG=(\d+|X|Y):(g\..*?);', info).groups()
            HGVSG = ':'.join([HGVSG[0], HGVSG[1]])

            anno_cosmic_list.append('$'.join([chrom,pos,cosv,ref,alt,match_cosm,match_cds,match_AA,HGVSC,HGVSP,HGVSG,matchgene]))

    else:
        pass

    return  anno_cosmic_list


def change_hgvsc(hgvsc):
    #c.598C>T  ====> c.C598T
    com = re.match('c.(\d+)([ATCG])>([ATCG])', hgvsc)
    if com:
        hgvsc_change = 'c.'+com.group(2)+com.group(1)+com.group(3)
        return hgvsc_change
    else:
        return hgvsc

def main(input,output):
    # data = '/home/vip/variant_call_v1.2/var_call/spyder/del_fusion_redup_hgvsp.txt'
    # fh = open('/home/vip/variant_call_v1.2/var_call/spyder/20210309.oncokb.treat.diag.oncogene.parser_out.test', 'w')
    data = input
    fh = open(output, 'w')

    header = ['Gene', 'protein_change', 'alterations', '癌症部位', 'Cancer Type', 'sub_cancerType','外显子','碱基变化','氨基酸变化','Drugs(s)','药物中文名','Level','Isoncogenic','mutationEffect', 'FDAdesc', 'variantSummary','pmids','COSV','COSM','HGVSC','HGVSP','Chr','Pos','REF','ALT','匹配到的基因','hgvsg','updatetime']
    fh.write('\t'.join(header) + '\n')

    for line in myreadline(data):
        genename, protein_change,*_ = line.split('\t')
        print(genename, protein_change)
        anno_cosmic_list = anno_cosmic(gene=genename,protein_change=protein_change,cosmic_database='/home/vip/variant_call/database/cosmic/CosmicCodingMuts.vcf')
        print(anno_cosmic_list)
        out = getbyProteinChange(genename=genename, protein_change=protein_change)
        #out = getbyProteinChange(genename = 'BRAF',protein_change = 'V600E')
        out = json.loads(out)
        EntrezGeneID = out['query'].get('entrezGeneId','-')
        HugoSymbol = out['query'].get('hugoSymbol','-')
        Alteration= out['query'].get('alteration','-')
        oncogenic = out.get('oncogenic','-')
        mutationEffect = out['mutationEffect'].get('knownEffect','-')
        SensitiveLevel=out['highestSensitiveLevel']
        ResistanceLevel=out['highestResistanceLevel']
        ishotspot=out['hotspot']
        geneSummary=out['geneSummary']
        updatetime = out.get('lastUpdate','-')
        variantSummary = out.get('variantSummary','-')
        treatL = out.get('treatments',None)
        diagnosticImplications_list = out.get('diagnosticImplications',None)
        prognosticImplications_list = out.get('prognosticImplications',None)
        if diagnosticImplications_list :
            for diagnostic in diagnosticImplications_list:
                alterations = ''.join(diagnostic.get('alterations','-'))
                level = diagnostic.get('levelOfEvidence','-')
                #主要癌症类型
                cancerType = diagnostic['tumorType']['mainType'].get('name','-')
                sub_cancerType = diagnostic['tumorType'].get('name','-')
                location = diagnostic['tumorType'].get('tissue','-')
                pmids = ','.join(diagnostic.get('pmids','-'))
                ##对某些癌症只有诊断意义 drugname = 'NA'  Variabt_sumary= 'variantSummary',FDAdescription = 'NA'
                drug_name = '-'
                FDAdescription = '-'
                if anno_cosmic_list:
                    for i in range(len(anno_cosmic_list)):
                        chrom, pos, cosv, ref, alt, match_cosm, match_cds, match_AA, HGVSC, HGVSP, HGVSG,matchgene= anno_cosmic_list[i].split('$')
                        match_cds = change_hgvsc(hgvsc=match_cds)
                        print('\t'.join(
                            [HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType,
                             "exonnum?", match_cds, match_AA, drug_name, '药物中文名', level, oncogenic, mutationEffect,
                             FDAdescription, variantSummary, pmids, cosv, match_cosm, HGVSC, HGVSP, chrom, pos, ref,
                             alt, matchgene, HGVSG, updatetime]), file=fh)
                else:
                    chrom, pos, cosv, ref, alt, match_cosm, match_cds, match_AA, HGVSC, HGVSP, HGVSG ,matchgene= ('-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-','-')
                    print('\t'.join(
                        [HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType, "exonnum?",
                         match_cds, match_AA, drug_name, '药物中文名', level, oncogenic, mutationEffect, FDAdescription,
                         variantSummary, pmids, cosv, match_cosm, HGVSC, HGVSP, chrom, pos, ref, alt, matchgene,
                         HGVSG,updatetime]), file=fh)
        if treatL:
            print('treatL 有:', treatL)
            for treat in treatL:
                #alterations = ','.join(treat["alterations"])  ###分开写
                alt_L = treat["alterations"]
                for alterations in alt_L:
                    cancerType = treat["levelAssociatedCancerType"]['mainType'].get('name','-')
                    sub_cancerType = treat["levelAssociatedCancerType"].get('name','-')
                    location = treat["levelAssociatedCancerType"].get('tissue','-')
                    drug_name =','.join([i.get('drugName','no drug') for i in treat["drugs"]])
                    药物别名 = ','.join([','.join(i.get('synonyms','no')) for i in treat["drugs"]])
                    #[ drug.drugName for drug in treat["drugs"] ]
                    level = treat.get('level','-')

                    ##FDA 批准的适应症
                    FDAdescription = ",".join(treat.get('approvedIndications','-'))
                    pmids = ','.join(treat.get("pmids","-"))
                    if anno_cosmic_list:
                        for i in range(len(anno_cosmic_list)):
                            chrom, pos, cosv, ref, alt, match_cosm, match_cds, match_AA, HGVSC, HGVSP, HGVSG,matchgene = anno_cosmic_list[i].split('$')
                            match_cds = change_hgvsc(hgvsc=match_cds)
                            # print('\t'.join(
                            #     [HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType, oncogenic,
                            #      mutationEffect, drug_name, level, FDAdescription, variantSummary, pmids, updatetime, chrom,
                            #      pos, cosv, ref, alt, match_cosm, match_cds, match_AA, HGVSC, HGVSP, HGVSG]), file=fh)
                            print('\t'.join(
                                [HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType,
                                 "exonnum?", match_cds, match_AA, drug_name, '药物中文名', level, oncogenic, mutationEffect,
                                 FDAdescription, variantSummary, pmids, cosv, match_cosm, HGVSC, HGVSP, chrom, pos, ref,
                                 alt, matchgene, HGVSG, updatetime]), file=fh)

                    else:
                        chrom, pos, cosv, ref, alt, match_cosm, match_cds, match_AA, HGVSC, HGVSP, HGVSG,matchgene = ('-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-','-')
                        print('\t'.join(
                            [HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType, "exonnum?",
                             match_cds, match_AA, drug_name, '药物中文名', level, oncogenic, mutationEffect, FDAdescription,
                             variantSummary, pmids, cosv, match_cosm, HGVSC, HGVSP, chrom, pos, ref, alt, matchgene,
                             HGVSG,updatetime]), file=fh)
        if prognosticImplications_list:
            print('有预后')
            for prognostic in prognosticImplications_list:
                alterations = ''.join(prognostic.get('alterations','-'))
                level = prognostic.get('levelOfEvidence','-')
                cancerType = prognostic['tumorType']['mainType'].get('name','-')
                sub_cancerType = prognostic['tumorType'].get('name','-')
                location = prognostic['tumorType'].get('tissue','-')
                pmids = ','.join(prognostic.get('pmids','-'))
                drug_name = '-'
                description =prognostic.get('description','-')
                if anno_cosmic_list:
                    for i in range(len(anno_cosmic_list)):
                        chrom, pos, cosv, ref, alt, match_cosm, match_cds, match_AA, HGVSC, HGVSP, HGVSG,matchgene = anno_cosmic_list[i].split('$')
                        match_cds = change_hgvsc(hgvsc=match_cds)
                        #print('\t'.join([HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType, oncogenic, mutationEffect, drug_name,level,description, variantSummary, pmids, updatetime,chrom, pos, cosv, ref, alt, match_cosm, match_cds, match_AA, HGVSC, HGVSP, HGVSG]), file=fh)

                        print('\t'.join([HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType, "exonnum?",match_cds,match_AA, drug_name,'药物中文名',level,oncogenic,mutationEffect,description, variantSummary, pmids, cosv,match_cosm,HGVSC, HGVSP, chrom, pos,ref, alt,matchgene,HGVSG,updatetime]), file=fh)
                else:
                    chrom, pos, cosv, ref, alt, match_cosm, match_cds, match_AA, HGVSC, HGVSP, HGVSG,matchgene = ('-','-','-','-','-','-','-','-','-','-','-','-')
                    print('\t'.join([HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType, "exonnum?",match_cds, match_AA, drug_name, '药物中文名', level, oncogenic, mutationEffect, description,
                         variantSummary, pmids, cosv, match_cosm, HGVSC, HGVSP, chrom, pos, ref, alt,matchgene, HGVSG,
                         updatetime]), file=fh)

        if (not treatL) and (not diagnosticImplications_list):
            print("treatL=[],diagnosticImplications_list =[],treatL 没有,diagnostic 也没有，有药吗？")
            if oncogenic == 'Oncogenic' or oncogenic == 'Likely Oncogenic':
                pmids = ','.join(out['mutationEffect']['citations'].get('pmids','-'))
                alterations,location,cancerType,sub_cancerType,drug_name,level,FDAdescription = ('-','-','-','-','-','-','-')
                if anno_cosmic_list:
                    for i in range(len(anno_cosmic_list)):
                        chrom, pos, cosv, ref, alt, match_cosm, match_cds, match_AA, HGVSC, HGVSP, HGVSG,matchgene = anno_cosmic_list[i].split('$')
                        match_cds = change_hgvsc(hgvsc=match_cds)
                        print('\t'.join([HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType, "exonnum?",match_cds,match_AA, drug_name,'药物中文名',level,oncogenic,mutationEffect,FDAdescription, variantSummary, pmids, cosv,match_cosm,HGVSC, HGVSP, chrom, pos,ref, alt,matchgene,HGVSG,updatetime]), file=fh)

                else:
                    chrom, pos, cosv, ref, alt, match_cosm, match_cds, match_AA, HGVSC, HGVSP, HGVSG = ('-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-')
                    alterations, location, cancerType, sub_cancerType, drug_name, level, FDAdescription,matchgene = ('-', '-', '-', '-', '-', '-', '-','-')
                    print('\t'.join([HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType, "exonnum?",match_cds, match_AA, drug_name, '药物中文名', level, oncogenic, mutationEffect, FDAdescription,
                         variantSummary, pmids, cosv, match_cosm, HGVSC, HGVSP, chrom, pos, ref, alt, matchgene, HGVSG,
                         updatetime]), file=fh)
                #print('\t'.join([HugoSymbol, protein_change, alterations, location, cancerType, sub_cancerType, oncogenic,mutationEffect,drug_name,level, FDAdescription,variantSummary, pmids, updatetime]), file=fh)
        time.sleep(6)

    fh.close()


if  __name__ == "__main__":
    if len(sys.argv) != 3 :
        print("usage: python {} input_gene_protein_change_file ,outputfile".format(sys.argv[0]))

    else:
        main(sys.argv[1],sys.argv[2])




