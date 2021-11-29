#!/usr/bin/env python3
# -*- coding: utf-8 -*
# @Author  : yellow_huang

import subprocess

def getbyProteinChange(genename,protein_change,API_Token='90ec77a3-cfc6-44b7-9cd8-abd952a3749b'):
    cmd = (
        f'curl -X GET "https://oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol={genename}&alteration={protein_change}" '
        f'-H "Authorization: Bearer {API_Token}" '
        f'-H "accept: application/json" '
    )
    print(cmd)
    out = subprocess.check_output(cmd,shell=True)
    return  out



if __name__ == '__main__':
    #out = getbyProteinChange("KIT","552_558del")
    out = getbyProteinChange("PTEN","Y176X")
    print(out)