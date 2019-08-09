'''
This code converts the .tsv output from neoepiscope to a .bed which can be uploaded as a in the genome browser
'''

import pandas as pd
import os
import sys

input_file = str(sys.argv[1])

df1 = pd.read_csv(input_file, sep='\t', header=1, dtype=object)
select_columns_for_bed = ['Neoepitope', 'Chromsome', 'Pos', 'Ref', 'Alt',
                          'Mutation_type', 'VAF', 'Paired_normal_epitope',
                          'Warnings', 'Transcript_ID']
lista = []

for i in df1.keys():
    if i not in select_columns_for_bed:
        lista.append(i)

df1.dropna(subset=lista, inplace=True)

df1 = df1[[float(x) <= 500 for x in df1[df1.keys()[10]]]]
# (df1.drop(df1[df1.index[11] >= 500]))

new = df1.filter(
    ['Neoepitope', 'Chromsome', 'Pos', 'Ref', 'Alt', 'Paired_normal_epitope'])
new['chromStart'] = new.apply(lambda x: int(x['Pos']) - 1, axis=1)

new['chromEnd'] = new.apply(lambda x: x['chromStart'] + len(x['Alt']),
                            axis=1)
#print(new)
new = new.drop("Pos", axis=1)
new = new[['Chromsome', 'chromStart', 'chromEnd', 'Neoepitope']]
new = new.rename(columns={'Chromsome': 'chrom'})
new.dropna(inplace=True)

new.to_csv("new.txt", sep="\t", index=False, header=0)
with open("new.txt") as f:
    g = open("new1.txt", "w+")
    g.write(
        "#browser position chr22:20100000-20100900 \n#track name=neoepitopes description=\"From neoepiscope comp.bed\" color=0,0,255 \n#chrom\tchromStart\tchromEnd\tNeoepitope \n")
    g.write(f.read())
    g.close()

f = open("new1.txt", "r")
contents = f.read()
print(contents)
