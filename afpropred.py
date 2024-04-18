##############################################################################
#AfProPred is developed for predicting AFP and non-AFP      #
#peptides from their primary sequence. It is developed by Prof G. P. S.       #
#Raghava's group. Please cite : ToxinPred 3.0                                  #
# ############################################################################
import argparse  
import warnings
import pickle
import os
import re
import sys
import pandas as pd
import numpy as np
import glob
from Bio import SeqIO
import shutil
import requests
import zipfile


warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments.') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-o", "--output",type=str, default="outfile.csv", help="Output: File for saving results by default outfile.csv")
parser.add_argument("-t","--threshold", type=float, default=0.48, help="Threshold: Value between 0 to 1 by default 0.48")
parser.add_argument("-m","--model",type=int, default=2, choices = [1, 2], help="Model: 1: AAC feature based ExtraTrees Classifier , 2: AAC + PSSM feature based ExtraTrees Classifier, by default 1")
parser.add_argument("-d","--display", type=int, choices = [1,2], default=2, help="Display: 1:AFP, 2: All peptides, by default 2")
parser.add_argument("-wd", "--working",type=str, default='.', help="Working Directory: Temporary directory to write files")

args = parser.parse_args()


nf_path = os.path.dirname(__file__)
std = list("ACDEFGHIKLMNPQRSTVWY")


if not os.path.exists(nf_path + '/model'):
    response = requests.get('https://webs.iiitd.edu.in/raghava/afpropred/model.zip')
    if response.status_code == 200:
        with open(os.path.join(nf_path + '/model.zip'), 'wb') as f:
                f.write(response.content)
    with zipfile.ZipFile(nf_path + '/model.zip', 'r') as zip_ref:
            zip_ref.extractall(nf_path + '/' )
    print("ZIP file contents extracted successfully.")
    os.remove(nf_path + '/model.zip')
    
if not os.path.exists(nf_path + '/swissprot'):
    response = requests.get('https://webs.iiitd.edu.in/raghava/afpropred/swissprot.zip')
    if response.status_code == 200:
        with open(os.path.join(nf_path + '/swissprot.zip'), 'wb') as f:
                f.write(response.content)
    with zipfile.ZipFile(nf_path + '/swissprot.zip', 'r') as zip_ref:
            zip_ref.extractall(nf_path + '/' )
    print("ZIP file contents extracted successfully.")
    os.remove(nf_path + '/swissprot.zip')
    

def fasta_to_dataframe(fasta_file):
    sequences = {'ID': [], 'Sequence': []}
    
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequence_id = record.id.lstrip('>')  # Remove '>' from the ID
        sequences['ID'].append(sequence_id)
        sequences['Sequence'].append(str(record.seq))
    
    # Convert to DataFrame
    df = pd.DataFrame(sequences)
    
    return df

def aac_comp(file,out):
    filename, file_extension = os.path.splitext(file)
    f = open(out, 'w')
    sys.stdout = f
    df = fasta_to_dataframe(file)
    zz = df['Sequence']
    print("AAC_A,AAC_C,AAC_D,AAC_E,AAC_F,AAC_G,AAC_H,AAC_I,AAC_K,AAC_L,AAC_M,AAC_N,AAC_P,AAC_Q,AAC_R,AAC_S,AAC_T,AAC_V,AAC_W,AAC_Y,")
    for j in zz:
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            print("%.2f"%composition, end = ",")
        print("")
    f.close()
    sys.stdout = sys.__stdout__




def readseq(file):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    seqid = []
    seq = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '', ''.join(array[1:]).upper())
        seqid.append('>'+name)
        seq.append(sequence)
    if len(seqid) == 0:
        f=open(file,"r")
        data1 = f.readlines()
        for each in data1:
            seq.append(each.replace('\n',''))
        for i in range (1,len(seq)+1):
            seqid.append(">Seq_"+str(i))

    final_df = pd.concat([pd.DataFrame(seqid),pd.DataFrame(seq)], axis=1)
    final_df.columns = ['ID','Seq']
    return final_df

def file_split(file,path):
    df1 = readseq(file)
    for i in range(len(df1)):
        df1.loc[i].to_csv(path+'/'+df1['ID'][i].replace('>','')+'.fasta', index=None,header=False,sep="\n")

def gen_pssm(fasta_path, pssm_path):
    os.makedirs(pssm_path+'/pssm_raw1', exist_ok = True)
    os.makedirs(pssm_path+'/pssm_raw', exist_ok = True)

    listdir = glob.glob(fasta_path+'/*.fasta')
    for i in listdir:
        filename = i.split('/')[-1].rsplit('.', 1)[0]
        cmd = nf_path + "/ncbi_blast_2.15/bin/psiblast -out "+pssm_path+"/pssm_raw1/"+filename+".homologs -outfmt 7 -query "+fasta_path+"/"+filename+".fasta -db "+nf_path+"/swissprot/swissprot -evalue 0.001 -word_size 3 -max_target_seqs 6000 -num_threads 10 -gapopen 11 -gapextend 1 -matrix BLOSUM62 -comp_based_stats T -num_iterations 3 -out_pssm "+pssm_path+"/pssm_raw1/"+filename+".cptpssm -out_ascii_pssm "+pssm_path+"/pssm_raw/"+filename+".pssm"
        os.system(cmd)

def feat_gen_aac(file,out):
    aac_comp(file, wd + '/aac_temp')
    df = pd.read_csv(wd + '/aac_temp')
    df = df.iloc[:,:-1]
    df.to_csv(out, index=None)
    os.remove(wd + '/aac_temp')

def feat_gen_acc_pssm(file, out1, out2):
    feat_gen_aac(file, wd + '/aac_temp1')
    aac_df = pd.read_csv(wd + '/aac_temp1')
    file_split(file, wd + '/fasta_files')
    gen_pssm(wd + '/fasta_files', wd + '/pssm_files')
    df = fasta_to_dataframe(file)
    aac_df['ID'] = df['ID']
    aac_df.to_csv(out1, index=None)
    folder_path = wd + '/pssm_files/pssm_raw/'
    # List all items in the folder
    folder_items = os.listdir(folder_path)    
    df1 = pd.DataFrame({'ID': folder_items})
    df1['ID'] = df1['ID'].str.replace('.pssm', '')
    pssm1 = pd.merge(df1, df, on="ID", how="left")
    pssm1['ID'] = '>' + pssm1['ID']
    pssm1[['ID','Sequence']].to_csv(wd + '/pssm_input.fasta' ,index= None, header=False, sep="\n")
    os.system('python ' + nf_path + '/possum/possum.py -i ' + wd + '/pssm_input.fasta' + ' -o ' + wd + '/pssm_temp1 -t pssm_composition -p ' + wd + '/pssm_files/pssm_raw')
    os.system('python ' + nf_path + '/possum/headerHandler.py -i ' + wd + '/pssm_temp1' + ' -o ' + wd + '/pssm_temp2 -p pssm_')
    df2 = pd.read_csv(wd + '/pssm_temp2')

    df2['ID'] = df1['ID']
    df3 = pd.merge(aac_df, df2, on="ID", how="right")
    columns = list(df3.columns)
    columns.remove('ID')
    columns.append('ID')
    df3 = df3[columns]
    df3.to_csv(out2)
    os.remove(wd + '/aac_temp1')
    os.remove(wd + '/pssm_input.fasta')
    os.remove(wd + '/pssm_temp1')
    os.remove(wd + '/pssm_temp2')



def prediction_aac(inputfile1, model,out):
    a=[]
    file_name = inputfile1
    file_name1 = out
    file_name2 = model
    with open(file_name2, 'rb') as file:
        clf1 = pickle.load(file)      
    data_test1 = pd.read_csv(inputfile1)    
    X_test = data_test1
    y_p_score1=clf1.predict_proba(X_test)
    y_p_s1=y_p_score1.tolist()
    df = pd.DataFrame(y_p_s1)
    df_1 = df.iloc[:,-1]
    df_1.to_csv(file_name1, index=None, header=False)

def prediction_aac_pssm(inputfile1, inputfile2, model1, model2, out1, out2):
    a=[]
    with open(model1, 'rb') as file:
        clf1 = pickle.load(file)    
    with open(model2, 'rb') as file:
        clf2 = pickle.load(file)   
     
    data_test1 = pd.read_csv(inputfile1)
    data_test2 = pd.read_csv(inputfile2, index_col=0)

    X_test1 = np.array(data_test1.iloc[:,:-1])
    X_test2 = np.array(data_test2.iloc[:,:-1])
    y_p_score1=clf1.predict_proba(X_test1)
    y_p_score2=clf2.predict_proba(X_test2)

    y_p_s1=y_p_score1.tolist()
    y_p_s2=y_p_score2.tolist()

    df_dict1 = {'ML Score' : pd.DataFrame(y_p_s1).iloc[:,1], 'ID' : list(data_test1.iloc[:,-1])}
    df_1 = pd.DataFrame(data=df_dict1)
    df_1.to_csv(out1, index=None)

    df_dict2 = {'ML Score' : pd.DataFrame(y_p_s2).iloc[:,1], 'ID' : list(data_test2.iloc[:,-1])}
    df_2 = pd.DataFrame(data=df_dict2)
    df_2.to_csv(out2, index=None)

def class_assignment_aac(file1, thr, out):
    df1 = pd.read_csv(file1, header=None)
    df1.columns = ['ML Score']
    cc = []
    for i in range(0,len(df1)):
        if df1['ML Score'][i]>=float(thr):
            cc.append('AFP')
        else:
            cc.append('Non-AFP')
    df1['Prediction'] = cc
    df1 =  df1.round(3)
    df1.to_csv(out, index=None)

def class_assignment_aac_pssm(file1, file2, thr, out):
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    merged_df = pd.merge(df1, df2, on='ID', suffixes=('_aac', '_pssm'), how='left')
    merged_df['ML Score_aac'] = merged_df['ML Score_pssm'].fillna(merged_df['ML Score_aac'])
    merged_df.drop(columns=['ID', 'ML Score_pssm'], inplace=True)
    merged_df.rename(columns={'ML Score_aac': 'ML Score'}, inplace=True)

    cc = []
    for i in range(0,len(merged_df)):
        if merged_df['ML Score'][i]>=float(thr):
            cc.append('AFP')
        else:
            cc.append('Non-AFP')
    df1['Prediction'] = cc
    df1.drop(columns=['ID'], inplace=True)
    df1 =  df1.round(3)
    df1.to_csv(out, index=None)

print('##############################################################################')
print('# The program AfProPred is developed for predicting AFP and non-AFP #')
print("# from their primary sequence, developed by Prof G. P. S. Raghava's group. #")
print('# ############################################################################')

# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 
 
# Output file 
result_filename = args.output
         
# Threshold 
Threshold= float(args.threshold)

# Model
Model = int(args.model)

# Display
dplay = int(args.display)

# Working Directory
wd = args.working

print('Summary of Parameters:')
print('Input File: ',Sequence,'; Model: ',Model,'; Threshold: ', Threshold)
print('Output File: ',result_filename,'; Display: ',dplay)

#------------------ Read input file ---------------------
f=open(Sequence,"r")
len1 = f.read().count('>')
f.close()

with open(Sequence) as f:
        records = f.read()
records = records.split('>')[1:]
seqid = []
seq = []
for fasta in records:
    array = fasta.split('\n')
    name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(array[1:]).upper())
    seqid.append(name)
    seq.append(sequence)
if len(seqid) == 0:
    f=open(Sequence,"r")
    data1 = f.readlines()
    for each in data1:
        seq.append(each.replace('\n',''))
    for i in range (1,len(seq)+1):
        seqid.append("Seq_"+str(i))

seqid_1 = list(map(">{}".format, seqid))

#======================= Prediction Module start from here =====================
if Model==1:
    feat_gen_aac(Sequence, wd + '/seq.aac')
    prediction_aac(wd + '/seq.aac', nf_path + '/model/model_aac', wd + '/seq.pred')
    class_assignment_aac(wd +'/seq.pred',Threshold, wd + '/seq.out')
    df1 = pd.DataFrame(seqid)
    df2 = pd.DataFrame(seq)
    df3 = pd.read_csv(wd + "/seq.out")
    df3 = round(df3,3)
    df4 = pd.concat([df1,df2,df3],axis=1)
    df4.columns = ['ID','Sequence','ML Score','Prediction']
    df4.loc[df4['ML Score'] > 1, 'ML Score'] = 1
    df4.loc[df4['ML Score'] < 0, 'ML Score'] = 0
    if dplay == 1:
        df4 = df4.loc[df4.Prediction=="AFP"]
    df4.to_csv(result_filename, index=None)
    os.remove(wd + '/seq.aac')
    os.remove(wd + '/seq.pred')
    os.remove(wd + '/seq.out')
else:
    os.makedirs(wd + '/fasta_files', exist_ok=True)
    os.makedirs(wd + '/pssm_files', exist_ok=True)
    feat_gen_acc_pssm(Sequence, wd + '/seq.aac', wd + '/seq.pssm')
    prediction_aac_pssm(wd + '/seq.aac', wd + '/seq.pssm', nf_path + '/model/model_aac', nf_path + '/model/model_pssm', wd + '/seq_aac.pred', wd + '/seq_pssm.pred')
    class_assignment_aac_pssm(wd +'/seq_aac.pred', wd +'/seq_pssm.pred', Threshold, wd + '/seq.out')
    df1 = pd.DataFrame(seqid)
    df2 = pd.DataFrame(seq)
    df3 = pd.read_csv(wd + "/seq.out")
    df3 = round(df3,3)
    df4 = pd.concat([df1,df2,df3],axis=1)
    df4.columns = ['ID','Sequence','ML Score','Prediction']
    df4.loc[df4['ML Score'] > 1, 'ML Score'] = 1
    df4.loc[df4['ML Score'] < 0, 'ML Score'] = 0
    if dplay == 1:
        df4 = df4.loc[df4.Prediction=="AFP"]
    df4.to_csv(result_filename, index=None)
    os.remove(wd + '/seq.aac')
    os.remove(wd + '/seq.pssm')
    os.remove(wd + '/seq_aac.pred')
    os.remove(wd + '/seq_pssm.pred')
    os.remove(wd + '/seq.out')
    shutil.rmtree(wd + '/fasta_files')
    shutil.rmtree(wd + '/pssm_files')

print('\n======= Thanks for using AfProPred. Your results are stored in file :',result_filename,' =====\n\n')
print('Please cite: AfProPred\n')
