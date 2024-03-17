#!/usr/bin/env python
# Written by Niraj K. Nepal, Ph.D.
"""
Install matminer to extract composition based magpie features and structure based Jarvis features.
"""
import os
import zipfile
import glob
import json
import numpy as np
import pandas as pd
from matminer.featurizers.conversions import StrToComposition
from matminer.featurizers.composition import ElementProperty
from matminer.featurizers.structure import JarvisCFID
from sklearn.metrics import mean_absolute_error as mae
from sklearn import metrics as m
from mp_api.client import MPRester
from htepc import MpConnect
try:
    PWD = os.getcwd()
    if os.path.isfile(PWD+"/htepc.json"):
        JSONFILE = PWD+"/htepc.json"
    else:
        JSONFILE = "../../htepc.json"
    with open(JSONFILE, "r") as readjson:
        input_data = json.load(readjson)
except FileNotFoundError:
    print("htepc.json file not found\n")
class MlProcess:
    """
    This class processes output files from HTEPC package.
    Command "mainprogram.py 21" extracts the result file "result.csv"
    from EPC calculation and stores relaxed
    structure in .cif format inside "cif" folder.
    Useful in Machine Learning work, especially for CGCNN, and ALIGNN models.

    Parameters
    ----------
    in_file : str, optional
        The input .csv file with at least an ID (material id). Default is 'result.csv'.

    Attributes
    ----------
    data : pandas.DataFrame
        DataFrame containing the data loaded from the input file.
    prop : None
        Placeholder for the property. Not initialized.
    prop_list : list
        List to store properties.
    """
    def __init__(self,in_file='result.csv'):
        """
        Initialize the MlProcess object.

        Parameters
        ----------
        in_file : str, optional
            The input .csv file with at least an ID (material id). Default is 'result.csv'.
        """
        self.data = pd.DataFrame(pd.read_csv(in_file))
        self.prop = None
        self.prop_list = []
    def add_property(self,prop="energy_above_hull"):
        """
        Function to add a property in a new column.

        Parameters
        ----------
        prop : str, optional
            The name of the property as named in the Materials Project database.
            Default is "energy_above_hull".

        Returns
        -------
        tuple
            Shape of the updated data.

        Notes
        -----
        This method sets the specified property attribute, appends the property name to
        the property list, retrieves property data from Materials Project database
        using MpConnect, and adds the property data as a new column to the DataFrame.

        If the property is "symmetry", it retrieves the symmetry number from the Materials Project.
        For other properties, it retrieves the property value using MpConnect.

        """
        self.prop = prop
        self.prop_list.append(prop)
        obj = MpConnect()
        ydata = []
        for i in range(self.data.shape[0]):
            obj.setting(self.data["ID"][i])
            if prop == 'symmetry':
                ydata.append(obj.data['symmetry']['symbol'])
            else:
                ydata.append(obj.property(prop))
        self.data[self.prop] = ydata
        return self.data.shape
    def magpie_feature(self,outfile="magpie_featurize.csv"):
        """
        Function to add composition-based Magpie features from the matminer package.

        Parameters
        ----------
        outfile : str, optional
            The name of the .csv file to save the dataframe with columns containing Magpie features.
            Default is "magpie_featurize.csv".

        Notes
        -----
        This method featurizes the composition data in the dataframe using Magpie features
        from the matminer package. It utilizes the StrToComposition class to convert
        composition strings into Composition objects, and then adds Magpie features
        using ElementProperty.from_preset.

        The resulting dataframe with Magpie features is saved to the specified outfile.
        """
        df_i = StrToComposition().featurize_dataframe(self.data, "compound")
        ep_feat = ElementProperty.from_preset(preset_name="magpie")
        df_i = ep_feat.featurize_dataframe(df_i, col_id="composition")
        df_i.to_csv(outfile,index=False)
    def jarvis_structure_feature(self,outfile="jarvis_featurize.csv"):
        """
        Function to add structure-based Jarvis features from the matminer package.

        Parameters
        ----------
        outfile : str, optional
            The name of the .csv file to save the dataframe with columns containing Jarvis features.
            Default is "jarvis_featurize.csv".

        Notes
        -----
        This method retrieves structures using Materials Project API based on material IDs.
        It then adds Jarvis features using JarvisCFID from the matminer package.

        If the file 'htepc.json' or its relative path exists, it assumes it contains an API key
        for the Materials Project database. It retrieves the key and uses it to authenticate with
        MPRester. If the file does not exist, it prints a message indicating that 'htepc.json' is
        not found.

        The resulting dataframe with Jarvis features is saved to the specified outfile.
        """
        if os.path.isfile("htepc.json") or os.path.isfile("../../htepc.json"):
            key = input_data['mpi_key']['API_KEY']
            mpr = MPRester(key['key'])
        else:
            print("htepc.json not found")
            #print(" Please, replace the XXX with your key\n")
            #with open("mpi_key.py", "w") as k:
            #    api = {'key':"XXX"}
            #    k.write("API_KEY={}".format(api) + "\n")
        structure = []
        for i in range(self.data.shape[0]):
            mpid = self.data['ID'][i]
            structure.append(mpr.get_structure_by_material_id(mpid))
        self.data['structure'] = structure
        jcfid = JarvisCFID()
        df_i = jcfid.featurize_dataframe(self.data, 'structure')
        df_i = df_i.drop(['structure'], axis=1)
        df_i.to_csv(outfile, index=False)
    def class_convert(self,prop='Tc',cutoff=1.0):
        """
        Function to change a property into a label (0 or 1) used to train and test CGCNN classification models.

        Parameters
        ----------
        prop : str, optional
            The property to be used for labeling. Default is 'Tc'.
        cutoff : float, optional
            The cutoff value to separate data into binary labels. Default is 1.0.

        Notes
        -----
        This method creates a new column 'target' in the dataframe and assigns binary labels (0 or 1)
        based on whether the property value in the specified column 'prop' is greater than the cutoff value.

        """
        dcol = self.data[prop] > cutoff
        dcol = dcol.astype('int')
        self.data['target'] = dcol
    def write_output(self,columns,outfile="result-new.csv"):
        """
        Function to write new file
        parameters
        ---------------------
        columns : (list of str) columns to be extracted.
        outfile : .csv file
        Returns
        ----------------------
        shape of current state of data
        """
        self.data[columns].to_csv(outfile, index=False)
        return self.data.shape
    def write_id_prop_csv(self,mode='alignn',prop='Tc',outfile="id_prop.csv"):
        """
        function to write id_prop.csv file required to work with CGCNN and ALIGNN models.
        parameters
        -------------------------
        mode : (str) name of model. Default: 'alignn' or 'ALIGNN'
        prop : (float) property
        outfile : id_prop.csv file
        """
        new_data = pd.DataFrame(columns=["ID",prop])
        if mode in ('alignn', 'ALIGNN'):
            new_data['ID'] = self.data['ID'] + '.cif'
        else:
            new_data['ID'] = self.data['ID']
        new_data[prop] = self.data[prop]

        new_data.to_csv(outfile,index=False)
        os.system("sed -i '1d' {}".format(outfile))
def ml_tc(lam='lambda.csv',wlog='wlog.csv',calc_tc='Tc_calc.csv'):
    """
    Function to compute critical temperature Tc from ML trained lambda and ML trained wlog
    parameters
    ------------------
    lam : .csv file for lambda
    wlog : .csv file for wlog
    calc_tc : .csv file to store calculated Tc
    Each .csv file has 3 columns, ID,target,prediction.
    "ID.cif" be the structure file for particular "ID"
    """
    lam = lam.set_index('ID')
    lam = lam.reindex(index=wlog['ID'])
    lam = lam.reset_index()
    data1 = pd.DataFrame(columns=["ID","target","prediction"])
    data1['target'] = np.round((wlog.target/1.2) * np.exp(-1.04*(1.0 + lam.target)/np.abs(lam.target - (1.0 + 0.62*lam.target)*0.16)),5)
    data1['prediction'] = np.round((wlog.prediction/1.2) * np.exp(-1.04*(1.0 + lam.prediction)/np.abs(lam.prediction - (1.0 + 0.62*lam.prediction)*0.16)),5)
    data1['ID'] = lam['ID']
    data1.to_csv(calc_tc,index=False)
def mae_compute(outfile="mae.in"):
    """
    function to calculate mean-absolute-error (mae).
    parameters
    ------------------
    outfile : .in file to store mae values for each file
    """
    print("Make sure, you have .csv files with ID,target, and prediction columns\n")
    filename = glob.glob("*.csv")
    with open(outfile, "w") as gfile:
        for filei in filename:
            data = pd.DataFrame(pd.read_csv(filei))
            gfile.write(filei + ": MAE {}".format(mae(data.target, data.prediction))+ "\n")
def class_accuracy(outfile="classification_score.in"):
    """
    Function to calculate accuracy measures for classification models.
    Accuracy: TP+TN/(TP+TN+FP+FN)
    Precision: TP/(TP+FP)
    Recall : TP/(TP+FN)
    F1-score : 2*(Precision*Recall)/(Precision+Recall)
    parameters
    ------------------
    outfile : .in file to store mae values for each file
    """
    print("Make sure, you have .csv files with target, and prediction columns\n")
    filename = glob.glob("*.csv")
    with open(outfile, "w") as gfile:
        for filei in filename:
            data = pd.DataFrame(pd.read_csv(filei))
            recall = round(m.recall_score(data.target,data.prediction),2)
            precision = round(m.precision_score(data.target,data.prediction),2)
            acc = round(m.accuracy_score(data.target,data.prediction),2)
            f1_score = round(m.f1_score(data.target,data.prediction),2)
            gfile.write("file: {}, accuracy_score: {}, recall_score: {}, precision_score: {}, f1_score: {}".format(filei,acc,precision,recall,f1_score)+ "\n")
def zip_alignn_output(folder,outfile,epochs=3000):
    """
    Function to compress output file of ALIGNN model to use as model for new prediction.
    parameters
    -----------------
    folder : output folder to compress
    outfile : name of .zip file to be used as a ALIGNN model.
    epochs : (int) number of epochs for training
    """
    files = [folder+'/checkpoint_{}.pt'.format(epochs-1),folder+'/checkpoint_{}.pt'.format(epochs),folder+'/config.json',folder+'/history_train.json',folder+'/history_val.json',folder+'/ids_train_val_test.json',folder+'/mad',folder+'/prediction_results_test_set.csv',folder+'/prediction_results_train_set.csv',folder+'/test_data_data_range',folder+'/train_data_data_range',folder+'/val_data_data_range']
    files.append(folder)
    zout = zipfile.ZipFile(outfile+".zip", "w") # <--- this is the change you need to make
    for fname in files:
        zout.write(fname)
    zout.close()
#if __name__ == "__main__":
