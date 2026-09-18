#!/usr/bin/env python
# Written by Niraj K. Nepal, Ph.D.
"""
Install matminer to extract composition based magpie features and structure based Jarvis features.

FIX(15): ``matminer`` and ``scikit-learn`` are optional extras, but they used
to be imported at module scope, so ``import htesp.ml_processing`` -- and any
package-wide import that reached it -- failed with ImportError on an
installation that only has the core dependencies.  They are now imported
inside the functions that use them, behind :func:`_require`, which explains
how to install them.
"""
import warnings
import zipfile
import glob
import numpy as np
import pandas as pd
from htesp.htepc import MpConnect, mprester
from htesp.check_json import config, require_api_key

#: Allen-Dynes Coulomb pseudopotential used by :func:`ml_tc`
MUSTAR = 0.16


def _require(module_name, extra="ml"):
    """Import ``module_name`` or raise a message naming the extra to install.

    FIX(15).
    """
    import importlib
    try:
        return importlib.import_module(module_name)
    except ImportError as exc:
        raise ImportError(
            "{} is required for this function but is not installed.\n"
            "    pip install htesp[{}]".format(module_name.split(".")[0], extra)
        ) from exc
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
        # FIX(15): optional dependency imported at the point of use
        conversions = _require("matminer.featurizers.conversions")
        composition = _require("matminer.featurizers.composition")
        df_i = conversions.StrToComposition().featurize_dataframe(self.data, "compound")
        ep_feat = composition.ElementProperty.from_preset(preset_name="magpie")
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

        If the file 'config.json' or its relative path exists, it assumes it contains an API key
        for the Materials Project database. It retrieves the key and uses it to authenticate with
        MPRester. If the file does not exist, it prints a message indicating that 'config.json' is
        not found.

        The resulting dataframe with Jarvis features is saved to the specified outfile.
        """
        # FIX(15): optional dependency imported at the point of use
        structure_featurizers = _require("matminer.featurizers.structure")
        # The key comes from $MP_API_KEY / ~/.config/htesp/credentials /
        # config.json, and require_api_key() raises a message that says how to
        # set it rather than leaving ``mpr`` unbound.
        mpr = mprester(require_api_key(config()))
        structure = []
        for i in range(self.data.shape[0]):
            mpid = self.data['ID'][i]
            structure.append(mpr.get_structure_by_material_id(mpid))
        self.data['structure'] = structure
        jcfid = structure_featurizers.JarvisCFID()
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

        # FIX(all): os.system("sed -i '1d' ...") replaced by a pure-Python
        # header strip, so ``outfile`` is never handed to a shell.
        new_data.to_csv(outfile,index=False)
        with open(outfile, "r") as read_idprop:
            body = read_idprop.readlines()
        with open(outfile, "w") as write_idprop:
            write_idprop.writelines(body[1:])
def allen_dynes(wlog, lam, mustar=MUSTAR):
    """Allen-Dynes Tc for array-like ``wlog`` (K) and ``lam``.

    FIX(14): the denominator ``lam - mustar*(1 + 0.62*lam)`` used to be wrapped
    in ``np.abs()``.  That turns the pole at lambda ~= 0.178 into a mirror and
    manufactures a finite, *positive* Tc for every lambda below it, where the
    Allen-Dynes formula in fact has no superconducting solution.  The
    denominator is used as-is and a non-positive one yields NaN, which is the
    physically correct answer, with a warning.
    """
    wlog = np.asarray(wlog, dtype=float)
    lam = np.asarray(lam, dtype=float)
    denom = lam - (1.0 + 0.62 * lam) * mustar
    bad = ~(denom > 0.0)
    if np.any(bad):
        warnings.warn(
            "Allen-Dynes denominator lambda - mu*(1 + 0.62 lambda) is "
            "non-positive for {} of {} entries (lambda below ~{:.3f} for "
            "mu*={}); Tc is NaN there, not a mirrored finite value.".format(
                int(np.count_nonzero(bad)), lam.size,
                mustar / (1.0 - 0.62 * mustar), mustar),
            RuntimeWarning, stacklevel=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_c = (wlog / 1.2) * np.exp(-1.04 * (1.0 + lam) / denom)
    return np.where(bad, np.nan, t_c)
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
    # FIX(14): np.abs() around the Allen-Dynes denominator removed
    data1['target'] = np.round(allen_dynes(wlog.target, lam.target),5)
    data1['prediction'] = np.round(allen_dynes(wlog.prediction, lam.prediction),5)
    data1['ID'] = lam['ID']
    data1.to_csv(calc_tc,index=False)
def mae_compute(outfile="mae.in"):
    """
    function to calculate mean-absolute-error (mae).
    parameters
    ------------------
    outfile : .in file to store mae values for each file
    """
    # FIX(15): optional dependency imported at the point of use
    mae = _require("sklearn.metrics").mean_absolute_error
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
    # FIX(15): optional dependency imported at the point of use
    m = _require("sklearn.metrics")
    print("Make sure, you have .csv files with target, and prediction columns\n")
    filename = glob.glob("*.csv")
    with open(outfile, "w") as gfile:
        for filei in filename:
            data = pd.DataFrame(pd.read_csv(filei))
            recall = round(m.recall_score(data.target,data.prediction),2)
            precision = round(m.precision_score(data.target,data.prediction),2)
            acc = round(m.accuracy_score(data.target,data.prediction),2)
            f1_score = round(m.f1_score(data.target,data.prediction),2)
            # FIX(13): the recall_score and precision_score labels were
            # swapped -- the format arguments were (acc, precision, recall)
            # against the labels (accuracy, recall, precision).
            gfile.write("file: {}, accuracy_score: {}, recall_score: {}, precision_score: {}, f1_score: {}".format(filei,acc,recall,precision,f1_score)+ "\n")
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
