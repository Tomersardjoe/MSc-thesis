"""
rf_module.py

Module used for the analysis of pangenomes using random forests.
"""
import sys
import re
import ast
import random
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report

def read_fasta(filename):
    """Reads a fasta file returning the seqs as a dictionary."""
    import re
    import collections
    lines = get_file_data(filename)
    seqs = collections.OrderedDict()
    key = ""
    value = ""
    for line in lines:
        line = line.rstrip("\n")
        if re.search(">", line):
            if key:
                seqs[key] = value
                key = line[1:]
            else:
                key = line[1:]
            value = ""
        else:
            value = value + line
    seqs[key] = value
    return(seqs)


def get_file_data(filename):
    """Stores the lines of the program with name filename as a list."""
    with open(filename, encoding = "utf8") as in_file:
        lines = []
        for line in in_file:
            lines.append(line.rstrip("\n"))
    return lines

def init_tables(table):
    """Initialise the performance and importance tables."""
    n_g = table.shape[0] #number of genes (gene families)
    imp = pd.DataFrame(0.0, index = np.arange(n_g),
                       columns = table.index.values)
    imp.index = table.index.values
    metrics = ['count', 'TPte', 'FPte', 'FNte', 'TNte', 'Ete', 'Ate', 'P1te',
               'P0te', 'Pte', 'R1te', 'R0te', 'Rte', 'F1te', 'F0te', 'Fte',
               'TPtr','FPtr', 'FNtr', 'TNtr', 'Etr', 'Atr', 'P1tr', 'P0tr',
               'Ptr', 'R1tr', 'R0tr', 'Rtr', 'F1tr', 'F0tr', 'Ftr']
    performance = pd.DataFrame(0.0, index = np.arange(n_g), columns = metrics)
    performance.index = table.index.values

    # Normalize indices/columns to consistent flattened labels
    imp.index = flatten_index_like(imp.index)
    imp.columns = flatten_index_like(imp.columns)
    performance.index = flatten_index_like(performance.index)

    return imp, performance

def preprocess_df(table, null_h, min_missing, min_present):
    """Modify the matrix so it's ready for random forest."""
    #new_rows = []
    #for row in table.index:
    #    new_row = []
    #    for field in row:
    #        new_row.append(str(field))
    #    new_rows.append(",".join(new_row))
    #table.index = new_rows #convert row header to strings.
    table = table.fillna(0) #replace absent with 0
    table = table.replace({'.{2,}': '1'}, regex=True) #convert to 1
    table = table.astype(int) #convert to numeric
    #remove genes with too much present
    table = table[table.sum(axis=1) < table.shape[1] - (min_missing -1)]
    #remove genes with too much absent
    table = table[table.sum(axis=1) > min_present - 1]
    if null_h:
        print("Randomising!!!!")
        for i in range(table.shape[0]):
            table.iloc[i] = random.sample(list(table.iloc[i]),
                                       len(list(table.iloc[i])))
    return table

def update_performance(table, i, y_sets):
    """Update the performance table."""
    
    # Number of genomes present
    table.iloc[i, table.columns.get_loc('count')] = sum(list(y_sets[0]))

    cm_train = confusion_matrix(y_sets[1], y_sets[3])
    cm_test = confusion_matrix(y_sets[2], y_sets[4])

    # Training metrics
    table.iloc[i, table.columns.get_loc('TPtr')] = cm_train[1,1]
    table.iloc[i, table.columns.get_loc('FPtr')] = cm_train[0,1]
    table.iloc[i, table.columns.get_loc('TNtr')] = cm_train[0,0]
    table.iloc[i, table.columns.get_loc('FNtr')] = cm_train[1,0]

    # Testing metrics
    table.iloc[i, table.columns.get_loc('TPte')] = cm_test[1,1]
    table.iloc[i, table.columns.get_loc('FPte')] = cm_test[0,1]
    table.iloc[i, table.columns.get_loc('TNte')] = cm_test[0,0]
    table.iloc[i, table.columns.get_loc('FNte')] = cm_test[1,0]

    # Reports for recall, precision, f1, accuracy
    train_report = classification_report(y_sets[1], y_sets[3], output_dict=True, zero_division=0)
    test_report = classification_report(y_sets[2], y_sets[4], output_dict=True, zero_division=0)

    # Accuracy and error
    table.iloc[i, table.columns.get_loc('Etr')] = 1 - train_report['accuracy']
    table.iloc[i, table.columns.get_loc('Ete')] = 1 - test_report['accuracy']
    table.iloc[i, table.columns.get_loc('Atr')] = train_report['accuracy']
    table.iloc[i, table.columns.get_loc('Ate')] = test_report['accuracy']

    # Precision
    table.iloc[i, table.columns.get_loc('P1tr')] = train_report['1']['precision']
    table.iloc[i, table.columns.get_loc('P0tr')] = train_report['0']['precision']
    table.iloc[i, table.columns.get_loc('Ptr')] = train_report['macro avg']['precision']
    table.iloc[i, table.columns.get_loc('P1te')] = test_report['1']['precision']
    table.iloc[i, table.columns.get_loc('P0te')] = test_report['0']['precision']
    table.iloc[i, table.columns.get_loc('Pte')] = test_report['macro avg']['precision']

    # Recall
    table.iloc[i, table.columns.get_loc('R1tr')] = train_report['1']['recall']
    table.iloc[i, table.columns.get_loc('R0tr')] = train_report['0']['recall']
    table.iloc[i, table.columns.get_loc('Rtr')] = train_report['macro avg']['recall']
    table.iloc[i, table.columns.get_loc('R1te')] = test_report['1']['recall']
    table.iloc[i, table.columns.get_loc('R0te')] = test_report['0']['recall']
    table.iloc[i, table.columns.get_loc('Rte')] = test_report['macro avg']['recall']

    # F1
    table.iloc[i, table.columns.get_loc('F1tr')] = train_report['1']['f1-score']
    table.iloc[i, table.columns.get_loc('F0tr')] = train_report['0']['f1-score']
    table.iloc[i, table.columns.get_loc('Ftr')] = train_report['macro avg']['f1-score']
    table.iloc[i, table.columns.get_loc('F1te')] = test_report['1']['f1-score']
    table.iloc[i, table.columns.get_loc('F0te')] = test_report['0']['f1-score']
    table.iloc[i, table.columns.get_loc('Fte')] = test_report['macro avg']['f1-score']

def flatten_index_like(obj):
    def clean_entry(entry):
        # If it's a tuple, prefer the first non-empty element
        if isinstance(entry, tuple):
            for part in entry:
                if pd.notna(part) and str(part) != "":
                    return str(part)
            return "_".join(str(i) for i in entry if pd.notna(i))

        # If it's a string, extract a compact token (letters, digits, _ or -)
        if isinstance(entry, str):
            match = re.search(r"([A-Za-z0-9\-_]+)", entry)
            if match:
                return match.group(1)
            return entry

        # Fallback
        return str(entry)

    return [clean_entry(e) for e in obj]

def fit_classifiers(table, results, params, output, checkpoint):
    """
    Fit a random forest classifier for all genes.
    results = [imp, performance]
    params = [ntrees, depth, purity, nthreads]
    """
    n_g = table.shape[0]
    gene_labels = flatten_index_like(table.index.values)
    table = table.transpose()
    start = 0
    if checkpoint != 0:
        results[0] = pd.read_csv(output + "/imp.csv", header=0, index_col=0)
        results[1] = pd.read_csv(output + "/performance.csv", header=0, index_col=0)
        # Normalize immediately after loading
        results[0].index = flatten_index_like(results[0].index)
        results[0].columns = flatten_index_like(results[0].columns)
        results[1].index = flatten_index_like(results[1].index)
        results[1] = results[1].reindex(gene_labels)

        start = checkpoint

    for i in range(start, n_g):
        print("gene number\t" + str(i+1) + "\tout of\t" + str(n_g))
        y_all = table[table.columns[i]]
        x_all = table.drop([table.columns[i]], axis = 1)
        #now split the dataset into test and train - train with 75% in each
        #class
        # datasets = x_train, x_test, y_train, y_test
        datasets = train_test_split(x_all, y_all,
                                    test_size=0.25, stratify=y_all)
        #random forest
        if params[2]:
            if params[1]:
                model = RandomForestClassifier(n_estimators = params[0],
                                               min_impurity_decrease =\
                                                       params[2],
                                               max_depth = params[1],
                                               max_features='sqrt',
                                               min_samples_split=2,
                                               n_jobs=params[3])
            else:
                model = RandomForestClassifier(n_estimators = params[0],
                                               min_impurity_decrease =\
                                                       params[2],
                                               max_features='sqrt',
                                               min_samples_split=2,
                                               n_jobs=params[3])
 
        elif params[1]:
            model = RandomForestClassifier(n_estimators = params[0],
                                           max_depth = params[1],
                                           max_features='sqrt',
                                           min_samples_split=2,
                                           n_jobs=params[3])

            
        else:
            print("either depth or purity (or both) must be set")
            sys.exit()
        model.fit(datasets[0], datasets[2])
        y_pred_train = model.predict(datasets[0])
        y_pred_test = model.predict(datasets[1])
        #assess performance
        update_performance(results[1], i,
                           [y_all, datasets[2], datasets[3],
                            y_pred_train, y_pred_test])
        #update importances
        results[0][results[0].columns[i]] = \
                np.insert(model.feature_importances_, i, [0])

        # At checkpoints
        if i % 1000 == 0 and i != 0:
            print("\nWriting importance and performance matrices...\n")
            results[0].index = flatten_index_like(results[0].index)
            results[0].columns = flatten_index_like(results[0].columns)
            results[1].index = gene_labels
            results[0].round(5).to_csv(output + "/imp.csv")
            results[1].round(5).to_csv(output + "/performance.csv")

    # Ensure final state is clean even if no checkpoint triggered
    results[0].index = flatten_index_like(results[0].index)
    results[0].columns = flatten_index_like(results[0].columns)
    results[1].index = gene_labels
    results[0].round(5).to_csv(output + "/imp.csv")
    results[1].round(5).to_csv(output + "/performance.csv")

    return results
